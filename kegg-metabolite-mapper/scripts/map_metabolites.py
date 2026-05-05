#!/usr/bin/env python3
"""Map metabolite names or IDs to KEGG compounds and pathways.

Uses the official KEGG REST API and only Python standard-library modules.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


BASE = "https://rest.kegg.jp"
DELAY_SECONDS = 0.36


class KeggClient:
    def __init__(self, delay: float = DELAY_SECONDS):
        self.delay = delay
        self.cache: Dict[str, str] = {}
        self.last_request = 0.0

    def get_text(self, endpoint: str) -> str:
        endpoint = endpoint.lstrip("/")
        url = f"{BASE}/{endpoint}"
        if url in self.cache:
            return self.cache[url]
        elapsed = time.time() - self.last_request
        if elapsed < self.delay:
            time.sleep(self.delay - elapsed)
        req = urllib.request.Request(url, headers={"User-Agent": "Codex KEGG metabolite mapper"})
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                text = resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                text = ""
            else:
                raise
        self.last_request = time.time()
        self.cache[url] = text
        return text


class PubChemClient:
    def __init__(self, delay: float = 0.2):
        self.delay = delay
        self.cache: Dict[str, str] = {}
        self.last_request = 0.0

    def get_text(self, url: str) -> str:
        if url in self.cache:
            return self.cache[url]
        elapsed = time.time() - self.last_request
        if elapsed < self.delay:
            time.sleep(self.delay - elapsed)
        req = urllib.request.Request(url, headers={"User-Agent": "Codex KEGG metabolite mapper"})
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                text = resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError:
            text = ""
        self.last_request = time.time()
        self.cache[url] = text
        return text

    def sid_to_cid(self, sid: str) -> str:
        sid = sid.strip()
        if not sid:
            return ""
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{urllib.parse.quote(sid)}/cids/TXT"
        text = self.get_text(url)
        cids = [line.strip() for line in text.splitlines() if line.strip().isdigit()]
        return ";".join(cids)


def read_table(path: Path, name_column: str | None) -> Tuple[List[dict], str]:
    if not path.exists():
        raise FileNotFoundError(path)
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    sample = text[:4096]
    if "\t" in sample and sample.count("\t") >= sample.count(","):
        dialect = csv.excel_tab
    else:
        dialect = csv.excel
    rows = list(csv.DictReader(text.splitlines(), dialect=dialect))
    if rows:
        columns = list(rows[0].keys())
        if name_column:
            if name_column not in columns:
                raise ValueError(f"Name column '{name_column}' not found. Available columns: {columns}")
            return rows, name_column
        for candidate in ("metabolite", "name", "compound", "Metabolite", "Name", "Compound"):
            if candidate in columns:
                return rows, candidate
        return rows, columns[0]

    # Fallback for one-name-per-line text files.
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    return [{"metabolite": line} for line in lines], "metabolite"


def normalize_name(value: str) -> str:
    value = (value or "").strip()
    value = re.sub(r"\s+", " ", value)
    value = re.sub(r"\s*\([^)]*(?:internal standard|is|qc|rt|mz)[^)]*\)\s*", " ", value, flags=re.I)
    return value.strip(" ;,")


def parse_tsv(text: str) -> List[Tuple[str, str]]:
    rows = []
    for line in text.splitlines():
        if not line.strip() or "\t" not in line:
            continue
        left, right = line.split("\t", 1)
        rows.append((left.strip(), right.strip()))
    return rows


def parse_kegg_flat(text: str) -> Dict[str, object]:
    data: Dict[str, object] = {"raw": text}
    current = None
    for line in text.splitlines():
        if line.startswith("///"):
            break
        key = line[:12].strip()
        value = line[12:].strip()
        if key:
            current = key
            if key in data:
                if isinstance(data[key], list):
                    data[key].append(value)
                else:
                    data[key] = [data[key], value]
            else:
                data[key] = value
        elif current and value:
            if current in data:
                if isinstance(data[current], list):
                    data[current].append(value)
                else:
                    data[current] = f"{data[current]} {value}"
    return data


def parse_dblinks(raw_text: str) -> Dict[str, List[str]]:
    links: Dict[str, List[str]] = defaultdict(list)
    in_dblinks = False
    for line in raw_text.splitlines():
        key = line[:12].strip()
        value = line[12:].strip()
        if key == "DBLINKS":
            in_dblinks = True
        elif key and in_dblinks:
            break
        if not in_dblinks or not value:
            continue
        match = re.match(r"([^:]+):\s*(.+)$", value)
        if not match:
            continue
        db = match.group(1).strip()
        vals = [v.strip() for v in re.split(r"\s+", match.group(2).strip()) if v.strip()]
        links[db].extend(vals)
    return {db: sorted(set(vals)) for db, vals in links.items()}


def first_existing(row: dict, names: Iterable[str]) -> str:
    lowered = {str(k).lower(): k for k in row.keys()}
    for name in names:
        key = lowered.get(name.lower())
        if key is not None and str(row.get(key, "")).strip():
            return str(row.get(key, "")).strip()
    return ""


def names_from_find_description(desc: str) -> List[str]:
    desc = desc.split("\t")[-1]
    return [part.strip() for part in desc.split(";") if part.strip()]


def score_candidate(query: str, candidate_id: str, desc: str) -> Tuple[int, str]:
    q = query.lower().strip()
    names = [n.lower() for n in names_from_find_description(desc)]
    if candidate_id.lower().endswith(q) or q in (candidate_id.lower(), candidate_id.lower().replace("cpd:", "")):
        return 100, "exact_id"
    if any(q == n for n in names):
        return 95, "exact_name"
    if any(q.replace("-", " ") == n.replace("-", " ") for n in names):
        return 90, "normalized_exact_name"
    if any(q in n or n in q for n in names):
        return 70, "partial_name"
    tokens = [tok for tok in re.split(r"[^a-z0-9]+", q) if len(tok) > 2]
    if tokens and any(all(tok in n for tok in tokens) for n in names):
        return 60, "token_match"
    return 40, "weak"


def confidence_from_score(score: int, candidate_count: int) -> str:
    if score >= 95 and candidate_count <= 3:
        return "high"
    if score >= 80:
        return "medium"
    if score >= 60:
        return "low"
    return "low"


def clean_cpd_id(value: str) -> str:
    value = value.strip()
    if value.startswith("compound:"):
        value = value.replace("compound:", "cpd:", 1)
    if re.fullmatch(r"C\d{5}", value, flags=re.I):
        return "cpd:" + value.upper()
    if re.fullmatch(r"cpd:C\d{5}", value, flags=re.I):
        return "cpd:" + value.split(":", 1)[1].upper()
    return value


def map_identifier(client: KeggClient, query: str) -> List[Tuple[str, str, int, str]]:
    q = query.strip()
    cpd = clean_cpd_id(q)
    if re.fullmatch(r"cpd:C\d{5}", cpd):
        return [(cpd, "KEGG ID supplied", 100, "exact_id")]

    chebi = re.search(r"(?:CHEBI:|chebi:)?(\d{2,8})", q) if "chebi" in q.lower() else None
    if chebi:
        text = client.get_text(f"conv/compound/chebi:{chebi.group(1)}")
        rows = parse_tsv(text)
        if rows:
            return [(clean_cpd_id(right), f"Converted from ChEBI:{chebi.group(1)}", 95, "chebi_conv") for _, right in rows]

    pubchem = re.search(r"(?:pubchem:|sid:|SID:)\s*(\d+)", q)
    if pubchem:
        text = client.get_text(f"conv/compound/pubchem:{pubchem.group(1)}")
        rows = parse_tsv(text)
        if rows:
            return [(clean_cpd_id(right), f"Converted from PubChem SID:{pubchem.group(1)}", 95, "pubchem_sid_conv") for _, right in rows]

    encoded = urllib.parse.quote(q, safe="")
    rows = parse_tsv(client.get_text(f"find/compound/{encoded}"))
    candidates = []
    for left, right in rows:
        score, reason = score_candidate(q, left, right)
        candidates.append((clean_cpd_id(left), right, score, reason))
    candidates.sort(key=lambda item: item[2], reverse=True)
    return candidates


def get_compound(client: KeggClient, cpd_id: str) -> Dict[str, object]:
    kid = cpd_id.split(":", 1)[-1]
    return parse_kegg_flat(client.get_text(f"get/{kid}"))


def get_pathways_for_compound(client: KeggClient, cpd_id: str) -> List[str]:
    text = client.get_text(f"link/pathway/{cpd_id}")
    pathways = []
    for _, right in parse_tsv(text):
        if right.startswith("path:"):
            pathways.append(right.replace("path:", ""))
    return sorted(set(pathways))


def list_pathway_names(client: KeggClient, pathway_ids: Iterable[str]) -> Dict[str, str]:
    ids = sorted(set(pid.replace("path:", "") for pid in pathway_ids if pid))
    out = {}
    for i in range(0, len(ids), 10):
        chunk = ids[i:i + 10]
        text = client.get_text("list/" + "+".join(chunk))
        for left, right in parse_tsv(text):
            out[left.replace("path:", "")] = right
    return out


def get_pathway_details(client: KeggClient, pathway_id: str) -> Dict[str, object]:
    pid = pathway_id.replace("path:", "")
    return parse_kegg_flat(client.get_text(f"get/{pid}"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Map metabolites to KEGG compounds and pathways.")
    parser.add_argument("--input", required=True, help="CSV/TSV/TXT input file.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix, e.g. output/kegg_results.")
    parser.add_argument("--name-column", default=None, help="Column containing metabolite names or IDs.")
    parser.add_argument("--direction-column", default=None, help="Optional up/down direction column.")
    parser.add_argument("--top-n", type=int, default=5, help="Number of candidates to keep in JSON.")
    parser.add_argument("--include-pathway-details", action="store_true", help="Also fetch KEGG pathway flat-file details. Slower; use for final interpretation.")
    parser.add_argument("--max-detail-pathways", type=int, default=30, help="Maximum pathways for --include-pathway-details.")
    parser.add_argument("--skip-pubchem-cid", action="store_true", help="Skip PubChem SID-to-CID lookup.")
    args = parser.parse_args()

    in_path = Path(args.input)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    rows, name_col = read_table(in_path, args.name_column)
    client = KeggClient()
    pubchem = PubChemClient()

    mappings = []
    pathway_hits = defaultdict(list)
    all_pathways = set()

    for idx, row in enumerate(rows, start=1):
        original = row.get(name_col, "")
        query = normalize_name(original)
        direction = row.get(args.direction_column, "") if args.direction_column else ""
        result = {
            "input_index": idx,
            "input_name": original,
            "query": query,
            "direction": direction,
            "matched_kegg_compound": "",
            "matched_name": "",
            "confidence": "unmapped",
            "match_reason": "",
            "pathways": [],
            "external_ids": {"HMDB": "", "ChEBI": "", "PubChem_SID": "", "PubChem_CID": ""},
            "candidates": [],
            "input_row": row,
        }
        if not query:
            mappings.append(result)
            continue
        try:
            candidates = map_identifier(client, query)
            result["candidates"] = [
                {"kegg_compound": cpd, "description": desc, "score": score, "reason": reason}
                for cpd, desc, score, reason in candidates[: args.top_n]
            ]
            if candidates:
                cpd, desc, score, reason = candidates[0]
                result["matched_kegg_compound"] = cpd
                result["matched_name"] = "; ".join(names_from_find_description(desc)) or desc
                result["confidence"] = confidence_from_score(score, len(candidates))
                result["match_reason"] = reason
                compound = get_compound(client, cpd)
                if compound.get("NAME"):
                    result["matched_name"] = str(compound["NAME"]).rstrip(";")
                dblinks = parse_dblinks(str(compound.get("raw", "")))
                hmdb = first_existing(row, ("HMDB", "hmdb", "HMDB ID", "hmdb_id"))
                chebi = first_existing(row, ("ChEBI", "chebi", "ChEBI ID", "chebi_id"))
                pubchem_cid = first_existing(row, ("PubChem CID", "pubchem_cid", "CID", "cid"))
                pubchem_sid = first_existing(row, ("PubChem SID", "pubchem_sid", "SID", "sid"))
                if not hmdb:
                    hmdb = ";".join(dblinks.get("HMDB", []))
                if not chebi:
                    chebi = ";".join(dblinks.get("ChEBI", []))
                if not pubchem_sid:
                    pubchem_sid = ";".join(dblinks.get("PubChem", []))
                if not pubchem_cid and pubchem_sid and not args.skip_pubchem_cid:
                    cid_values = []
                    for sid in pubchem_sid.split(";"):
                        cid = pubchem.sid_to_cid(sid)
                        if cid:
                            cid_values.extend(cid.split(";"))
                    pubchem_cid = ";".join(sorted(set(cid_values)))
                result["external_ids"] = {
                    "HMDB": hmdb,
                    "ChEBI": chebi,
                    "PubChem_SID": pubchem_sid,
                    "PubChem_CID": pubchem_cid,
                }
                pathways = get_pathways_for_compound(client, cpd)
                result["pathways"] = pathways
                all_pathways.update(pathways)
                for pid in pathways:
                    pathway_hits[pid].append(result)
        except Exception as exc:
            result["confidence"] = "error"
            result["match_reason"] = str(exc)
        mappings.append(result)

    pathway_names = list_pathway_names(client, all_pathways)
    pathway_details = {}
    if args.include_pathway_details:
        ranked_pids = [pid for pid, _ in sorted(pathway_hits.items(), key=lambda kv: (-len(kv[1]), kv[0]))]
        for pid in ranked_pids[: args.max_detail_pathways]:
            try:
                pathway_details[pid] = get_pathway_details(client, pid)
            except Exception:
                pathway_details[pid] = {}

    mapping_csv = out_prefix.with_suffix(".metabolite_mappings.csv")
    with mapping_csv.open("w", newline="", encoding="utf-8-sig") as fh:
        fieldnames = [
            "input_name", "query", "matched_kegg_compound", "matched_name", "confidence",
            "match_reason", "direction", "kegg_pathway_ids", "kegg_pathway_names"
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for item in mappings:
            pids = item["pathways"]
            writer.writerow({
                "input_name": item["input_name"],
                "query": item["query"],
                "matched_kegg_compound": item["matched_kegg_compound"],
                "matched_name": item["matched_name"],
                "confidence": item["confidence"],
                "match_reason": item["match_reason"],
                "direction": item["direction"],
                "kegg_pathway_ids": ";".join(pids),
                "kegg_pathway_names": "; ".join(pathway_names.get(pid, pid) for pid in pids),
            })

    summary_csv = out_prefix.with_suffix(".pathway_summary.csv")
    with summary_csv.open("w", newline="", encoding="utf-8-sig") as fh:
        fieldnames = ["pathway_id", "pathway_name", "pathway_class", "hit_count", "up_count", "down_count", "matched_metabolites"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for pid, hits in sorted(pathway_hits.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            dirs = Counter((h.get("direction") or "").lower() for h in hits)
            details = pathway_details.get(pid, {})
            writer.writerow({
                "pathway_id": pid,
                "pathway_name": pathway_names.get(pid, pid),
                "pathway_class": details.get("CLASS", ""),
                "hit_count": len(hits),
                "up_count": dirs.get("up", 0) + dirs.get("upregulated", 0),
                "down_count": dirs.get("down", 0) + dirs.get("downregulated", 0),
                "matched_metabolites": "; ".join(f"{h['input_name']} ({h['matched_kegg_compound']})" for h in hits),
            })

    final_csv = out_prefix.with_suffix(".final_table.csv")
    with final_csv.open("w", newline="", encoding="utf-8-sig") as fh:
        fieldnames = [
            "metabolite", "matched name", "HMDB", "ChEBI", "PubChem CID",
            "KEGG ID", "KEGG pathway", "confidence"
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for item in mappings:
            ids = item.get("external_ids", {})
            pids = item["pathways"]
            pathway_text = "; ".join(f"{pid} {pathway_names.get(pid, pid)}" for pid in pids)
            writer.writerow({
                "metabolite": item["input_name"],
                "matched name": item["matched_name"],
                "HMDB": ids.get("HMDB", ""),
                "ChEBI": ids.get("ChEBI", ""),
                "PubChem CID": ids.get("PubChem_CID", ""),
                "KEGG ID": item["matched_kegg_compound"],
                "KEGG pathway": pathway_text,
                "confidence": item["confidence"],
            })

    json_path = out_prefix.with_suffix(".json")
    payload = {
        "input": str(in_path),
        "name_column": name_col,
        "mapping_csv": str(mapping_csv),
        "pathway_summary_csv": str(summary_csv),
        "final_table_csv": str(final_csv),
        "mappings": mappings,
        "pathway_names": pathway_names,
        "pathway_details": pathway_details,
    }
    json_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"Wrote {mapping_csv}")
    print(f"Wrote {summary_csv}")
    print(f"Wrote {final_csv}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
