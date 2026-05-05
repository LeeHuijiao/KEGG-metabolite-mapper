---
name: kegg-metabolite-mapper
description: Map differential metabolites to KEGG Compound IDs and KEGG pathways using the official KEGG REST API, with optional ChEBI and PubChem SID identifier conversion, pathway hit summaries, confidence flags, and biological pathway interpretation. Use when the user asks for differential metabolite annotation, KEGG Compound ID lookup, metabolite-to-pathway mapping, pathway biology explanation, metabolomics pathway interpretation, or Chinese metabolomics results such as 差异代谢物、KEGG ID、代谢通路、通路富集、通路生物学解释.
---

# KEGG Metabolite Mapper

Use this skill to turn a differential metabolite list into KEGG Compound IDs, KEGG pathway links, pathway hit summaries, and concise biological interpretation.

Primary source: official KEGG REST API at `https://rest.kegg.jp/`. KEGG states that the REST API is for academic use by academic users and should be limited to 3 calls per second. Always rate-limit requests.

## Default Workflow

1. Inspect the input table.
   - Accept metabolite names, KEGG IDs (`C00031` or `cpd:C00031`), ChEBI IDs, PubChem SIDs, HMDB IDs, m/z, retention time, fold change, p-value, FDR, VIP, and direction if present.
   - If only m/z values are provided, explain that KEGG name matching is not enough; ask for names or formula/adduct information.
2. Run `scripts/map_metabolites.py` for deterministic KEGG mapping when a file or list is provided.
3. Review low-confidence matches manually before interpretation.
4. Summarize pathways by number of matched differential metabolites and direction if available.
5. If the user provides a measured background universe, run pathway enrichment with `--enrichment --background-input <file>`.
   - Use hypergeometric over-representation over the measured/detectable background, not all KEGG compounds.
   - Report p-value, BH-FDR, fold enrichment, overlap count, background pathway size, and overlapping KEGG compounds.
   - For large exploratory backgrounds, use `--skip-pubchem-cid` unless PubChem CID output is needed.
6. If the user provides an organism code such as `hsa`, `mmu`, `rno`, or `ath`, use `--organism <code> --pathway-scope organism` for organism-specific pathway filtering. Use `map` reference pathways by default.
7. Write biological interpretation from pathway context:
   - pathway class and KEGG pathway name,
   - which differential metabolites hit the pathway,
   - up/down direction if available,
   - plausible biological process represented by the pathway,
   - caveats such as ambiguous compound matching, generic pathways, or missing organism specificity.
6. Use complementary skills only when needed:
   - `life-science-research:hmdb-skill` for HMDB-oriented metabolite details,
   - `life-science-research:chebi-skill` for synonym and ontology checks,
   - `life-science-research:pubchem-pug-skill` for PubChem names/properties and CID/SID clarification,
   - `life-science-research:reactome-skill` for non-KEGG pathway context.

## Script

Use the bundled script first for batch work:

```powershell
python "C:\Users\lenovo\.codex\skills\kegg-metabolite-mapper\scripts\map_metabolites.py" `
  --input "metabolites.csv" `
  --out-prefix "output/kegg_metabolites" `
  --name-column "metabolite" `
  --direction-column "direction"
```

For final pathway interpretation, rerun with `--include-pathway-details --max-detail-pathways 20`
to fetch detailed KEGG flat files for the top pathway hits. Do not use detailed fetching for every
pathway in a large exploratory table unless the user explicitly wants slower exhaustive output.

Outputs:

- `<out-prefix>.metabolite_mappings.csv`: one row per metabolite with matched KEGG compound and pathways.
- `<out-prefix>.pathway_summary.csv`: pathway-level hit summary.
- `<out-prefix>.pathway_enrichment.csv`: enrichment results when `--enrichment` is used.
- `<out-prefix>.json`: full mapping details and candidate lists.

If `--name-column` is omitted, the script tries `metabolite`, `name`, `compound`, `Metabolite`, `Compound`, then the first column.

## Confidence Rules

Use these confidence labels:

- `high`: exact KEGG ID supplied, exact name match, or only one strong KEGG candidate.
- `medium`: first KEGG candidate is plausible but synonyms or isomers could conflict.
- `low`: multiple candidates, generic name, partial match, salt/adduct/form ambiguity, or no clear exact hit.
- `unmapped`: no KEGG compound found.

Never silently choose a low-confidence mapping as final. Include it, but flag it for user review.

## Biological Interpretation Template

For each important pathway, write:

```text
Pathway: [KEGG ID] [name]
Matched metabolites: [metabolites with KEGG IDs and direction if available]
Biological meaning: [2-4 sentences explaining the process represented by the pathway and why these metabolites matter]
Interpretation caution: [matching ambiguity, generic map, organism issue, small hit count, or missing direction]
```

Prefer restrained language. A metabolite hit does not prove pathway activation by itself. Use wording such as "suggests altered involvement of", "is consistent with changes in", or "points to possible perturbation of" unless there is quantitative enrichment or independent evidence.

## Official KEGG REST Endpoints Used

- `find/compound/<query>`: name or keyword search for KEGG compounds.
- `conv/compound/pubchem:<sid>` and `conv/compound/chebi:<id>`: convert outside chemical identifiers when provided. KEGG's PubChem conversion uses PubChem SID, not CID.
- `get/<compound>`: retrieve compound names, formula, exact mass, DBLINKS, and pathway cross-references.
- `link/pathway/<compound>`: retrieve pathways linked to a KEGG compound.
- `list/<pathway IDs>` or `get/<pathway IDs>`: retrieve pathway names and flat-file details.

## References

Open these only when needed:

- `references/kegg-rest-notes.md`: endpoint examples, rate limit, and parsing notes.
