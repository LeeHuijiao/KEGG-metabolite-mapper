# KEGG Metabolite Mapper Skill

A Codex skill for mapping differential metabolites to KEGG Compound IDs and KEGG pathways using the official KEGG REST API, then producing pathway hit summaries and biological interpretation-ready tables.

## What It Does

- Maps metabolite names or KEGG IDs to KEGG Compound IDs.
- Retrieves KEGG pathways linked to each compound.
- Extracts available ChEBI and PubChem cross-references from KEGG records.
- Converts KEGG PubChem SIDs to PubChem CIDs through PubChem PUG REST when possible.
- Preserves user-provided HMDB / ChEBI / PubChem CID columns.
- Outputs:
  - metabolite-level mapping table,
  - pathway-level summary table,
  - final report-style table,
  - JSON with candidates and details.

## Repository Layout

```text
kegg-metabolite-mapper-skill-repo/
├── README.md
├── LICENSE
├── .gitignore
├── examples/
│   └── metabolites_example.csv
└── kegg-metabolite-mapper/
    ├── SKILL.md
    ├── agents/
    │   └── openai.yaml
    ├── references/
    │   └── kegg-rest-notes.md
    └── scripts/
        └── map_metabolites.py
```

The actual Codex skill is the `kegg-metabolite-mapper/` folder.

## Install

Copy `kegg-metabolite-mapper/` into your Codex skills directory:

```powershell
Copy-Item -Recurse ".\kegg-metabolite-mapper" "C:\Users\lenovo\.codex\skills\kegg-metabolite-mapper"
```

Restart Codex after installation.

## Quick Script Usage

```powershell
python ".\kegg-metabolite-mapper\scripts\map_metabolites.py" `
  --input ".\examples\metabolites_example.csv" `
  --out-prefix ".\output\kegg_metabolites" `
  --name-column "metabolite" `
  --direction-column "direction"
```

For final interpretation of top pathways, fetch additional KEGG flat-file details:

```powershell
python ".\kegg-metabolite-mapper\scripts\map_metabolites.py" `
  --input ".\examples\metabolites_example.csv" `
  --out-prefix ".\output\kegg_metabolites_final" `
  --name-column "metabolite" `
  --direction-column "direction" `
  --include-pathway-details `
  --max-detail-pathways 20
```

## Input Format

Minimum:

```csv
metabolite
Glucose
L-Lactate
Citric acid
```

Recommended:

```csv
metabolite,direction,log2FC,FDR
Glucose,up,1.23,0.01
L-Lactate,down,-0.86,0.03
Citric acid,up,0.72,0.04
```

Optional identifier columns such as `HMDB`, `ChEBI`, and `PubChem CID` are preserved in the final output table if present.

## Outputs

Given `--out-prefix output/kegg_metabolites`, the script writes:

- `output/kegg_metabolites.metabolite_mappings.csv`
- `output/kegg_metabolites.pathway_summary.csv`
- `output/kegg_metabolites.final_table.csv`
- `output/kegg_metabolites.json`

`final_table.csv` uses this report-friendly schema:

```text
metabolite | matched name | HMDB | ChEBI | PubChem CID | KEGG ID | KEGG pathway | confidence
```

## Data Sources and Policy

- KEGG REST API: https://www.genome.jp/kegg/rest/
- KEGG API manual: https://www.genome.jp/kegg/rest/keggapi.html
- PubChem PUG REST is used only to convert KEGG-provided PubChem SIDs to PubChem CIDs.

KEGG REST should be used according to KEGG's usage policy. This skill rate-limits KEGG calls to avoid exceeding 3 requests per second.

## Important Caveats

- Name-only matching can confuse isomers, salts, stereochemistry, conjugates, and generic class names.
- m/z-only input is not enough for reliable KEGG mapping without formula, adduct, ion mode, and mass tolerance.
- KEGG's PubChem conversion uses PubChem SID, not PubChem CID.
- A pathway hit is not statistical enrichment by itself; enrichment requires a background universe and statistical method.
- A metabolite hit does not prove pathway activation or inhibition without additional evidence.

