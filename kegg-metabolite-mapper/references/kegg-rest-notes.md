# KEGG REST Notes

Official documentation:

- KEGG API landing page: https://www.genome.jp/kegg/rest/
- KEGG API manual: https://www.genome.jp/kegg/rest/keggapi.html

Important policy note:

- KEGG REST API at `rest.kegg.jp` is available for academic use by academic users.
- Limit calls to at most 3 requests per second. Use a delay of at least 0.34 seconds between uncached requests.

Endpoint patterns:

```text
https://rest.kegg.jp/find/compound/glucose
https://rest.kegg.jp/get/C00031
https://rest.kegg.jp/link/pathway/cpd:C00031
https://rest.kegg.jp/list/map00010+map00020
https://rest.kegg.jp/conv/compound/pubchem:<PubChem-SID>
https://rest.kegg.jp/conv/compound/chebi:17234
```

Parsing notes:

- `find`, `list`, `conv`, and `link` return tab-delimited text.
- `get` returns KEGG flat-file text.
- Compound entries may include multiple names separated by semicolons.
- `PATHWAY` fields in compound flat files and `link/pathway/<compound>` are both useful; prefer `link` for batch consistency.
- Reference pathway IDs usually look like `path:map00010`; organism-specific pathway IDs look like `path:hsa00010`. For metabolite-only pathway context, reference `map` pathways are usually acceptable unless the user asks for a specific organism.

Interpretation caveats:

- Name-only metabolite matching can confuse isomers, salts, stereochemistry, conjugates, and generic class names.
- KEGG's PubChem conversion is based on PubChem SID. Do not treat PubChem CID as equivalent to SID.
- m/z-only matching is not reliable without formula, adduct, ion mode, and mass tolerance.
- A pathway hit is not enrichment by itself; pathway enrichment requires a background universe and statistical method.
