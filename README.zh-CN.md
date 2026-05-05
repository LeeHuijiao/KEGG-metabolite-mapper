# KEGG Metabolite Mapper Skill 中文说明

这是一个用于 Codex 的代谢组学 skill，目标是把差异代谢物列表整理成：

```text
差异代谢物 -> KEGG Compound ID -> KEGG pathway -> 通路生物学解释
```

它主要调用 KEGG 官方 REST API，并在需要时通过 PubChem PUG REST 把 KEGG 记录里的 PubChem SID 转成常见 PubChem CID。

## 功能概览

- 根据代谢物名称或 KEGG ID 查询 KEGG Compound ID。
- 查询每个 KEGG Compound 关联的 KEGG pathway。
- 从 KEGG 记录中提取 ChEBI、PubChem 等交叉引用。
- 如果输入表中已有 HMDB / ChEBI / PubChem CID，会保留到最终结果表。
- 输出逐代谢物映射表、通路汇总表、最终展示表和 JSON 详情。
- 支持后续写通路生物学解释，例如能量代谢、氨基酸代谢、脂质代谢等方向的解释。

## 仓库结构

```text
kegg-metabolite-mapper-skill-repo/
├── README.md
├── README.zh-CN.md
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

真正的 Codex skill 是 `kegg-metabolite-mapper/` 这个目录。

## 安装方式

把 `kegg-metabolite-mapper/` 复制到 Codex skills 目录：

```powershell
Copy-Item -Recurse ".\kegg-metabolite-mapper" "C:\Users\lenovo\.codex\skills\kegg-metabolite-mapper"
```

复制完成后，重启 Codex，让新 skill 生效。

## 输入文件格式

最少只需要一列代谢物名称：

```csv
metabolite
Glucose
L-Lactate
Citric acid
```

推荐格式：

```csv
metabolite,direction,log2FC,FDR
Glucose,up,1.23,0.01
L-Lactate,down,-0.86,0.03
Citric acid,up,0.72,0.04
```

推荐字段说明：

| 字段 | 是否必需 | 说明 |
| --- | --- | --- |
| `metabolite` | 必需 | 代谢物名称、KEGG ID、ChEBI ID 或 PubChem SID |
| `direction` | 推荐 | 上调或下调，例如 `up`、`down` |
| `log2FC` | 可选 | 差异倍数 |
| `pvalue` / `FDR` | 可选 | 显著性信息 |
| `HMDB` / `ChEBI` / `PubChem CID` | 可选 | 如果你已有这些 ID，脚本会保留到最终表 |
| `mz` / `rt` | 可选 | 辅助信息；只有 m/z 不能可靠映射 KEGG |

## 运行示例

基础批量映射：

```powershell
python ".\kegg-metabolite-mapper\scripts\map_metabolites.py" `
  --input ".\examples\metabolites_example.csv" `
  --out-prefix ".\output\kegg_metabolites" `
  --name-column "metabolite" `
  --direction-column "direction"
```

如果要为最终通路解释抓取更详细的 KEGG pathway flat-file 信息：

```powershell
python ".\kegg-metabolite-mapper\scripts\map_metabolites.py" `
  --input ".\examples\metabolites_example.csv" `
  --out-prefix ".\output\kegg_metabolites_final" `
  --name-column "metabolite" `
  --direction-column "direction" `
  --include-pathway-details `
  --max-detail-pathways 20
```

## 输出文件

假设 `--out-prefix` 设置为 `output/kegg_metabolites`，脚本会生成：

| 文件 | 内容 |
| --- | --- |
| `output/kegg_metabolites.metabolite_mappings.csv` | 每个代谢物的 KEGG compound 匹配结果、置信度和 pathway 列表 |
| `output/kegg_metabolites.pathway_summary.csv` | 按 pathway 汇总命中数、上下调数量和命中的代谢物 |
| `output/kegg_metabolites.final_table.csv` | 适合论文/汇报展示的最终表 |
| `output/kegg_metabolites.json` | 完整候选匹配、pathway 名称和详细信息 |

`final_table.csv` 的列为：

```text
metabolite | matched name | HMDB | ChEBI | PubChem CID | KEGG ID | KEGG pathway | confidence
```

## 通路生物学解释模板

建议解释时使用谨慎表达，不要把“命中 pathway”直接写成“通路激活”或“通路抑制”。

模板：

```text
Pathway: [KEGG ID] [pathway name]
Matched metabolites: [代谢物名称 + KEGG ID + 上下调方向]
Biological meaning: [2-4 句解释该通路代表的生物过程，以及这些代谢物变化为何重要]
Interpretation caution: [匹配歧义 / 泛通路 / 小命中数 / 缺少背景集 / 缺少物种特异性]
```

中文表述示例：

```text
[通路名] 中检测到 [代谢物A]、[代谢物B] 等差异代谢物，提示该通路相关的 [能量代谢/氨基酸代谢/脂质代谢/氧化应激等] 过程可能受到影响。结合这些代谢物的变化方向，该结果更符合 [简要生物学解释]。需要注意的是，单个或少数代谢物命中 KEGG pathway 不能单独证明通路激活或抑制，仍需结合通路富集、背景集、定量趋势和实验背景进行判断。
```

## 重要注意事项

- 只有 m/z 不足以可靠映射 KEGG，最好提供代谢物名称、分子式、adduct、ion mode 和质量误差范围。
- 名称匹配可能混淆异构体、盐形式、构型、同分异构脂质和泛化合物名称。
- KEGG 的 PubChem 转换接口使用的是 PubChem SID，不是常见的 PubChem CID。
- pathway hit summary 不是严格通路富集分析；富集分析需要背景集和统计方法。
- 代谢物命中某条 pathway 只能提示相关代谢过程可能改变，不能单独证明通路激活或抑制。

## 官方数据源

- KEGG REST API: https://www.genome.jp/kegg/rest/
- KEGG API Manual: https://www.genome.jp/kegg/rest/keggapi.html
- PubChem PUG REST: 用于把 KEGG 记录中的 PubChem SID 转换为 PubChem CID。

使用 KEGG REST API 时请遵守 KEGG 使用政策。本 skill 脚本内置了请求限速，避免超过每秒 3 次请求。

## 在 Codex 中的推荐启动语

```text
请用 kegg-metabolite-mapper 分析我的差异代谢物表。
输入文件：[CSV/TSV 文件路径]
代谢物列名：[metabolite/name/compound/...]
方向列名：[direction，如果没有就写无]
统计列：[log2FC、pvalue、FDR、VIP，如果有]
物种：[human/mouse/plant/bacteria/不限定]
目标：输出 KEGG Compound ID、KEGG pathway 汇总、低置信度匹配清单，并对主要通路写生物学解释。
```
