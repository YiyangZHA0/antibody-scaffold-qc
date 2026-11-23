# 抗体骨架提取与风险检查

从 FASTA、PDB 或 mmCIF 提取抗体可变域，按指定方案划分 FR1–FR4 和 CDR1–CDR3，生成 BoltzGen YAML 与中英文骨架风险报告。支持 VHH、Fab 和显式拆分的 scFv。Fab 输出重链、轻链的可变域（VH、VL）。

## 安装

在 Windows PowerShell 中运行，需已安装 Git 和 Conda。

```powershell
git clone https://github.com/YiyangZHA0/antibody-scaffold-qc.git
cd antibody-scaffold-qc
.\scripts\setup.ps1
conda activate scaffold
scaffold --version
scaffold doctor
```

安装脚本创建或复用 `scaffold` 环境，安装 Windows CPU 依赖并下载 BoltzGen 校验所需的分子词典。已有词典时，可用 `.\scripts\setup.ps1 -SkipMoleculeDownload` 跳过下载，转换时用 `--moldir D:/data/mols.zip` 指定位置。

## 单文件转换

```powershell
scaffold convert tests/data/5y9k.fasta --scheme imgt --out outputs/my_fab
scaffold convert tests/data/5y9k.pdb --scheme chothia
scaffold convert tests/data/5y9k.cif --scheme kabat
```

默认情况下，FASTA 直接生成序列约束 YAML，无需预测三维结构。PDB/CIF 默认保留骨架坐标，也可以用 `--mode sequence` 只提取序列约束。FASTA 要生成结构约束时，使用 `--structure` 提供匹配的结构。指定的输出目录应为新目录或空目录。

| 参数 | 用途 |
|---|---|
| `--scheme` | 选择 `imgt`、`kabat`、`chothia`、`abm`、`maccallum`（别名 `contact`）或 `aho`。AHo 默认采用 IMGT 的 CDR 边界；可用 `--cdr-definition` 另行指定。 |
| `--type` | 选择 `auto`、`fab`、`vhh` 或 `scfv`。 |
| `--heavy-chain`、`--light-chain`、`--vhh-chain` | 指定结构中的输入链。默认使用作者链 ID，可用 `--chain-namespace label` 切换。多对候选链时先运行 `inspect`。 |
| `--records` | 用逗号分隔的记录 ID 选择 FASTA 记录。 |
| `--cdrs`、`--design-chains`、`--cdr-length` | 指定待设计的环、链和替换后的完整环长。支持固定长度和范围。 |

```powershell
scaffold inspect tests/data/5j13.cif --scheme imgt
scaffold convert tests/data/5j13.cif --scheme imgt --type fab --heavy-chain C --light-chain B
scaffold convert tests/data/7eow.pdb --scheme imgt --type vhh --vhh-chain B
scaffold convert tests/data/5y9k.fasta --scheme chothia --design-chains H --cdrs CDR3 --cdr-length H-CDR3=12
scaffold convert tests/data/1x9q.fasta --scheme imgt --type scfv --split-scfv
```

## 混合批量转换

```powershell
scaffold batch tests/data/5y9k.fasta tests/data/5j13.cif tests/data/7eow.pdb --scheme imgt
scaffold batch my_inputs --recursive --scheme chothia --jobs 2 --out outputs/batch_run
```

每个输入产生独立结果目录。需要逐文件指定方案或链时，建立 `jobs.csv`，路径相对于清单所在目录。

```csv
id,input,scheme,type,heavy_chain,light_chain,vhh_chain
fab,tests/data/5j13.cif,chothia,fab,C,B,
vhh,tests/data/7eow.pdb,imgt,vhh,,,B
```

```powershell
scaffold batch --manifest jobs.csv --out outputs/from_manifest
scaffold batch --manifest outputs/from_manifest/failed_jobs.json --out outputs/retry
```

## 如何读风险报告

成功转换后，先看 `framework_check.md`。平台可读 `framework_check.json` 中的 `risk_summary.status`，再按链、FR、`message_zh`、`message_en` 和证据展示具体风险。

| 状态 | 含义 |
|---|---|
| `no_flag` | 基于内置参考库检索及本次可执行的检查，未发现骨架异常。请对照原始构建记录进一步核对。 |
| `review` | 发现需复核的项目。报告列出相关链、FR、位置或参考差异，以及中英文提示。 |
| `insufficient_evidence` | 已完成的检查未发现明确异常，部分位置证据不足，需核对相应端点或原始文件。 |

内置参考包括已知构建体可变域和人源、部分羊驼的 V/J 胚系蛋白。对能核验的原始 FASTA，程序依据 IMGT 对齐位置辅助检查 FR1 的 N 端与 FR4 的 C 端。遇到单 S／双 S 等末端歧义时，需核实具体构建版本。对结构文件，程序还检查声明残基的坐标覆盖和主链原子。风险检查保留输入序列。

## 常用命令与输出

```powershell
scaffold annotate tests/data/5y9k.fasta --scheme imgt --out outputs/annotation
scaffold validate outputs/my_fab
scaffold postcheck outputs/my_fab
scaffold audit tests/data/5y9k.pdb --out outputs/audit_structure
scaffold audit input.fasta --reference intended_construct.fasta --out outputs/audit_reference
scaffold compare outputs/result_a outputs/result_b
scaffold --help
scaffold convert --help
scaffold batch --help
```

- `annotate` 仅导出编号及区域。独立的 `audit` 命令使用 IMGT 检查序列与坐标。已有转换目录可用 `postcheck` 按当时所选方案刷新双语报告。
- `design.yaml` 是 BoltzGen 的设计入口。结构模式下 `scaffold.yaml` 引用同目录的 `scaffold.cif`。编号和残基映射见 `numbering.*` 与 `residue_map.*`。
- `validation.json` 和 `boltzgen_check.log` 记录官方解析和 CLI 检查。YAML 校验通过说明设计约束能被读取；序列完整性请查看骨架风险报告。
- 测试样本及来源位于 `tests/data/`。在已安装的环境运行 `python -m pytest` 可执行开发测试。
