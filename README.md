# MutationsCounter

用于对 mpileup 结果（`.mpileup.cns.filter.xls`，制表符分隔）进行突变/缺失统计、位置归一化分析并生成图表与 HTML 报告的脚本集合。

## 功能一览
- 统计总突变率、目标碱基突变率（toA/toT/toC/toG）、缺失率等指标
- 统计整体 substitution 突变谱（Template->Alt 的 4x4 matrix，off-diagonal 总和为 100%）并注明总 substitution 突变率
- 统计整体 deletion（总 deletion 率 + A/T/C/G 各自的 deletion 率）
- 输出每个 position 上、按 Template=A/T/C/G 分组的 deletion 率（不做长度归一化）
- 按序列长度归一化到 1-100 的位置区间，便于跨序列比较
- 生成多张 PNG 图表与 CSV 统计文件
- 自动生成 HTML 报告（可在浏览器中打印为 PDF）
- 支持单文件与批量处理（Windows 还提供交互式 `.bat`）

## 依赖与环境
- Python 3.x
- 依赖库：`pandas`、`numpy`、`matplotlib`
- GUI 目录选择需要 `tkinter`（Windows 通常自带；Linux 可能需要单独安装）

安装依赖：

```bash
pip install pandas numpy matplotlib
```

## 输入文件格式
脚本读取的是“制表符分隔”的文本文件（扩展名为 `.xls`，但**不是** Excel 二进制格式）。

文件名后缀历史上常见两类，脚本批处理会同时兼容：
- `*.mpileup.cns.filter.xls`
- `*.mpileup*.cns*.xls`（例如 `*.mpileup.1Dseq.cns.xls`）

必须至少包含以下列（区分大小写）：

- `Chrom`：序列/染色体名称
- `Pos`：位置（数字）
- `Template`：模板碱基（A/T/C/G）
- `Depths`：测序深度
- `Muts`：突变数（不含缺失）
- `Acount`/`Tcount`/`Ccount`/`Gcount`：突变为对应碱基的数量
- `delcount`：缺失数量

示例表头：

```text
Chrom	Pos	Template	Depths	Muts	Acount	Tcount	Ccount	Gcount	delcount
```

## 使用方法

### 1) 单文件分析（命令行）

```bash
python analyze_mutations_single.py path/to/sample.mpileup.cns.filter.xls
```

- 如果文件名不包含 `.mpileup.cns`，脚本会提示是否继续。

### 2) 批量分析（选择目录）

```bash
python analyze_mutations.py
```

- 会弹出目录选择窗口，自动处理目录内所有 `*.mpileup.cns.filter.xls` 文件。

### 3) Windows 交互式入口（推荐）

```bat
run_mutation_analysis.bat
```

支持三种模式：
1. 选择目录批处理
2. 当前目录批处理
3. 手动指定单个文件（支持拖拽文件到窗口）

## 输出内容
输出目录位于**输入文件所在目录**，结构如下：

```text
<输入文件名>_analysis/
  images/   # PNG 图表
  data/     # CSV 统计数据
  reports/  # HTML 报告
```

主要文件：
- `data/normalized_position_mutation_rates.csv`：归一化位置（1-100）的统计结果
- `reports/<文件名>_report.html`：可直接打开的报告页面
- `images/` 内包含如下图表（部分示例）：
  - `normalized_position_mutation_rate.png`
  - `normalized_position_total_rate.png`
  - `normalized_position_to_base_specific_rates.png`
  - `normalized_position_from_base_specific_rates.png`
  - `normalized_position_base_mut_vs_del_comparison.png`
  - `normalized_position_base_specific_deletion_rates.png`
  - `normalized_position_deletion_rate.png`
  - `normalized_position_error_types_comparison.png`

## 生成 PDF 报告
HTML 报告不会自动打开。你可以：
1. 用浏览器打开 `reports/<文件名>_report.html`
2. 点击页面顶部“打印报告为PDF”按钮或使用浏览器打印功能导出 PDF

## 常见问题

- **缺少依赖报错**：请先安装 `pandas`/`numpy`/`matplotlib`。
- **中文乱码或中文不显示**：脚本默认使用 `SimHei` 字体，非 Windows 环境可能缺失。可在脚本顶部改为你系统可用字体。
- **提示未找到文件**：确认输入路径正确，且为制表符分隔文件。
- **批处理未找到文件**：仅匹配 `*.mpileup.cns.filter.xls`，请确认文件名后缀。

## 项目结构
- `analyze_mutations.py`：核心分析逻辑与批量处理入口（带目录选择）
- `analyze_mutations_single.py`：单文件命令行入口
- `run_mutation_analysis.bat`：Windows 交互式入口

## 许可
当前仓库未包含 `LICENSE` 文件。如需发布或共享，请补充相应许可证说明。
