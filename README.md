# R_Enrichment

项目代码提取的 `src/R/aucell_scores.R` 提供 `aucell_score_block(expr, gene_sets, auc_max_rank)` 和 `mask_aucell_zero_ranks()`：使用调用者指定的基因集与 cutoff，显式排除稀疏矩阵中的零表达排名，不继承特定项目的生物学标签。须在新 workflow 中显式 source，依赖 AUCell/BiocParallel；掩码检查可独立运行：`Rscript --vanilla tests/test_aucell_zero_ranks.R`。GSVA 默认 CON/EXP 颜色也已统一为 NPG。

跨项目复用的 R 公共分析代码库。仓库管理公共方法、默认配置、环境和项目示例；实际项目在 `projects/<项目名>/` 使用这些方法，数据和运行产物不提交 Git。

本仓库由 `R-Python-scores-enrichment-workflow-template` 按语言拆分；另一语言版本见 [Python_Enrichment](https://github.com/Minkakaku/Python_Enrichment)。保留拆分前提交历史；`legacy/` 和 `source_catalog/` 是原样来源档案，可能包含另一语言，不参与正式执行。

## 目录

```text
src/R/       公共分析方法与配置加载器
config/            公共默认参数与参考配置
envs/              共享环境定义
projects/project1/ 示例项目（未运行）
  config/          项目覆盖、样本表、QC/注释等人工决策
  workflow/        编号短脚本与分阶段 runner
  data/processed/  中间对象
  data/reference/  项目参考资源
  results/         按分析点保存图表与中文说明
  logs/            运行记录
tests/             路径和项目配置检查
```

## 新项目

在仓库根目录执行：

```bash
conda env create -f envs/environment.yaml
cp -R projects/project1 projects/my_project
```

修改 `projects/my_project/config/config.yaml` 和样本/注释表，准备输入对象与参考文件。示例中的物种、reference、root、分组和文件名必须按真实项目核对，不能直接作为生物学结论。

每个编号脚本可单独执行；脚本从自身位置定位项目，因此复制后不会继续写入 project1。Python 包项目先完成 editable install。

```bash
Rscript projects/my_project/workflow/10_enrichment.R
```

可用步骤：
- `10_enrichment.R`
- `11_gsva.R`

`run_all.sh` 对主流程只运行到 QC 检查点；独立分析方法需显式选择步骤，避免把不同候选方法混为同一条流程。已运行脚本不改写，重跑新增编号并使用新输出目录。

## 配置与路径

公共 `config/default.yaml` → `config/reference.yaml` → 项目 `config/config.yaml`，后者覆盖前者。嵌套映射按键合并，列表和标量整体替换。共享方法接收配置/参数，不包含实际项目路径与样本名。项目相对输入输出路径以项目目录为基准；参考文件路径写在配置中。

## 结果与解释

最终图保存同名 PNG/PDF，配套 CSV/TSV 和 `README_results_CN.md` 放在分析点结果目录。中间对象放在 `data/processed/`。细胞或空间 bins 不作为独立生物学重复；条件比较用样本级汇总。方法特有约束见 [AGENTS.md](AGENTS.md) 和公共代码注释。

## 验证

```bash
Rscript tests/test_project.R
```

此检查验证配置、路径或安全导入，不运行真实组学分析。完整分析需要对应生信依赖、真实输入数据和人工 QC/注释。

检查独立入口参数而不运行分析：`Rscript projects/my_project/workflow/10_enrichment.R --dry-run`。
