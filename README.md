# Induced-Gravity Cosmology（EPJC 投稿仓库）

论文 **Induced-Gravity Cosmology: Inflation, Dark Energy, and Dark Matter from a Common Scalar Origin**  
目标期刊：**Eur. Phys. J. C**（Springer Nature `sn-jnl`；SCOAP3 全刊资助，APC=0）

## 目录结构

```text
paper_prd_merged.tex     # 唯一正文（EPJC / sn-jnl）
paper_prd_merged.pdf     # 本地编译产物（pdflatex ×2）
sn-jnl.cls               # Springer Nature 期刊类
sn-mathphys-num.bst      # 物理 numbered 参考文献样式
figures/*.pdf            # 正文插图（fig1–fig3）
scripts/*.py             # 数值验证与作图脚本
README.md
.gitignore
```

## 编译

```powershell
cd D:\work\papers\llun
pdflatex -interaction=nonstopmode paper_prd_merged.tex
pdflatex -interaction=nonstopmode paper_prd_merged.tex
```

依赖：`sn-jnl.cls`、`sn-mathphys-num.bst`、`figures/*.pdf`。  
参考文献已内嵌 `thebibliography`，**无需** bibtex/biber。

## 文档类与格式

- `\documentclass[pdflatex,sn-mathphys-num]{sn-jnl}`
- EPJC：摘要 150–250 词；`Declarations`（Funding / Competing interests / Data / Code / Author contributions）
- 无机构作者登记 city/country（Guiyang, Guizhou, China）

## 数值脚本

环境（本地，勿提交 venv）：

```powershell
python -m venv scripts\.venv
scripts\.venv\Scripts\python.exe -m pip install numpy scipy matplotlib mpmath
$env:PYTHONUNBUFFERED=1
& scripts\.venv\Scripts\python.exe scripts\run_all.py
```

主链（`scripts/run_all.py`）：锁定 $N$ 约定、KG 积分与再加热、$\psi$ Bogoliubov 指数审计、跨跃迁模方程丰度、残余精质、三张图等。

关键结果脚本：

| 脚本 | 内容 |
|---|---|
| `background_and_reheating.py` | 精确 KG 积分、$N$ 窗、再加热通道 |
| `psi_production_bogoliubov.py` | de Sitter 指数 $2\pi$ 审计、$g$ 匹配 |
| `psi_abundance_oscillating.py` | 跨跃迁模方程（幂律谱） |
| `dm_gap_closure_test.py` | 小 $g$ 轻支与自由流长度 |
| `residual_quintessence.py` | 两流体积分、$\Delta w$ |

## 暗物质口径（正文 §V）

真实模方程谱下，丰度匹配落在**轻支**  
$g\simeq1.0\times10^{-7}$，$m_\psi\simeq7.5\times10^{10}$ GeV（冷暗物质）；  
指数闭式重支在模型自身 $T_{\rm reh}$ 下超产。绝对归一化待 lattice/Floquet。

## 投稿打包

上传 zip 建议包含：`paper_prd_merged.tex`、`sn-jnl.cls`、`sn-mathphys-num.bst`、`figures/*.pdf`。  
封面信可写明 SCOAP3 资助（EPJC APC=0）。

**Funding**：当前声明为无外部资助（独立研究者）。若有项目支持请自行改写 `Declarations` 中 Funding 句。

## Git

勿提交：`scripts/.venv/`、`__pycache__/`、LaTeX 中间文件（见 `.gitignore`）。
