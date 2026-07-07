# torc1predict — 植物 TORC1 结构预测流程

基于 **Boltz-2** 对 *拟南芥* TORC1 复合物（TOR + RAPTOR1B + LST8-1）进行 **AI 驱动的从头结构预测**。

## 概述

本流程实现了 Protocol 1，在没有任何实验模板、Cryo-EM 密度图或先验结构模型的情况下预测植物 TORC1 复合物的三维结构。输入仅包含氨基酸序列；AI 模型从 MSA 中的共进化信号和学习的蛋白质结构知识推断复合物结构、链间界面和置信度指标。

### 目标蛋白

| 链 | 蛋白 | 基因 ID | UniProt | 长度 (aa) | 结构域 |
|-------|---------|---------|---------|-------------|---------|
| A | TOR | AT1G50030 | Q9FR53 | ~2,481 | HEAT repeats, FAT, FRB, Kinase, FATC |
| B | RAPTOR1B | AT3G08850 | Q93ZJ0 | ~1,326 | RNC, HEAT repeats, WD40 |
| C | LST8-1 | AT3G18140 | Q9LSM2 | ~333 | 7× WD40 β-propeller |

### 流程阶段

```
Stage 1 → Stage 2 → Stage 3 → Stage 4 → Stage 5
  │          │         │         │         │
  │          │         │         │         └── ConSurf + FoldX + MD + 界面分析
  │          │         │         └── pLDDT, PAE, 界面, 氢键
  │          │         └── Boltz-2 预测
  │          └── 多链 YAML/FASTA 输入构建
  └── 序列获取与验证 (TAIR/UniProt)
```

## 安装

```bash
# 克隆并进入项目
cd torc1predict

# 创建 conda 环境（推荐）
conda create -n torc1predict python=3.10
conda activate torc1predict

# 安装依赖
pip install -r requirements.txt

# 安装结构预测引擎：
# Boltz-2（推荐使用 CUDA 加速）：
pip install boltz[cuda] -U
```

> **建议**：PyMOL/ChimeraX 只用于打开流程生成的 `.pml` 可视化脚本，不是预测和分析的硬依赖。为避免 `pymol-open-source` 与预测环境中的 NumPy 版本冲突，主环境不要安装 PyMOL；需要可视化时，使用外部 PyMOL/ChimeraX，或另建单独环境安装 `requirements-visualization.txt`。

### 外部依赖（可选）

| 工具 | 用途 | 安装方式 |
|------|---------|-------------|
| MMseqs2 | MSA 生成 | `conda install -c bioconda mmseqs2` |
| IUPred3 | 无序区域预测 | Web API（自动）或 `pip install iupred3` |
| ConSurf | 保守性分析 | https://consurf.tau.ac.il |
| FoldX | ΔΔG 预测 | http://foldxsuite.crg.eu（需学术许可） |
| OpenMM | 分子动力学模拟 | `conda install -c conda-forge openmm` |
| PyMOL / ChimeraX | 可视化 | 外部安装，或单独环境中 `pip install -r requirements-visualization.txt` |
| FreeSASA | 溶剂可及表面积计算 | `pip install freesasa` |

## 快速开始

```bash
# 1. 编辑 Boltz 配置
vim boltz/config/config.yaml

# 2. 运行完整流程
python pipeline.py

# 或使用 Boltz 专用入口
bash boltz/scripts/run_torc1.sh

# 3. 逐步运行（首次运行推荐）：
python pipeline.py --stage 1        # 获取并验证序列
python pipeline.py --stage 2        # 构建多链输入
python pipeline.py --stage 3        # 运行预测（耗时数小时）
python pipeline.py --stage 4        # 分析结果
python pipeline.py --stage 5        # 下游分析

# 4. 干运行（仅验证配置，不执行实际计算）：
python pipeline.py --dry-run
```

## 配置

编辑 [`boltz/config/config.yaml`](boltz/config/config.yaml) 设置以下内容：

- **目标序列** — 基因 ID、UniProt 登录号
- **复合物化学计量** — 链标签及映射关系
- **序列来源** — `uniprot`（默认）或 `tair`
- **预测参数** — 采样数、MSA 模式、回收次数等
- **分析阈值** — pLDDT、PAE、距离截断值
- **MD 参数** — 力场、温度、模拟时长
- **亲和力预测** — 结合物概率阈值、亲和力截断值

> **AlphaFold3 支持**：当前暂未部署 AF3。如需使用，请参考 `af3/config/config.yaml` 和 `af3/docs/ALPHAFOLD3_TORC1_CN.md`。

## 目录分类

项目首先按预测引擎分组，避免混淆输入、输出和中间文件：

- `boltz/` — Boltz-2 的配置、脚本、文档、输入、输出和日志。
- `af3/` — AlphaFold 3 的配置、脚本、文档、输入、输出和日志。
- `modules/`、`pipeline.py` — 两个引擎共享的流程代码和通用入口。

推荐结构：

```
boltz/
├── config/config.yaml
├── scripts/
├── docs/
├── input_data/
├── output/
└── logs/

af3/
├── config/config.yaml
├── scripts/
├── docs/
├── input_data/
├── output/
└── logs/
```

## 小分子结合亲和力预测（Boltz-2）

Boltz-2 可以同时预测复合物结构和**小分子-蛋白结合亲和力**，适用于药物发现流程。官方 affinity 模块要求 binder 是小分子 ligand；它不应用于 TOR/RAPTOR/LST8 这类纯蛋白-蛋白复合物界面亲和力。对当前 TORC1 预测，应主要使用 `protein_iptm`、`pair_chains_iptm`、PAE、pLDDT 和界面接触分析判断蛋白-蛋白界面可信度。

如果后续加入小分子 ligand，并在 YAML 的 `properties.affinity` 中指定 ligand chain，流程会输出两个亲和力指标：

| 字段 | 描述 | 取值范围 | 使用场景 |
|-------|-------------|-------|----------|
| `affinity_probability_binary` | 小分子配体为结合物的预测概率 | 0–1 | 命中发现（从诱饵中检测结合物） |
| `affinity_pred_value` | 结合亲和力，以 log₁₀(IC₅₀) 表示，由 μM 单位的 IC₅₀ 导出 | ~ -3 到 +3 | 配体优化（命中到先导、先导优化） |

这两个指标基于不同的大规模数据集训练，使用不同的监督方式，应在不同的场景中使用。详见 [Boltz-2 预测说明](boltz/docs/prediction.md)。

## 输出结构

```
boltz/output/
├── models/                          # PDB/mmCIF 结构（按置信度排序）
│   ├── seed_1.pdb
│   ├── seed_2.pdb
│   └── ...
├── rankings/                        # 模型置信度排名
│   └── ranking_scores.json
├── pae/                             # 预测对齐误差矩阵
│   └── pae.npy
├── plddt/                           # 逐残基置信度值
│   └── plddt.npy
├── affinity/                        # 亲和力预测（仅 Boltz-2）
│   ├── affinity_results.json        # affinity_pred_value + affinity_probability_binary
│   └── affinity_table.csv           # 表格格式（如有）
├── plots/                           # 可视化
│   ├── seed_1_plddt.png             # 逐残基 pLDDT 图
│   ├── seed_1_pae.png               # PAE 热力图
│   ├── seed_1_view.pml              # PyMOL 会话脚本
│   └── seed_1_analysis.json         # 完整分析结果
├── downstream/                      # 下游分析
│   ├── consurf/                     # ConSurf 保守性
│   ├── foldx/                       # FoldX 突变扫描
│   ├── md/                          # MD 轨迹 + 日志
│   └── interface_profile/           # 界面性质分析
└── pipeline_summary.json            # 汇总摘要
```

流程日志单独写入 `boltz/logs/`。

## 关键指标

本流程使用以下指标评估模型质量：

| 指标 | 描述 | 优良值 |
|--------|-------------|------------|
| **pLDDT** | 逐残基置信度（来自 B 因子） | ≥ 70（可信），≥ 90（非常高） |
| **PAE** | 残基对之间的预测对齐误差 | ≤ 5 Å（链间预测可信） |
| **ipTM** | 界面预测 TM-score | ≥ 0.7（界面预测可信） |
| **pTM** | 整体预测 TM-score | ≥ 0.7（折叠良好） |
| **界面接触** | 链间残基对距离 ≤ 5 Å | 每个界面 ≥ 10 对 |
| **氢键** | 链间氢键 | 每个界面多个 |
| **ΔΔG (FoldX)** | 突变对稳定性的影响 | > 2 kcal/mol = 去稳定化 |

## 模块参考

| 模块 | 文件 | 功能 |
|--------|------|----------|
| 序列准备 | [`modules/input_prep.py`](modules/input_prep.py) | `prepare_sequences()` |
| 输入构建 | [`modules/input_builder.py`](modules/input_builder.py) | `build_inputs()` |
| 结构预测 | [`modules/prediction.py`](modules/prediction.py) | `run_prediction()` |
| 结果分析 | [`modules/analysis.py`](modules/analysis.py) | `analyze_structure()` |
| 下游分析 | [`modules/downstream.py`](modules/downstream.py) | `run_downstream_analyses()` |

每个模块也可独立运行进行测试：

```bash
python modules/input_prep.py
python modules/input_builder.py
python modules/prediction.py
python modules/analysis.py
python modules/downstream.py
```

## 核心设计原则

1. **从头预测** — 不使用模板、实验密度图或先验 PDB 结构。AI 模型完全从序列 + MSA 推断一切。

2. **避免模板偏差** — 仅以序列作为输入，流程能够发现可能不同于酵母/人类 TORC1 结构的植物特有结构特征。

3. **严格验证** — 多个正交指标（pLDDT、PAE、ipTM、界面接触、氢键、FoldX ΔΔG、MD 稳定性）为模型质量提供汇聚证据。

4. **可复现性** — 所有参数集中在一个配置文件中。流程记录每一步日志。随机种子是确定性的。

## 局限性

- **预测精度**：AI 模型对大型复合物（总长 >3,000 aa）的置信度可能较低。TORC1 复合物约 4,140 aa，已达当前方法的极限。
- **MSA 深度**：植物特有序列的 MSA 可能比酵母/人类的更浅，影响植物特有区域的预测质量。
- **构象动态**：单一静态模型无法捕捉 TORC1 的构象灵活性（如激活态 vs. 非激活态）。
- **翻译后修饰**：磷酸化、泛素化等 PTM 未被建模。
- **MD 资源需求**：对溶剂化的约 4,000 残基复合物进行 100 ns 模拟需要 GPU 资源（单 GPU 需数天）。

## 常见问题与排错

> **说明**：以下方案均不需要 root/sudo 权限，只用 conda + pip 即可解决。

### 1. Blackwell GPU 特别说明（RTX 6000 / B100 / B200）

Blackwell 架构 + CUDA 13.2 太新，cuEquivariance 内核尚未适配。**必须用 `--no_kernels`。**

两步操作：

**（1）安装普通版 boltz**（不装 cuEquivariance）：

```bash
pip uninstall boltz[cuda] boltz cuequivariance-cuda-kernels -y
pip install boltz
```

**（2）boltz/config/config.yaml 中设置**：

```yaml
prediction:
  no_kernels: true    # ← Blackwell 必须设为 true
```

推理仍然用 GPU，只是不依赖 cuEquivariance 加速。速度稍慢但不影响精度。

### 2. `nvcc fatal: Unsupported gpu architecture 'compute_120'`

同上，跳过 cuEquivariance 即可。

### 3. GCC 版本太新（nvcc 不兼容）

CUDA 的 `nvcc` 编译器不支持新版 GCC。

**症状**：
```
error: unsupported GNU version! gcc versions later than 13 are not supported!
```

**解决方案**（不需要 sudo）：

```bash
# 用 conda 安装兼容的 GCC（conda 装包不需要 root 权限）
conda install -c conda-forge gcc=12 gxx=12 -y

# 设置环境变量，告诉 CUDA 用 conda 里的编译器
export CC=$CONDA_PREFIX/bin/gcc
export CXX=$CONDA_PREFIX/bin/g++
export CUDAHOSTCXX=$CONDA_PREFIX/bin/g++

# 重新安装
pip install cuequivariance-cuda-kernels --no-cache-dir
```

如果仍然失败，直接用方案 1（跳过内核，装普通版 boltz）。

### 3. `torch.cuda.is_available()` 返回 `False`

CUDA 13.2 需要 PyTorch 2.7+ nightly 版本。

**解决方案**（不需要 sudo）：

```bash
# 卸载旧版 torch
pip uninstall torch torchvision torchaudio -y

# 安装 CUDA 13.x 兼容的 nightly 版本
pip install torch torchvision --index-url https://download.pytorch.org/whl/nightly/cu128
```

验证：
```bash
python -c "import torch; print(torch.cuda.is_available())"
```

### 4. `out of memory`（显存不足）

TORC1 复合物约 4,140 aa，单次推理峰值显存约 30–50 GB。RTX 6000 Blackwell（97 GB VRAM）通常不会 OOM，但多个 sample 叠加可能溢出。

**解决方案**：

编辑 `boltz/config/config.yaml`，减少采样数：

```yaml
prediction:
  diffusion_samples: 1    # 先跑 1 个，确认不 OOM 再逐步增加
```

用 `nvidia-smi` 监控显存使用情况。

### 5. MSA 服务器连接失败

Boltz-2 默认连接远程 MSA 服务器生成 MSA。

**解决方案**：

- 确认云服务器有外网访问权限
- 检查是否需要认证令牌（在 config 中设置 `msa_server_auth`）
- 备选：用 conda 安装本地 MMseqs2，切到离线模式：

```bash
conda install -c bioconda mmseqs2 -y
```

然后修改 `boltz/config/config.yaml`：

```yaml
prediction:
  use_msa_server: false
  msa_mode: "mmseqs2"
```

### 6. 序列获取失败（UniProt/TAIR 超时）

Stage 1 可能因为网络原因无法从外部 API 获取序列。

**解决方案**：

手动下载 FASTA 序列文件，放到 `boltz/input_data/fasta/` 目录下，命名为 `{uniprot_id}.fasta`（如 `Q9FR53.fasta`），然后从 Stage 2 开始：

```bash
python pipeline.py --stage 2
```

### 7. pip 安装包时权限错误

如果遇到 `Permission denied` 错误，说明 pip 试图写入系统目录。

**解决方案**：

```bash
# 用 --user 装到用户目录
pip install --user <包名>

# 或者更推荐：用 conda 环境（包全装在用户自己的环境里，天然不需要 root）
conda activate torc1predict
pip install <包名>
```

确保所有操作都在 conda 环境内进行，永远不会碰到权限问题。

### 8. Blackwell 架构 + CUDA 13.2 兼容性速查表

| 包 | 状态 | 无需 sudo 的安装方式 |
|---|---|---|
| PyTorch | nightly 已支持 | `pip install torch --index-url https://download.pytorch.org/whl/nightly/cu128` |
| Flash Attention | 可能不支持 | `pip install flash-attn --no-build-isolation` |
| cuEquivariance | **未适配** | 跳过，`pip install boltz` 并在 config 中设 `no_kernels: true` |
| boltz（普通版） | ✅ 可用 | `pip install boltz` |
| NCCL | 通常已装好 | 无需操作 |
| triton | 需 3.0+ | `pip install triton --pre` |
| MMseqs2 | ✅ 可用 | `conda install -c bioconda mmseqs2` |

### 9. 快速环境验证

部署完成后，几条命令确认一切正常：

```bash
# CUDA 是否可用
python -c "import torch; print('CUDA:', torch.cuda.is_available()); print('VRAM:', torch.cuda.get_device_properties(0).total_memory / 1e9, 'GB')"

# boltz 可导入
python -c "import boltz; print('Boltz OK')"

# 干跑验证
python pipeline.py --dry-run
```

## 引用

如果在研究中使用本流程，请同时引用本流程和底层工具 Boltz-2。

## 许可证

MIT
