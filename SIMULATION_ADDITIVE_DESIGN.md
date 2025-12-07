# 加性型（Additive）LLA 模拟设计方案

## 目标
构建加性（additive）型模拟数据生成机制，统一支持全局与局部调节窗口，用于评估 LLA 在不同效应大小下的检测力（Recall / FPR）、定位精度（窗口 IoU、边界误差）以及（可选）延迟恢复能力。取代旧“混合（mixing）”方式的优势：显式分离“基线噪声”与“信号”，信噪比 (SNR) 由参数直接解释，避免旧法中噪声与信号同阶导致的相关强度歧义。

## 共同设定
- 序列长度：`n`（基准常用 `n=100`，亦可扩展 `{20,40,60,80,100}`）。
- 调节窗口 \(I \subseteq \{1,\dots,n\}\)：
  - 全局：\(I = \{1,\dots,n\}\)（`is_global = True`）。
  - 局部：需满足 \(|I| \ge \max(10, n/2)\)；代码中强制校验。
  - 局部窗口按给定占比 `window_fraction=f` 居中：
    \[ w = \max(10, \lfloor n f \rfloor),\quad k = \left\lfloor \frac{n-w}{2} \right\rfloor,\quad I = \{k, k+1, \dots, k+w-1\} \]
- 调节指示变量：
  \[ Z_i = \mathbf{1}_{\{i \in I\}} \]
- 延迟（可选）：\( d = \text{delay}_{YZ} \ge 0 \)，表示窗口内信号引用 \(Y_{i+d}\)（当前多数实验设 \(d=0\)）。
- 复制次数：每个条件（效应大小 × 窗口类型 × 延迟值）进行 `n_replicates`（默认 50）实验组与对照组各自重复。

## 随机结构（加性方法核心）
1. 基线生成：
   \[ Y^{(\text{base})}_i \sim N(0, \sigma_Y^2),\quad X^{(\text{base})}_i \sim N(0, \sigma_{\epsilon}^2) \]
   - 参数：`sigma_y`（默认 1.0），`noise_sd`（默认 0.1）。
2. 信号幅度（效应大小）：`effect_size = E ≥ 0`（典型扫掠区间 `[0.8,1.0]`）。
3. 信号叠加（仅实验组且限窗口内）：
   \[ X_i = X^{(\text{base})}_i + E \cdot Y^{(\text{base})}_{i+d},\quad i \in I,\ 0 \le i+d < n \]
4. 窗口外：\(X_i = X^{(\text{base})}_i\)。
> 越界索引 \(i+d\) 不加信号（代码中边界检查）。

## 对照组生成
- 与实验组共享同一个 \(Z\)（窗口结构），但不添加信号：\(X_i = X^{(\text{base})}_i\)。
- \(Y\) 仍按 \(N(0, \sigma_Y^2)\) 独立生成。
- 目的：估计假阳性率 (FPR)，与实验组区分。

## 信噪比与参数解释
- 噪声标准差：\(\sigma_{\epsilon} = \text{noise\_sd}\)。
- 信号幅度：\(E = \text{effect\_size}\)。
- 近似瞬时信号对噪声比：\( \text{SNR} \approx E / \sigma_{\epsilon} \)。
  - 默认 `E=1.0, noise_sd=0.1` ⇒ \(\text{SNR} \approx 10\)。
- 与 mixing 法比较：旧法 \(X_i = \alpha Y_i + (1-\alpha)\epsilon_i\) 中噪声与信号同阶；新法幅度与噪声分离，可独立调节。

## 延迟处理
- 若 \(d>0\)：窗口内引用偏移后的 \(Y^{(\text{base})}_{i+d}\)。
- 当前基准：`delay_values=[0]`；可扩展测试 \(d \in \{0,1,2,\dots,D\}\)，并在 LLA 搜索使用 `--delay_limit D`。
- 越界时跳过信号叠加，避免边缘伪差。

## 参数推荐与扫掠策略
| 参数 | 默认 | 建议范围 | 说明 |
| ---- | ---- | -------- | ---- |
| `effect_size` | 0.8–1.0 | 0.4–1.2 | 扫描检测/定位曲线 |
| `noise_sd` | 0.1 | 0.05–0.2 | 控制 SNR |
| `sigma_y` | 1.0 | 0.5–1.5 | 基线波动规模 |
| `window_fraction` | 0.8 | 0.5–1.0 | 1.0=全局；局部长度下界 10 |
| `n_replicates` | 50 | 10 / 50 / 100 | 快速 / 正式 / 高精度 |
| `precision` | 1000 | 200–2000 | 置换次数影响 p 值稳定 |

## 文件与接口
- 生成器：`gen_unified_triplets.py`
  - 函数：`generate_unified_triplet(n, method='additive', effect_size, noise_sd, sigma_y, window_start, window_end, delay_yz, is_control, seed)` → `(X,Y,Z,metadata)`。
- 基准脚本：`benchmark_lla_unified.py`
  - 遍历 `(n, effect_size, delay, window_fraction)` 生成实验 + 对照，调用 `lla/lla_compute.py`。
  - 输出：
    - `benchmark_raw_results.csv`
    - `detection_summary.csv`（Recall, FPR, LLA 分数统计）
    - `localization_summary.csv`（`window_overlap`, `start_error`, `end_error`）
- 可视化：`visualize_effect_sweep.py` 与 `visualize_combined_results.py`（可合并多目录）。

## 评估指标
1. 检测性能：
   - 单次：`p < 0.05` 视为检出。
   - 汇总：\( \text{Recall} = \frac{\#\text{exp 检出}}{\#\text{exp 总数}} \), \( \text{FPR} = \frac{\#\text{ctrl 检出}}{\#\text{ctrl 总数}} \)。
2. 定位精度：
   - 边界误差：\(\text{start\_error} = \hat{s} - s\), \(\text{end\_error} = \hat{e} - e\)。
   - IoU：\( \text{IoU} = \frac{|I_{true} \cap I_{det}|}{|I_{true} \cup I_{det}|} \)。
   - 判据：IoU ≥ 0.9 或 |误差| ≤ 2 良好。
3. 延迟恢复（\(d>0\)）：\( \text{delay\_error} = \hat{d} - d \)。
4. 计算效率：每条件批次总耗时 `batch_time_seconds`。
5. 窗口膨胀：\( \text{inflation} = \hat{L} / L,\ \hat{L} \approx L + (-start\_error + end\_error) \)。目标 ≈ 1.0。

## 与旧 mixing 法差异
| 方面 | Mixing | Additive |
| ---- | ------ | -------- |
| 信号表达 | \(\alpha Y + (1-\alpha)\epsilon\) | 基线 + \(E Y\) |
| 噪声与信号尺度 | 同阶 | 可分离 \(E/\sigma_{\epsilon}\) |
| 参数可解释性 | \(\alpha\) 混合影响 | SNR \(\approx E/\sigma_{\epsilon}\) 直观 |
| 扩展到计数分布 | 不自然 | 仅替换基线分布 |
| 控制组一致性 | 需同时调整 \(\alpha\) | 固定噪声，去信号即可 |

## 推荐流程
1. 快速验证：
```bash
python benchmark_lla_unified.py --effect_size_values 0.8 1.0 --n_replicates 10 --precision 200 --output_dir quick_test
python visualize_effect_sweep.py --dir quick_test --save quick_test/plot.png
```
2. 正式基准：
```bash
python benchmark_lla_unified.py --effect_sweep --effect_start 0.8 --effect_end 1.0 --effect_step 0.02 \
  --n_replicates 50 --precision 1000 --output_dir additive_results
python visualize_effect_sweep.py --dir additive_results --save additive_results/effect_sweep.png
```
3. 合并多次运行：
```bash
python visualize_combined_results.py --dirs additive_results additive_results2 --save combined_effect_sweep.png
```

## 潜在扩展
- 非高斯：Poisson / NegBin 基线，保持加性结构。
- 多段窗口：\(Z\) 为多区间并集。
- 非线性：信号项换成 \(E \cdot g(Y_{i+d})\)（多项式 / 阈值）。

## 核心公式汇总
\[\begin{aligned}
Y^{(\text{base})}_i &\sim N(0, \sigma_Y^2)\\
X^{(\text{base})}_i &\sim N(0, \sigma_{\epsilon}^2)\\
Z_i &= \mathbf{1}_{\{i \in I\}}\\
X_i &= \begin{cases}
X^{(\text{base})}_i + E\,Y^{(\text{base})}_{i+d}, & i \in I,\ 0 \le i+d < n,\ \text{(实验组)}\\
X^{(\text{base})}_i, & \text{其他情况}
\end{cases}
\end{aligned}\]

## 结语
该加性设计改进了信噪分离与参数解释性，为后续扩展至非高斯分布和更复杂调节模式打下基础，可作为新版模拟基准的推荐标准。
