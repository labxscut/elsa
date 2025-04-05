# Theo_Sim Table

**问题总结**：
生成表格的时候，实际值的结果总为0，也就是P(LA＞2)=0。经过进一步测试，用`lla_sim.py`生成随机数据后放入`lla_compute.py`进行分析，发现当序列数为100或200时，虽然产生了大量的三元组（100万多个）,但是只有不到百分之0.01的三元组LLA值≥1，更是0个三元组数值达到2.

## 2025-3-30

- 后台挂载 elsa-test，运行 @test_simulated_pvalues.py 检测200条序列的时候（想看是不是随着序列数增加LA的上限也增加，增加速度怎样）的LLA值。(已删除容器)

''' ls -l logs/test_simulated_pvalues_20250328_*.log & cat logs/ test_simulated_pvalues_20250328_20250328_024644.log 

分析：
1. **数据生成**
   - 时间点：10个
   - 变量数：200条
   - 生成时间：2025-03-28 02:46:44

2. **计算过程**
   - 总三元组数：1,313,400个（C(200,3)）
   - x = 1 时过滤后有效三元组：151个
   - 过滤率：约99.99%
   - 方法：置换法（默认）

3. **计算结果**
   - LA值范围：*[0.0415, 1.1979]* **仍然离2还有很远，但似乎比100条序列时的1.0458大一点。**
   - x=1.0时的p值：0.00011
   - 计算耗时：*46.6小时（167,754.34秒）*



## 2025-3-27
运行 @test_simulated_pvalues.py 从测试结果中，一些重要的信息：

1. LA 值的分布范围：
   ```
   LA value range: [0.0579, 1.0458]
   ```
   - 最小值：0.0579
   - 最大值：1.0458
   - 这个范围解释了为什么 x ≥ 1.5 时得到 0.0 的结果，因为 LA 值的最大值只有 1.0458

2. 三元组数量：
   ```
   Total triplets: 19600
   Filtered triplets: 2 (for x=1.0)
   Filtered triplets: 0 (for x≥1.5)
   ```
   - 总共有 19600 个三元组
   - 当 x=1.0 时，只有 2 个三元组的 LA 值大于等于阈值
   - 当 x≥1.5 时，没有三元组满足条件

3. 计算时间：
   ```
   Time elapsed ~2040 seconds
   ```
   - 每次计算大约需要 34 分钟
   - 这是正常的，因为 `lla_compute` 需要处理大量三元组

4. 参数设置：
   ```
   -m 0.5 (minOccur)
   -p perm (pvalueMethod)
   -x 1000 (precision)
   ```

问题在于：
1. LA 值的范围（0.0579-1.0458）太小，导致无法得到较大的 x 值的结果
2. 即使对于 x=1.0，也只有 2 个三元组满足条件，这可能表明：
   - 样本量（n_points=10）可能不够
   - 序列数（n_series=50）可能不够
   - 模拟数据的生成参数可能需要调整


- n_points=10，n_series=50时LA 值的范围（0.0579-1.0458）19,600个三元组，计算大约 34 分钟。当 x=1.0 时，只有 2 个三元组的 LA 值大于等于阈值；当 x≥1.5 时，没有三元组满足条件

- 改成了n_points=10, nserious=100。后台挂容器生成了161,700个三元组，计算时间约4.7小时 LA值范围：[0.0574, 1.1475]，因此结果仍为0 =>@comparison.log详细记录了整个生成表格的过程。

## 2025-3-26
发现 generate_comparison_table.py 期望从 lla_sim 生成的输出文件中读取包含'Score'列的数据但实际上输出文件中的列名是 'LA' 而不是 'Score'。已在generate中修复此问题。

## 2025-3-15
container name: theo_sim_new
command: docker run -d --name theo_sim_new 8b80cb8db986 python generate_comparison_table.py --xmin 2.0 --xmax 3.0 --xstep 0.5 --sample_sizes 10 20
UUID: 23c0648668c60115a5d06178b0b6b2eee8c00cc82d6b08eb57faf80952e37e99
This is not volumned! So we create a new volumned container:

container name: theo_sim_table_new
command: docker run -d --name theo_sim_table_new -v /work1/wyx/aaa/elsa/new_output:/tmp 8b80cb8db986 python generat
e_comparison_table.py --xmin 2.0 --xmax 3.0 --xstep 0.5 --sample_sizes 10 20
UUID: 5669a7d895088e3e9278f74ec3f69bdb6576ba2e9efddd30adfe912e1c70537f

## 2024-11-26
# If @sim_data.txt is what @lla_compute.py expects:
Let me analyze the format of your `sim_data.txt` and compare it with what `lla_compute.py` expects.

Your `sim_data.txt` has this format:
```
#	T1R1	T1R2	T1R3	T2R1	T2R2	T2R3	T3R1	T3R2	T3R3	T4R1	T4R2	T4R3
S1	-0.243721	-0.435952	-0.349710	-0.908051	...
S2	1.202188	1.314801	1.220819	0.412809	...
S3	0.225368	0.358246	0.226262	-1.238503	...
```

Breaking this down:
- 3 series (S1, S2, S3) - good for X, Y, Z variables
- 4 timepoints (T1-T4)
- 3 replicates per timepoint (R1-R3)

Looking at `lla_compute.py`, it expects input data in this exact format:
1. First row starts with '#' as header
2. Tab-delimited text
3. Format is m variables × (r replicates × s spots)
4. Data arranged as TiRj (Timepoint i, Replicate j)

Your `sim_data.txt` matches all these requirements:
1. ✓ Has '#' header
2. ✓ Is tab-delimited
3. ✓ Has 3 variables × (3 replicates × 4 spots)
4. ✓ Data arranged as TiRj format

Therefore, your `sim_data.txt` is already in the correct format for `lla_compute.py`. You can use it directly with `lla_compute.py`.

To run the analysis, you would use:
```bash
lla_compute sim_data.txt results.txt
```

The format matches because:
1. The data dimensions are correct:
   - factorNum = 3 (S1, S2, S3)
   - repNum = 3 (R1, R2, R3)
   - spotNum = 4 (T1, T2, T3, T4)

2. The data layout matches what `lla_compute.py` expects:
   ```python
   # In lla_compute.py, it reads data as:
   firstData = np.genfromtxt(args.dataFile, comments='#', delimiter='\t',
                          missing_values=['na','','NA'],
                          filling_values=np.nan,
                          usecols=list(range(1, args.spotNum*args.repNum+1)),
                          dtype=float)
   ```

This matches your data format exactly. The only potential adjustment needed would be if you want to use any optional parameters like:
- delayLimit (-d)
- pvalueMethod (-p)
- precision (-x)
- bootNum (-b)

But the data format itself is already correct for analysis.

# Key points about parameters in @lla_sim.py:
- lengthSeries: Total number of measurements. Must be divisible by repNum. Default: 50
- numSeries: Number of series to generate (min 3 for X,Y,Z variables). Default: 3
- repNum: Number of replicates per timepoint. Default: 1
- trendThresh: Threshold for trend analysis. If specified, generates trend series. Default: None
- simMethod: Simulation method: idn,mean,variance for independent normal distribution. Default: idn,0,1
- timeSeriesOutput: Output file path for generated time series data (tab-delimited). Default: None  

# Key points about output format of @lla_compute.py:
- Output is tab-delimited with:
  - Header row: # T1R1 T1R2 ... T2R1 T2R2 ... (Timepoint/Replicate labels)
  - Data rows: Si value1 value2 ... (Series i values)
  - Values are formatted to 6 decimal places

# Key points about input parameters of @lla_sim.py:
- lengthSeries: Total number of measurements. Must be divisible by repNum. Default: 50
- numSeries: Number of series to generate (min 3 for X,Y,Z variables). Default: 3
- repNum: Number of replicates per timepoint. Default: 1
- trendThresh: Threshold for trend analysis. If specified, generates trend series. Default: None
- simMethod: Simulation method: idn,mean,variance for independent normal distribution. Default: idn,0,1
- timeSeriesOutput: Output file path for generated time series data (tab-delimited). Default: None    

## 2024-11-29
# Issues and Recommendations
a. Code Duplication:
Many utility functions are duplicated between LSA and LLA
Recommendation: Create a shared utilities module
b. R Dependencies:
Both LSA and LLA query tools have commented out R dependencies:
   in @lsa_query.py
   Recommendation: Replace R dependencies with pure Python implementations
c. File Organization:
LLA package structure mirrors LSA but is less complete
Recommendation: Standardize the package structure between LSA and LLA
# Action Items:
1. Create shared utility module for common functions
2. Update Python 3 compatibility
3. Standardize error handling between LSA and LLA
4. Remove R dependencies or provide alternative implementations
5. Update documentation to reflect relationship between LSA and LLA
6. Add proper type hints and docstrings
7. Implement consistent testing framework