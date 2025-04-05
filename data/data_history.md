# LLA 数据以及计算说明

**总结**：理论法有问题，算不了；置换法无法处理有replicate的情况。

## marine bacteria OTU data files info
- `convert_otu.ipynb`: 将原始OTU数据转换为ELSA可处理的格式。
- `otu_matrix_thre20.txt`:  经过Z-score标准化处理后的以20%为阈值筛选OTU序列的文件。（可省略）
- `data_distribution.png`: 对于上面文件的数据可视化。直观判断其分布是否合理（符合正则化）。
- `data.txt`: 处理了表头、列行名的文件，可直接放入 lla_compute 进行LLA分析。
- `otu-timepoints.txt`: 存储时间点信息的文件。


## compute LLA of OTU data
### 2025-3-28
- 命令：见 `otu_analysis.sh`
- keypoints: 使用*理论法*，测试延迟0-3。结果和日志见data/otu_results，运行时间在每个日志文件的最后一行。
- 未计算出来：日志中都出现报错 :
> Local Liquid Association Analysis Tool
> delayLimit	fillMethod	pvalueMethod	dataFile	resultFile	repNum	spotNum	bootNum	transFunc	normMethod	precision
> 0	linear	theo	/app/elsa/data/data.txt	/app/elsa/data/lla_results/lla_results_d0.txt	1	101	0	simple	pnz	1000
> caution: q-value estimation error
> Error during analysis: local variable 'rpvalues' referenced before assignment
- 错误原因：当 p-values 数组为空或所有值都是 NaN 时，*rpvalues* 变量没有被正确赋值。因此`lsalib.py` 中的 storeyQvalue函数出现问题。
- 为什么会出现此错误？1.lsalib缺少处理p值为空数组的情况（暂时不改）2.原数据为什么会出现p值计算为空？
- 考虑方案：
    1. 添加输入文件检测，改变minOccur
    2. *如果不计算q值行不行？*

### 2025-3-30
- 对于3.28p值为空的分析
    -  `llalib.py` 中的minOccur为0.5。如果超过 50% 的数据是缺失值或零值，该变量就会被过滤掉
        - 创建`check_data.py`添加输入文件检测。

- 尝试用 `lla_compute.py` 计算 `sim_data.txt`。却发现*脚本无法正确计算有replicate的序列*。
    - 命令：'''docker run -d     --name "elsa_test"     -v /work1/wyx/aaa/elsa:/app/elsa     elsa-debug     python /app/elsa/lla/lla_compute.py     /app/elsa/lla/sim_data.txt     /app/elsa/lla/test_results/test_results.txt -d 0 -p theo -x 1000 -r 3 -s 4 -m 0.2 -t SD -f linear -n pnz
    - 结果：Error processing triplet (0,1,2): Mask and data not compatible: data size is 4, mask size is 12.
    No valid triplets found for analysis
    Finishing up...
    Time elapsed 1.02 seconds
- 没有 replicate 的序列则可以 置换 计算出来。结果见`simulated.txt`的计算结果`perm_test_results.txt`
    - 但是用理论法仍然不行
    > '''docker run -d     --name "elsa_test"     
    -v /work1/wyx/aaa/elsa:/app/elsa     elsa-debug     python /app/elsa/lla/lla_compute.py     /app/elsa/lla/simulated.txt     /app/elsa/lla/test_results/test_results.txt -d 0 -x 1000 -p theo -r 1 -s 50 -m 0.2 -t SD -f linear -n pnz
    > docker logs elsa_test
        lla_compute (rev: v2.0.1) - copyright Li Charlie Xia, lcxia@scut.edu.cn
        ...
        caution: q-value estimation error
        Error during analysis: local variable 'rpvalues' referenced before assignment
    - 如果只是minOccur的问题，不应该对于simdata也出错。因为模拟的是没有零值的，也就是所有数据都不会被minOccur筛掉。

### 2025-4-1
- 针对理论法计算出错的问题进行修改：
  - 增加了 `lla_compute_test.py` 和 `llalib_test.py`
    - 添加 `skip_qvalue` 参数，使得 `-p theo` 时可以跳过 q 值计算
  - 修复了编译和导入问题：
    - 修改 `compcore.py` 和 `setup.py` 中的 `extra_link_args`，补充上面提到的两个新增文件。
    - 重新构建镜像 `elsa-debug`

- 进行两组测试：
  1. 使用原始 `otu_analysis.sh`：
     - 仍采用理论法计算，和30号一样的方法进行计算，再次检验。
     - 结果：出现相同错误，见日志文档 'otu_results'
     - 状态：由于出错(和30号一模一样)，未生成任何结果表格。

  2. 使用新建 `otu_analysis2.sh`：
     - 采用 `lla_compute_test.py` 进行理论法计算
     - 跳过 q 值计算以验证是否能解决错误
     - 结果和日志保存在 `lla_results_test2` 和 `otu_results_test2`
     - 状态：计算完成，成功生成结果表格。
     - 计算时间：d0:4922s; d1:4832s; d2:4849s; d3:4307s

**otu_lla结果分析** （data/lla_results_test2）
tbc