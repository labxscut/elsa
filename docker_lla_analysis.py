#!/usr/bin/env python3
"""
Docker环境下的完整LLA分析脚本
"""

import numpy as np
import pandas as pd
import subprocess
import os
import sys

def generate_simulation_data():
    """生成模拟数据"""
    # Set random seed for reproducibility
    np.random.seed(33)
    
    # Generate z as an alternating Bernoulli sequence
    n = 20
    z = np.random.randint(0, 2, size=n) * 2 - 1
    
    # generate x and y depending on z （no noise, local)
    x, y = np.zeros(n), np.zeros(n)  # 初始化x和y数组
    
    corr_start, corr_end = 3, 10 # 在（3，10）以外的区间非随机
    
    for i in range(n):
        if corr_start <= i <= corr_end:
            # 区间内：使用版本1的正负相关逻辑（无噪声）
            if z[i] == 1:
                base = np.random.normal(0, 1)
                x[i] = base
                y[i] = base  # 正相关
            else:
                base = np.random.normal(0, 1)
                x[i] = base
                y[i] = -base # 负相关
        else:
            # 区间外：x和y完全随机，无相关性
            x[i] = np.random.normal(0, 1)  # 随机值
            y[i] = np.random.normal(0, 1)  # 随机值（与x无关）
    
    return x, y, z

def create_lla_input_file(x, y, z, filename="lla_input.txt"):
    """创建LLA分析所需的输入文件"""
    
    # 创建头部行（以#开头）
    header = "#Factor\tT1\tT2\tT3\tT4\tT5\tT6\tT7\tT8\tT9\tT10\tT11\tT12\tT13\tT14\tT15\tT16\tT17\tT18\tT19\tT20"
    
    # 创建数据行
    data_lines = []
    data_lines.append(f"X\t{'\t'.join([f'{val:.6f}' for val in x])}")
    data_lines.append(f"Y\t{'\t'.join([f'{val:.6f}' for val in y])}")
    data_lines.append(f"Z\t{'\t'.join([f'{val:.6f}' for val in z])}")
    
    # 写入文件
    with open(filename, 'w') as f:
        f.write(header + '\n')
        for line in data_lines:
            f.write(line + '\n')
    
    print(f"✓ LLA输入文件已创建: {filename}")
    return filename

def run_lla_analysis(input_file, output_file="lla_results.txt"):
    """运行LLA分析"""
    
    print("开始运行LLA分析...")
    
    # 构建命令
    cmd = [
        "python3", "lla/lla_compute.py",
        input_file, output_file,
        "-d", "3",           # delayLimit
        "-p", "perm",        # pvalueMethod
        "-x", "1000",        # precision
        "-r", "1",           # repNum
        "-s", "20",          # spotNum
        "-m", "0.5",         # minOccur
        "-t", "simple",      # transFunc
        "-f", "linear",      # fillMethod
        "-n", "pnz"          # normMethod
    ]
    
    try:
        # 运行命令
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        print("✓ LLA分析完成")
        print(f"输出文件: {output_file}")
        
        # 显示输出
        if result.stdout:
            print("标准输出:")
            print(result.stdout)
        
        if result.stderr:
            print("标准错误:")
            print(result.stderr)
            
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"✗ LLA分析失败: {e}")
        print(f"错误输出: {e.stderr}")
        return False
    except FileNotFoundError:
        print("✗ 找不到lla_compute.py文件")
        print("请确保在正确的Docker环境中运行")
        return False

def analyze_results(output_file):
    """分析LLA结果"""
    
    if not os.path.exists(output_file):
        print(f"✗ 结果文件不存在: {output_file}")
        return
    
    print(f"\n=== 分析结果文件: {output_file} ===")
    
    try:
        # 读取结果
        with open(output_file, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 2:
            print("结果文件为空或格式不正确")
            return
        
        # 显示头部
        header = lines[0].strip().split('\t')
        print(f"结果列: {header}")
        
        # 显示数据行
        print("\n结果数据:")
        for i, line in enumerate(lines[1:], 1):
            if line.strip():
                parts = line.strip().split('\t')
                print(f"第{i}行: {parts}")
        
        # 尝试解析为DataFrame进行进一步分析
        try:
            df = pd.read_csv(output_file, sep='\t')
            print(f"\n数据框形状: {df.shape}")
            print("\n数据框预览:")
            print(df.head())
            
            # 如果有LA列，显示统计信息
            if 'LA' in df.columns:
                print(f"\nLA分数统计:")
                print(f"平均LA: {df['LA'].mean():.6f}")
                print(f"最大LA: {df['LA'].max():.6f}")
                print(f"最小LA: {df['LA'].min():.6f}")
                
        except Exception as e:
            print(f"解析DataFrame时出错: {e}")
            
    except Exception as e:
        print(f"读取结果文件时出错: {e}")

def create_visualization_script(x, y, z):
    """创建可视化脚本"""
    
    script_content = f"""
import numpy as np
import matplotlib.pyplot as plt

# 数据
x = {x.tolist()}
y = {y.tolist()}
z = {z.tolist()}

# 创建图表
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

# 原始时间序列
axes[0, 0].plot(x, 'b-', label='X', marker='o')
axes[0, 0].plot(y, 'r-', label='Y', marker='s')
axes[0, 0].plot(z, 'g-', label='Z', marker='^')
axes[0, 0].set_title('原始时间序列')
axes[0, 0].set_xlabel('时间点')
axes[0, 0].set_ylabel('值')
axes[0, 0].legend()
axes[0, 0].grid(True)

# X-Y散点图
axes[0, 1].scatter(x, y, c=z, cmap='coolwarm')
axes[0, 1].set_title('X vs Y (颜色表示Z)')
axes[0, 1].set_xlabel('X')
axes[0, 1].set_ylabel('Y')
axes[0, 1].grid(True)

# 三元乘积
xyz = np.array(x) * np.array(y) * np.array(z)
axes[1, 0].plot(xyz, 'purple', marker='o')
axes[1, 0].set_title('三元乘积 (X×Y×Z)')
axes[1, 0].set_xlabel('时间点')
axes[1, 0].set_ylabel('X×Y×Z')
axes[1, 0].grid(True)

# 局部相关性分析
segment_length = 5
local_corrs = []
for i in range(len(x) - segment_length + 1):
    x_seg = x[i:i+segment_length]
    y_seg = y[i:i+segment_length]
    corr = np.corrcoef(x_seg, y_seg)[0, 1]
    local_corrs.append(corr)

axes[1, 1].plot(range(len(local_corrs)), local_corrs, 'orange', marker='s')
axes[1, 1].set_title('局部相关性 (段长=5)')
axes[1, 1].set_xlabel('段起始位置')
axes[1, 1].set_ylabel('相关系数')
axes[1, 1].grid(True)

plt.tight_layout()
plt.savefig('simulation_analysis.png', dpi=300, bbox_inches='tight')
plt.show()

print("可视化图表已保存为: simulation_analysis.png")
"""
    
    with open("visualize_results.py", 'w') as f:
        f.write(script_content)
    
    print("✓ 可视化脚本已创建: visualize_results.py")

def main():
    """主函数"""
    print("=== Docker环境下的LLA分析 ===")
    
    # 1. 生成模拟数据
    print("1. 生成模拟数据...")
    x, y, z = generate_simulation_data()
    
    # 2. 创建输入文件
    print("2. 创建LLA输入文件...")
    input_file = create_lla_input_file(x, y, z)
    
    # 3. 运行LLA分析
    print("3. 运行LLA分析...")
    output_file = "lla_results.txt"
    success = run_lla_analysis(input_file, output_file)
    
    if success:
        # 4. 分析结果
        print("4. 分析结果...")
        analyze_results(output_file)
        
        # 5. 创建可视化脚本
        print("5. 创建可视化脚本...")
        create_visualization_script(x, y, z)
        
        print("\n=== 分析完成 ===")
        print("生成的文件:")
        print("- lla_input.txt: 输入数据")
        print("- lla_results.txt: LLA分析结果")
        print("- visualize_results.py: 可视化脚本")
        print("\n运行可视化: python3 visualize_results.py")
    else:
        print("\n=== 分析失败 ===")
        print("请检查Docker环境和文件路径")

if __name__ == "__main__":
    main()

