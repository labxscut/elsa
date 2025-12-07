#!/usr/bin/env python3
"""
LLA vs LSA 局部性差异演示
设计测试序列来展示两种方法的"局部"概念差异
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import sys
import os

# 添加路径以导入我们的模块
sys.path.append('lla')
sys.path.append('lsa')

try:
    from lla import llalib
    from lsa import lsalib
    from lsa import compcore
except ImportError:
    print("Warning: Could not import lla/lsa modules. Using manual calculations.")
    llalib = None
    lsalib = None
    compcore = None

def create_test_sequences():
    """
    创建测试序列：
    - 序列1: 有明显局部模式的序列
    - 序列2: 与序列1有局部相似性但全局相关性较弱
    - 序列3: 调节序列
    
    设计原则：
    1. 简单易懂
    2. 局部相似性明显
    3. 全局相似性较弱
    4. delay=0的情况下也能看出差异
    """
    
    # 时间点
    n_points = 20
    
    # 序列1: 有两个明显的局部模式
    X = np.array([
        1, 1, 1, 1, 1,    # 局部模式1: 高值段
        -2, -2, -2,       # 过渡段
        2, 2, 2, 2, 2,    # 局部模式2: 更高值段
        -1, -1,           # 过渡段
        1, 1, 1, 1        # 局部模式3: 高值段
    ])
    
    # 序列2: 与X有局部相似性，但全局相关性较弱
    Y = np.array([
        0.8, 0.9, 1.1, 0.9, 0.8,  # 与X前5个点相似
        -1.5, -1.8, -2.1,         # 与X中间段相似
        1.8, 2.1, 1.9, 2.0, 1.8,  # 与X中间段相似
        -0.5, -0.8,               # 与X不同
        -1.2, -1.0, -0.8, -1.1   # 与X相反模式
    ])
    
    # 序列3: 调节序列，影响局部关联
    Z = np.array([
        1, 1, 1, 1, 1,    # 与X、Y都正相关
        -1, -1, -1,       # 与X、Y都负相关
        1, 1, 1, 1, 1,    # 与X、Y都正相关
        1, 1,             # 与X正相关，与Y负相关
        1, 1, 1, 1        # 与X正相关，与Y负相关
    ])
    
    return X, Y, Z

def calculate_global_la(X, Y, Z):
    """计算全局LA"""
    return np.sum(X * Y * Z) / len(X)

def calculate_local_la_segments(X, Y, Z, segment_length=5):
    """
    手动计算局部LA（基于极差的概念）
    将序列分成若干段，每段计算LA，然后取最大值
    """
    n = len(X)
    local_las = []
    
    for i in range(n - segment_length + 1):
        x_seg = X[i:i+segment_length]
        y_seg = Y[i:i+segment_length]
        z_seg = Z[i:i+segment_length]
        
        la_seg = np.sum(x_seg * y_seg * z_seg) / segment_length
        local_las.append(la_seg)
    
    return local_las

def calculate_local_similarity_segments(X, Y, segment_length=5):
    """
    手动计算局部相似性（基于极差的概念）
    """
    n = len(X)
    local_similarities = []
    
    for i in range(n - segment_length + 1):
        x_seg = X[i:i+segment_length]
        y_seg = Y[i:i+segment_length]
        
        # 计算该段的相关系数
        corr, _ = pearsonr(x_seg, y_seg)
        local_similarities.append(abs(corr))
    
    return local_similarities

def run_lla_analysis(X, Y, Z):
    """运行LLA分析"""
    if llalib is None or compcore is None:
        return None, "LLA模块不可用"
    
    try:
        # 构造LLA数据
        lla_data = compcore.LLA_Data(0, X.tolist(), Y.tolist(), Z.tolist())
        lla_result = compcore.DP_lla(lla_data)
        
        return lla_result.score, "成功"
    except Exception as e:
        return None, f"LLA分析失败: {str(e)}"

def run_lsa_analysis(X, Y):
    """运行LSA分析"""
    if lsalib is None or compcore is None:
        return None, "LSA模块不可用"
    
    try:
        # 构造LSA数据
        lsa_data = compcore.LSA_Data(0, X.tolist(), Y.tolist())
        lsa_result = compcore.DP_lsa(lsa_data, True)
        
        return lsa_result.score, "成功"
    except Exception as e:
        return None, f"LSA分析失败: {str(e)}"

def create_comparison_table(X, Y, Z):
    """创建对比表格"""
    
    # 计算各种指标
    global_la = calculate_global_la(X, Y, Z)
    local_las = calculate_local_la_segments(X, Y, Z)
    max_local_la = max(local_las)
    
    global_corr_xy, _ = pearsonr(X, Y)
    local_corrs = calculate_local_similarity_segments(X, Y)
    max_local_corr = max(local_corrs)
    
    # 运行算法
    lla_score, lla_status = run_lla_analysis(X, Y, Z)
    lsa_score, lsa_status = run_lsa_analysis(X, Y)
    
    # 创建结果表格
    results = {
        '方法': [
            '全局LA (手动计算)',
            '局部LA-最大值 (手动计算)', 
            'LLA算法结果',
            '全局相关系数 (X-Y)',
            '局部相关系数-最大值 (X-Y)',
            'LSA算法结果'
        ],
        '数值': [
            f"{global_la:.4f}",
            f"{max_local_la:.4f}",
            f"{lla_score:.4f}" if lla_score is not None else "N/A",
            f"{global_corr_xy:.4f}",
            f"{max_local_corr:.4f}",
            f"{lsa_score:.4f}" if lsa_score is not None else "N/A"
        ],
        '状态': [
            "成功",
            "成功",
            lla_status,
            "成功", 
            "成功",
            lsa_status
        ],
        '说明': [
            "整个序列的三元乘积平均值",
            "所有局部段中LA的最大值",
            "LLA动态规划算法结果",
            "整个序列的相关系数",
            "所有局部段中相关系数的最大值",
            "LSA动态规划算法结果"
        ]
    }
    
    return pd.DataFrame(results)

def create_visualization(X, Y, Z):
    """创建可视化图表"""
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    # 子图1: 原始序列
    axes[0, 0].plot(X, 'o-', label='序列X', linewidth=2, markersize=6)
    axes[0, 0].plot(Y, 's-', label='序列Y', linewidth=2, markersize=6)
    axes[0, 0].plot(Z, '^-', label='序列Z', linewidth=2, markersize=6)
    axes[0, 0].set_title('原始时间序列', fontsize=14, fontweight='bold')
    axes[0, 0].set_xlabel('时间点')
    axes[0, 0].set_ylabel('数值')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 子图2: 局部LA分析
    local_las = calculate_local_la_segments(X, Y, Z)
    axes[0, 1].plot(local_las, 'o-', color='red', linewidth=2, markersize=6)
    axes[0, 1].set_title('局部LA分析 (段长=5)', fontsize=14, fontweight='bold')
    axes[0, 1].set_xlabel('起始位置')
    axes[0, 1].set_ylabel('LA值')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].axhline(y=max(local_las), color='red', linestyle='--', alpha=0.7, 
                       label=f'最大值: {max(local_las):.3f}')
    axes[0, 1].legend()
    
    # 子图3: 局部相关性分析
    local_corrs = calculate_local_similarity_segments(X, Y)
    axes[1, 0].plot(local_corrs, 'o-', color='blue', linewidth=2, markersize=6)
    axes[1, 0].set_title('局部相关性分析 (段长=5)', fontsize=14, fontweight='bold')
    axes[1, 0].set_xlabel('起始位置')
    axes[1, 0].set_ylabel('|相关系数|')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].axhline(y=max(local_corrs), color='blue', linestyle='--', alpha=0.7,
                       label=f'最大值: {max(local_corrs):.3f}')
    axes[1, 0].legend()
    
    # 子图4: 三元乘积可视化
    triple_product = X * Y * Z
    axes[1, 1].bar(range(len(triple_product)), triple_product, 
                   color=['green' if x > 0 else 'red' for x in triple_product],
                   alpha=0.7)
    axes[1, 1].set_title('三元乘积 (X×Y×Z)', fontsize=14, fontweight='bold')
    axes[1, 1].set_xlabel('时间点')
    axes[1, 1].set_ylabel('乘积值')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].axhline(y=0, color='black', linewidth=1)
    
    plt.tight_layout()
    return fig

def generate_detailed_analysis(X, Y, Z):
    """生成详细分析报告"""
    
    print("="*60)
    print("LLA vs LSA 局部性差异分析报告")
    print("="*60)
    
    print("\n1. 测试序列设计:")
    print(f"   序列长度: {len(X)}")
    print(f"   X: {X}")
    print(f"   Y: {Y}")
    print(f"   Z: {Z}")
    
    print("\n2. 关键发现:")
    
    # 全局 vs 局部对比
    global_la = calculate_global_la(X, Y, Z)
    local_las = calculate_local_la_segments(X, Y, Z)
    max_local_la = max(local_las)
    
    print(f"   全局LA: {global_la:.4f}")
    print(f"   最大局部LA: {max_local_la:.4f}")
    print(f"   差异倍数: {max_local_la/abs(global_la):.2f}x" if global_la != 0 else "全局LA为0")
    
    global_corr, _ = pearsonr(X, Y)
    local_corrs = calculate_local_similarity_segments(X, Y)
    max_local_corr = max(local_corrs)
    
    print(f"   全局相关系数: {global_corr:.4f}")
    print(f"   最大局部相关系数: {max_local_corr:.4f}")
    print(f"   差异倍数: {max_local_corr/abs(global_corr):.2f}x" if global_corr != 0 else "全局相关为0")
    
    print("\n3. 局部段分析:")
    segment_length = 5
    for i, la in enumerate(local_las):
        start_idx = i
        end_idx = i + segment_length
        print(f"   段{i+1} (位置{start_idx}-{end_idx-1}): LA = {la:.4f}")
    
    print("\n4. 结论:")
    print("   - 全局分析可能掩盖局部强关联")
    print("   - 局部分析能发现时间序列中的'热点'区域")
    print("   - LLA和LSA的'局部'概念存在本质差异")

def main():
    """主函数"""
    
    print("开始LLA vs LSA局部性差异演示...")
    
    # 1. 创建测试序列
    X, Y, Z = create_test_sequences()
    
    # 2. 生成对比表格
    print("\n生成对比分析表格...")
    comparison_df = create_comparison_table(X, Y, Z)
    print(comparison_df.to_string(index=False))
    
    # 3. 创建可视化
    print("\n生成可视化图表...")
    fig = create_visualization(X, Y, Z)
    plt.savefig('lla_vs_lsa_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 4. 保存结果
    comparison_df.to_csv('comparison_results.csv', index=False)
    
    # 5. 生成详细报告
    generate_detailed_analysis(X, Y, Z)
    
    print(f"\n结果已保存:")
    print(f"- 对比表格: comparison_results.csv")
    print(f"- 可视化图表: lla_vs_lsa_analysis.png")

if __name__ == "__main__":
    main()
