#!/usr/bin/env python3
"""
简化版LLA vs LSA演示
不依赖外部模块，纯手动计算
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

def create_test_sequences():
    """创建测试序列"""
    
    # 序列1: 有明显局部模式的序列
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
    
    # 序列3: 调节序列
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
    """计算所有局部段的LA"""
    n = len(X)
    local_las = []
    
    for i in range(n - segment_length + 1):
        x_seg = X[i:i+segment_length]
        y_seg = Y[i:i+segment_length]
        z_seg = Z[i:i+segment_length]
        
        la_seg = np.sum(x_seg * y_seg * z_seg) / segment_length
        local_las.append(la_seg)
    
    return local_las

def calculate_local_correlation_segments(X, Y, segment_length=5):
    """计算所有局部段的相关系数"""
    n = len(X)
    local_corrs = []
    
    for i in range(n - segment_length + 1):
        x_seg = X[i:i+segment_length]
        y_seg = Y[i:i+segment_length]
        
        corr, _ = pearsonr(x_seg, y_seg)
        local_corrs.append(abs(corr))
    
    return local_corrs

def simulate_lsa_concept(X, Y, segment_length=5):
    """
    模拟LSA概念：寻找局部相似性极差
    这里我们计算每个局部段的相关性，然后找最大值
    """
    local_corrs = calculate_local_correlation_segments(X, Y, segment_length)
    return max(local_corrs), local_corrs

def simulate_lla_concept(X, Y, Z, segment_length=5):
    """
    模拟LLA概念：寻找局部关联性极差
    这里我们计算每个局部段的LA，然后找最大值
    """
    local_las = calculate_local_la_segments(X, Y, Z, segment_length)
    return max(local_las), local_las

def create_analysis_table(X, Y, Z):
    """创建分析结果表格"""
    
    # 全局指标
    global_la = calculate_global_la(X, Y, Z)
    global_corr_xy, _ = pearsonr(X, Y)
    
    # 局部指标
    max_local_la, all_local_las = simulate_lla_concept(X, Y, Z)
    max_local_corr, all_local_corrs = simulate_lsa_concept(X, Y)
    
    # 创建结果表格
    results = {
        '分析类型': [
            '全局LA',
            '最大局部LA', 
            '全局相关系数(X-Y)',
            '最大局部相关系数(X-Y)'
        ],
        '数值': [
            f"{global_la:.4f}",
            f"{max_local_la:.4f}",
            f"{global_corr_xy:.4f}",
            f"{max_local_corr:.4f}"
        ],
        '差异倍数': [
            "基准",
            f"{max_local_la/abs(global_la):.2f}x" if global_la != 0 else "N/A",
            "基准", 
            f"{max_local_corr/abs(global_corr_xy):.2f}x" if global_corr_xy != 0 else "N/A"
        ],
        '说明': [
            "整个序列的三元乘积平均值",
            "所有局部段中LA的最大值",
            "整个序列的相关系数",
            "所有局部段中相关系数的最大值"
        ]
    }
    
    return pd.DataFrame(results), all_local_las, all_local_corrs

def create_visualization(X, Y, Z, all_local_las, all_local_corrs):
    """创建可视化"""
    
    plt.style.use('default')
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    # 子图1: 原始序列
    axes[0, 0].plot(X, 'o-', label='序列X', linewidth=2, markersize=6, color='blue')
    axes[0, 0].plot(Y, 's-', label='序列Y', linewidth=2, markersize=6, color='red')
    axes[0, 0].plot(Z, '^-', label='序列Z', linewidth=2, markersize=6, color='green')
    axes[0, 0].set_title('原始时间序列', fontsize=14, fontweight='bold')
    axes[0, 0].set_xlabel('时间点')
    axes[0, 0].set_ylabel('数值')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 子图2: 局部LA分析
    axes[0, 1].plot(all_local_las, 'o-', color='purple', linewidth=2, markersize=6)
    axes[0, 1].set_title('局部LA分析 (段长=5)', fontsize=14, fontweight='bold')
    axes[0, 1].set_xlabel('起始位置')
    axes[0, 1].set_ylabel('LA值')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].axhline(y=max(all_local_las), color='purple', linestyle='--', alpha=0.7, 
                       label=f'最大值: {max(all_local_las):.3f}')
    axes[0, 1].legend()
    
    # 子图3: 局部相关性分析
    axes[1, 0].plot(all_local_corrs, 'o-', color='orange', linewidth=2, markersize=6)
    axes[1, 0].set_title('局部相关性分析 (段长=5)', fontsize=14, fontweight='bold')
    axes[1, 0].set_xlabel('起始位置')
    axes[1, 0].set_ylabel('|相关系数|')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].axhline(y=max(all_local_corrs), color='orange', linestyle='--', alpha=0.7,
                       label=f'最大值: {max(all_local_corrs):.3f}')
    axes[1, 0].legend()
    
    # 子图4: 三元乘积可视化
    triple_product = X * Y * Z
    colors = ['green' if x > 0 else 'red' for x in triple_product]
    axes[1, 1].bar(range(len(triple_product)), triple_product, color=colors, alpha=0.7)
    axes[1, 1].set_title('三元乘积 (X×Y×Z)', fontsize=14, fontweight='bold')
    axes[1, 1].set_xlabel('时间点')
    axes[1, 1].set_ylabel('乘积值')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].axhline(y=0, color='black', linewidth=1)
    
    plt.tight_layout()
    return fig

def generate_detailed_report(X, Y, Z, all_local_las, all_local_corrs):
    """生成详细报告"""
    
    print("="*70)
    print("LLA vs LSA 局部性概念差异分析报告")
    print("="*70)
    
    print(f"\n📊 测试序列设计:")
    print(f"   序列长度: {len(X)}")
    print(f"   X (蓝色): {X}")
    print(f"   Y (红色): {Y}")
    print(f"   Z (绿色): {Z}")
    
    # 计算关键指标
    global_la = calculate_global_la(X, Y, Z)
    max_local_la = max(all_local_las)
    global_corr, _ = pearsonr(X, Y)
    max_local_corr = max(all_local_corrs)
    
    print(f"\n🔍 关键发现:")
    print(f"   📈 全局LA: {global_la:.4f}")
    print(f"   📈 最大局部LA: {max_local_la:.4f}")
    print(f"   📊 差异倍数: {max_local_la/abs(global_la):.2f}x" if global_la != 0 else "   全局LA为0")
    
    print(f"   📈 全局相关系数: {global_corr:.4f}")
    print(f"   📈 最大局部相关系数: {max_local_corr:.4f}")
    print(f"   📊 差异倍数: {max_local_corr/abs(global_corr):.2f}x" if global_corr != 0 else "   全局相关为0")
    
    print(f"\n📋 局部段详细分析:")
    segment_length = 5
    for i, (la, corr) in enumerate(zip(all_local_las, all_local_corrs)):
        start_idx = i
        end_idx = i + segment_length
        print(f"   段{i+1} (位置{start_idx}-{end_idx-1}):")
        print(f"      LA = {la:.4f}, 相关系数 = {corr:.4f}")
        print(f"      X段 = {X[start_idx:end_idx]}")
        print(f"      Y段 = {Y[start_idx:end_idx]}")
        print(f"      Z段 = {Z[start_idx:end_idx]}")
    
    print(f"\n💡 重要结论:")
    print(f"   1. 全局分析可能掩盖局部强关联")
    print(f"   2. 局部分析能发现时间序列中的'热点'区域")
    print(f"   3. LLA的'局部'概念：基于全局平均的局部最优路径")
    print(f"   4. LSA的'局部'概念：基于极差的局部相似性区域")
    print(f"   5. 两种方法在'局部'定义上存在本质差异")

def main():
    """主函数"""
    
    print("🚀 开始LLA vs LSA局部性差异演示...")
    
    # 1. 创建测试序列
    X, Y, Z = create_test_sequences()
    
    # 2. 生成分析表格
    print("\n📊 生成对比分析表格...")
    analysis_df, all_local_las, all_local_corrs = create_analysis_table(X, Y, Z)
    print(analysis_df.to_string(index=False))
    
    # 3. 创建可视化
    print("\n📈 生成可视化图表...")
    fig = create_visualization(X, Y, Z, all_local_las, all_local_corrs)
    plt.savefig('lla_vs_lsa_simple_demo.png', dpi=300, bbox_inches='tight')
    print("   图表已保存: lla_vs_lsa_simple_demo.png")
    
    # 4. 保存结果
    analysis_df.to_csv('simple_comparison_results.csv', index=False)
    print("   表格已保存: simple_comparison_results.csv")
    
    # 5. 生成详细报告
    generate_detailed_report(X, Y, Z, all_local_las, all_local_corrs)
    
    print(f"\n✅ 分析完成！结果文件:")
    print(f"   📊 对比表格: simple_comparison_results.csv")
    print(f"   📈 可视化图表: lla_vs_lsa_simple_demo.png")

if __name__ == "__main__":
    main()
