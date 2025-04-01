#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

def check_data(data_file):
    """
    检查数据文件的质量和特征
    
    Args:
        data_file: 输入数据文件路径
    """
    print(f"\n检查数据文件: {data_file}")
    print("=" * 50)
    
    # 读取数据
    try:
        # 读取数据，保留第一列作为索引
        data = pd.read_csv(data_file, sep='\t', index_col=0)
        print("\n1. 数据基本信息:")
        print(f"数据维度: {data.shape}")
        print(f"行数(变量数): {data.shape[0]}")
        print(f"列数(时间点): {data.shape[1]}")
        
        # 检查缺失值
        missing_stats = data.isna().sum()
        missing_percent = (missing_stats / len(data)) * 100
        print("\n2. 缺失值统计:")
        print(f"总缺失值数量: {data.isna().sum().sum()}")
        print(f"总缺失值比例: {(data.isna().sum().sum() / data.size) * 100:.2f}%")
        print("\n缺失值最多的前5个变量:")
        print(missing_percent.sort_values(ascending=False).head())
        
        # 检查零值
        zero_stats = (data == 0).sum()
        zero_percent = (zero_stats / len(data)) * 100
        print("\n3. 零值统计:")
        print(f"总零值数量: {(data == 0).sum().sum()}")
        print(f"总零值比例: {((data == 0).sum().sum() / data.size) * 100:.2f}%")
        print("\n零值最多的前5个变量:")
        print(zero_percent.sort_values(ascending=False).head())
        
        # 检查数据分布
        print("\n4. 数据分布统计:")
        print(data.describe())
        
        # 计算每个变量的有效值比例
        valid_percent = ((~data.isna()) & (data != 0)).sum() / len(data) * 100
        print("\n5. 有效值比例统计:")
        print(f"平均有效值比例: {valid_percent.mean():.2f}%")
        print(f"最小有效值比例: {valid_percent.min():.2f}%")
        print(f"最大有效值比例: {valid_percent.max():.2f}%")
        print("\n有效值比例最低的前5个变量:")
        print(valid_percent.sort_values().head())
        
        # 生成可视化图表
        output_dir = "data_quality_plots"
        os.makedirs(output_dir, exist_ok=True)
        
        # 1. 缺失值热图
        plt.figure(figsize=(15, 10))
        sns.heatmap(data.isna(), cmap='YlOrRd', cbar=False)
        plt.title('Missing Values Heatmap')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'missing_values_heatmap.png'))
        plt.close()
        
        # 2. 数据分布箱线图
        plt.figure(figsize=(15, 10))
        data.boxplot()
        plt.title('Data Distribution Boxplot')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'data_distribution_boxplot.png'))
        plt.close()
        
        # 3. 有效值比例条形图
        plt.figure(figsize=(15, 10))
        valid_percent.plot(kind='bar')
        plt.title('Valid Values Percentage by Variable')
        plt.xlabel('Variables')
        plt.ylabel('Valid Values Percentage (%)')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'valid_values_percentage.png'))
        plt.close()
        
        print(f"\n可视化图表已保存到 {output_dir} 目录")
        
    except Exception as e:
        print(f"错误: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("使用方法: python check_data.py <data_file>")
        sys.exit(1)
    
    data_file = sys.argv[1]
    check_data(data_file)