import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import glob

# 设置中文显示
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

# 导入主分析函数
from analyze_mutations import analyze_mutation_distribution

def main():
    """
    Entry function for processing a single file
    """
    # Check command line arguments
    if len(sys.argv) != 2:
        print("Usage: python analyze_mutations_single.py <file_path>")
        return
    
    file_path = sys.argv[1]
    
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"ERROR: File {file_path} does not exist")
        return
    
    # Check if file has the correct format
    if '.mpileup.cns' not in file_path:
        print(f"WARNING: File {file_path} does not contain '.mpileup.cns' in the filename")
        response = input("Continue processing? (y/n): ")
        if response.lower() != 'y':
            return
    
    # Run analysis
    try:
        analyze_mutation_distribution(file_path)
    except Exception as e:
        print(f"Error processing file: {str(e)}")

if __name__ == "__main__":
    main()
