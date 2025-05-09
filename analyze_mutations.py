import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
import glob
import tkinter as tk
from tkinter import filedialog

# Set up font for displaying characters
plt.rcParams['font.sans-serif'] = ['SimHei']  # For displaying Chinese characters if needed
plt.rcParams['axes.unicode_minus'] = False  # For displaying negative sign correctly

def create_output_dirs(file_path):
    """
    Create output directory structure for the input file
    
    Parameters:
    file_path: Input file path
    
    Returns:
    Dictionary containing paths to output directories
    """
    # 获取输入文件名和目录
    file_name = os.path.basename(file_path)
    file_name_no_ext = os.path.splitext(file_name)[0]
    input_dir = os.path.dirname(file_path)
    
    # 在输入文件所在目录创建以文件名命名的文件夹
    output_base_dir = os.path.join(input_dir, file_name_no_ext + '_analysis')
    
    # 创建子目录
    images_dir = os.path.join(output_base_dir, 'images')
    data_dir = os.path.join(output_base_dir, 'data')
    reports_dir = os.path.join(output_base_dir, 'reports')
    
    # 确保目录存在
    for dir_path in [output_base_dir, images_dir, data_dir, reports_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    # 返回各目录路径
    return {
        'base_dir': output_base_dir,
        'images_dir': images_dir,
        'data_dir': data_dir,
        'reports_dir': reports_dir,
        'file_name': file_name,
        'file_name_no_ext': file_name_no_ext
    }

def analyze_mutation_distribution(file_path):
    """
    分析突变位点在序列中的分布
    
    参数:
    file_path: XLS文件路径
    
    返回:
    无, 但会生成分析结果和图表
    """
    print(f"正在分析文件: {file_path}")
    
    # 创建输出目录
    dirs = create_output_dirs(file_path)
    
    # 读取数据
    df = pd.read_csv(file_path, sep='\t')
    
    # 查看数据基本情况
    print(f"数据形状: {df.shape}")
    print(f"序列数量: {df['Chrom'].nunique()}")
    
    # 确保'Pos'和'Muts'列是数值型
    df['Pos'] = pd.to_numeric(df['Pos'])
    df['Muts'] = pd.to_numeric(df['Muts'])
    df['Depths'] = pd.to_numeric(df['Depths'])
    
    # 计算每个序列的长度
    sequence_lengths = df.groupby('Chrom')['Pos'].max().to_dict()
    print(f"序列长度范围: {min(sequence_lengths.values())} - {max(sequence_lengths.values())}")
    
    # 计算各种突变率
    df['MutationRate'] = df['Muts'] / df['Depths']  # 总突变率
    df['ARate'] = df['Acount'] / df['Depths']  # toA突变率 (突变为A的比例)
    df['TRate'] = df['Tcount'] / df['Depths']  # toT突变率 (突变为T的比例)
    df['CRate'] = df['Ccount'] / df['Depths']  # toC突变率 (突变为C的比例)
    df['GRate'] = df['Gcount'] / df['Depths']  # toG突变率 (突变为G的比例)
    df['DelRate'] = df['delcount'] / df['Depths']  # 总缺失率
    df['TotalMutDelRate'] = (df['Muts'] + df['delcount']) / df['Depths']  # 总突变+缺失率
    
    # 1. 按模板碱基(Template)分类的缺失率
    # 创建各碱基的缺失数据框
    df_A_del = df[df['Template'] == 'A'].copy()
    df_T_del = df[df['Template'] == 'T'].copy()
    df_C_del = df[df['Template'] == 'C'].copy()
    df_G_del = df[df['Template'] == 'G'].copy()
    
    # 计算各碱基位置的缺失率
    df['A_DelRate'] = 0  # 默认值为0
    df['T_DelRate'] = 0
    df['C_DelRate'] = 0
    df['G_DelRate'] = 0
    
    # 将数据框合并
    df.loc[df['Template'] == 'A', 'A_DelRate'] = df_A_del['delcount'] / df_A_del['Depths']
    df.loc[df['Template'] == 'T', 'T_DelRate'] = df_T_del['delcount'] / df_T_del['Depths']
    df.loc[df['Template'] == 'C', 'C_DelRate'] = df_C_del['delcount'] / df_C_del['Depths']
    df.loc[df['Template'] == 'G', 'G_DelRate'] = df_G_del['delcount'] / df_G_del['Depths']
    
    # 2. 按模板碱基(Template)分类的突变率
    # 创建各碱基的突变数据框
    df_A_mut = df[df['Template'] == 'A'].copy()
    df_T_mut = df[df['Template'] == 'T'].copy()
    df_C_mut = df[df['Template'] == 'C'].copy()
    df_G_mut = df[df['Template'] == 'G'].copy()
    
    # 计算各碱基位置的突变率 (Ato, Tto, Cto, Gto)
    df['Ato_MutRate'] = 0  # A位置的突变率
    df['Tto_MutRate'] = 0  # T位置的突变率
    df['Cto_MutRate'] = 0  # C位置的突变率
    df['Gto_MutRate'] = 0  # G位置的突变率
    
    # 将数据框合并
    df.loc[df['Template'] == 'A', 'Ato_MutRate'] = df_A_mut['Muts'] / df_A_mut['Depths']
    df.loc[df['Template'] == 'T', 'Tto_MutRate'] = df_T_mut['Muts'] / df_T_mut['Depths']
    df.loc[df['Template'] == 'C', 'Cto_MutRate'] = df_C_mut['Muts'] / df_C_mut['Depths']
    df.loc[df['Template'] == 'G', 'Gto_MutRate'] = df_G_mut['Muts'] / df_G_mut['Depths']
    
    # 按位置统计所有序列的各种平均突变率
    position_rates = df.groupby('Pos').agg({
        'MutationRate': 'mean',
        'ARate': 'mean',
        'TRate': 'mean',
        'CRate': 'mean',
        'GRate': 'mean',
        'DelRate': 'mean',
        'TotalMutDelRate': 'mean',
        'A_DelRate': 'mean',
        'T_DelRate': 'mean',
        'C_DelRate': 'mean',
        'G_DelRate': 'mean',
        'Ato_MutRate': 'mean',
        'Tto_MutRate': 'mean',
        'Cto_MutRate': 'mean',
        'Gto_MutRate': 'mean'
    })
    
    # 保存按位置的突变率统计结果
    # position_rates.to_csv('position_mutation_rates.csv')  # 注释掉未归一化的结果保存
    
    # 注释掉未归一化的可视化部分
    '''
    # 可视化各种突变率沿序列位置的分布
    # 1. 总突变率分布
    plt.figure(figsize=(15, 6))
    plt.bar(position_rates.index, position_rates['MutationRate'])
    plt.xlabel('序列位置')
    plt.ylabel('平均突变率')
    plt.title('总突变率(ATCG)在不同序列位置的分布')
    plt.savefig(os.path.join(dirs['images_dir'], 'mutation_distribution.png'))
    
    # 2. 各碱基突变率分布图
    plt.figure(figsize=(15, 8))
    plt.plot(position_rates.index, position_rates['ARate'], 'r-', label='A突变率')
    plt.plot(position_rates.index, position_rates['TRate'], 'g-', label='T突变率')
    plt.plot(position_rates.index, position_rates['CRate'], 'b-', label='C突变率')
    plt.plot(position_rates.index, position_rates['GRate'], 'y-', label='G突变率')
    plt.xlabel('序列位置')
    plt.ylabel('碱基特异性突变率')
    plt.title('不同碱基在各序列位置的突变率')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(dirs['images_dir'], 'base_specific_mutation_rates.png'))
    
    # 3. 总缺失率分布图
    plt.figure(figsize=(15, 6))
    plt.bar(position_rates.index, position_rates['DelRate'])
    plt.xlabel('序列位置')
    plt.ylabel('平均总缺失率')
    plt.title('总缺失率在不同序列位置的分布')
    plt.savefig(os.path.join(dirs['images_dir'], 'deletion_rate_distribution.png'))
    
    # 3.1 按碱基分类的缺失率分布图
    plt.figure(figsize=(15, 8))
    plt.plot(position_rates.index, position_rates['A_DelRate'], 'r-', label='A位置缺失率')
    plt.plot(position_rates.index, position_rates['T_DelRate'], 'g-', label='T位置缺失率')
    plt.plot(position_rates.index, position_rates['C_DelRate'], 'b-', label='C位置缺失率')
    plt.plot(position_rates.index, position_rates['G_DelRate'], 'y-', label='G位置缺失率')
    plt.xlabel('序列位置')
    plt.ylabel('碱基特异性缺失率')
    plt.title('不同碱基位置的缺失率')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(dirs['images_dir'], 'base_specific_deletion_rates.png'))
    
    # 3.2 对比总缺失率和各碱基缺失率相加
    # 计算各碱基缺失率之和，应该等于总缺失率
    position_rates['Sum_Base_DelRate'] = position_rates['A_DelRate'] + position_rates['T_DelRate'] + \
                                       position_rates['C_DelRate'] + position_rates['G_DelRate']
    
    plt.figure(figsize=(15, 6))
    plt.plot(position_rates.index, position_rates['DelRate'], 'b-', label='总缺失率')
    plt.plot(position_rates.index, position_rates['Sum_Base_DelRate'], 'r--', label='各碱基缺失率之和')
    plt.xlabel('序列位置')
    plt.ylabel('缺失率')
    plt.title('总缺失率与各碱基缺失率之和的对比')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(dirs['images_dir'], 'total_vs_sum_deletion_rates.png'))
    
    # 4. 突变+缺失总率分布图
    plt.figure(figsize=(15, 6))
    plt.bar(position_rates.index, position_rates['TotalMutDelRate'])
    plt.xlabel('序列位置')
    plt.ylabel('突变+缺失总率')
    plt.title('突变+缺失总率在不同序列位置的分布')
    plt.savefig(os.path.join(dirs['images_dir'], 'total_mutation_deletion_rate.png'))
    
    # 5. 组合对比图
    plt.figure(figsize=(15, 8))
    plt.plot(position_rates.index, position_rates['MutationRate'], 'b-', label='总突变率(ATCG)')
    plt.plot(position_rates.index, position_rates['DelRate'], 'r-', label='缺失率')
    plt.plot(position_rates.index, position_rates['TotalMutDelRate'], 'g-', label='突变+缺失总率')
    plt.xlabel('序列位置')
    plt.ylabel('错误率')
    plt.title('不同错误类型在各序列位置的分布对比')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(dirs['images_dir'], 'error_rate_comparison.png'))
    '''
    
    # 打印额外的分析信息
    print("\n位置特异性突变率统计已保存到 'position_mutation_rates.csv'")
    print("各种突变类型的位置分布图已生成。")
    
    
    # 位置归一化分析
    print("\n正在进行位置归一化分析...")
    
    # 将所有序列长度归一化为100个位置点
    normalized_data = []
    
    for seq_name, seq_df in df.groupby('Chrom'):
        seq_length = sequence_lengths[seq_name]
        
        for _, row in seq_df.iterrows():
            # 计算归一化位置 (0-100)
            norm_pos = int((row['Pos'] / seq_length) * 100)
            if norm_pos == 0:  # 确保位置至少为1
                norm_pos = 1
            
            normalized_data.append({
                'Sequence': seq_name,
                'NormalizedPos': norm_pos,
                'MutationRate': row['MutationRate'],
                'ARate': row['ARate'],  # toA突变率
                'TRate': row['TRate'],  # toT突变率
                'CRate': row['CRate'],  # toC突变率
                'GRate': row['GRate'],  # toG突变率
                'DelRate': row['DelRate'],
                'TotalMutDelRate': row['TotalMutDelRate'],
                'A_DelRate': row['A_DelRate'],
                'T_DelRate': row['T_DelRate'],
                'C_DelRate': row['C_DelRate'],
                'G_DelRate': row['G_DelRate'],
                'Ato_MutRate': row['Ato_MutRate'],  # A位置突变率
                'Tto_MutRate': row['Tto_MutRate'],  # T位置突变率
                'Cto_MutRate': row['Cto_MutRate'],  # C位置突变率
                'Gto_MutRate': row['Gto_MutRate']   # G位置突变率
            })
    
    norm_df = pd.DataFrame(normalized_data)
    
    # 按归一化位置计算各类型平均突变率
    norm_pos_rates = norm_df.groupby('NormalizedPos').agg({
        'MutationRate': 'mean',
        'ARate': 'mean',
        'TRate': 'mean',
        'CRate': 'mean',
        'GRate': 'mean',
        'DelRate': 'mean',
        'TotalMutDelRate': 'mean',
        'A_DelRate': 'mean',
        'T_DelRate': 'mean',
        'C_DelRate': 'mean',
        'G_DelRate': 'mean',
        'Ato_MutRate': 'mean',
        'Tto_MutRate': 'mean',
        'Cto_MutRate': 'mean',
        'Gto_MutRate': 'mean'
    })
    
    # 可视化归一化位置的各种突变率分布
    
    # 1. 总突变率
    plt.figure(figsize=(15, 6))
    plt.plot(norm_pos_rates.index, norm_pos_rates['MutationRate'], marker='o', linestyle='-')
    plt.xlabel('归一化序列位置 (%)')
    plt.ylabel('平均总突变率')
    plt.title('归一化位置的总突变率分布')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(dirs['images_dir'], 'normalized_position_mutation_rate.png'))
    
    # 2. 突变目标碱基率(toA, toT, toC, toG)
    plt.figure(figsize=(15, 8))
    plt.plot(norm_pos_rates.index, norm_pos_rates['ARate'], 'r-', marker='o', label='突变为A(toA)')
    plt.plot(norm_pos_rates.index, norm_pos_rates['TRate'], 'g-', marker='o', label='突变为T(toT)')
    plt.plot(norm_pos_rates.index, norm_pos_rates['CRate'], 'b-', marker='o', label='突变为C(toC)')
    plt.plot(norm_pos_rates.index, norm_pos_rates['GRate'], 'y-', marker='o', label='突变为G(toG)')
    plt.xlabel('归一化序列位置 (%)')
    plt.ylabel('突变率')
    plt.title('归一化位置的突变目标碱基分布(toA/T/C/G)')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(dirs['images_dir'], 'normalized_position_to_base_specific_rates.png'))
    
    # 2.1 模板碱基突变率(Ato, Tto, Cto, Gto)
    plt.figure(figsize=(15, 8))
    plt.plot(norm_pos_rates.index, norm_pos_rates['Ato_MutRate'], 'r-', marker='o', label='A位置突变率(Ato)')
    plt.plot(norm_pos_rates.index, norm_pos_rates['Tto_MutRate'], 'g-', marker='o', label='T位置突变率(Tto)')
    plt.plot(norm_pos_rates.index, norm_pos_rates['Cto_MutRate'], 'b-', marker='o', label='C位置突变率(Cto)')
    plt.plot(norm_pos_rates.index, norm_pos_rates['Gto_MutRate'], 'y-', marker='o', label='G位置突变率(Gto)')
    plt.xlabel('归一化序列位置 (%)')
    plt.ylabel('突变率')
    plt.title('归一化位置的模板碱基突变率分布(Ato/Tto/Cto/Gto)')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(dirs['images_dir'], 'normalized_position_from_base_specific_rates.png'))
    
    # 3. 总缺失率
    plt.figure(figsize=(15, 6))
    plt.plot(norm_pos_rates.index, norm_pos_rates['DelRate'], 'r-', marker='o')
    plt.xlabel('归一化序列位置 (%)')
    plt.ylabel('平均总缺失率')
    plt.title('归一化位置的总缺失率分布')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(dirs['images_dir'], 'normalized_position_deletion_rate.png'))
    
    # 3.1 归一化位置的碱基特异性缺失率
    plt.figure(figsize=(15, 8))
    plt.plot(norm_pos_rates.index, norm_pos_rates['A_DelRate'], 'r-', marker='o', label='A位置缺失率')
    plt.plot(norm_pos_rates.index, norm_pos_rates['T_DelRate'], 'g-', marker='o', label='T位置缺失率')
    plt.plot(norm_pos_rates.index, norm_pos_rates['C_DelRate'], 'b-', marker='o', label='C位置缺失率')
    plt.plot(norm_pos_rates.index, norm_pos_rates['G_DelRate'], 'y-', marker='o', label='G位置缺失率')
    plt.xlabel('归一化序列位置 (%)')
    plt.ylabel('碱基特异性缺失率')
    plt.title('归一化位置的各碱基缺失率分布')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(dirs['images_dir'], 'normalized_position_base_specific_deletion_rates.png'))
    
    # 3.2 模板碱基突变率与缺失率的对比
    plt.figure(figsize=(15, 12))
    # A的突变与缺失
    plt.subplot(2, 2, 1)
    plt.plot(norm_pos_rates.index, norm_pos_rates['Ato_MutRate'], 'r-', marker='o', label='A位置突变率')
    plt.plot(norm_pos_rates.index, norm_pos_rates['A_DelRate'], 'r--', marker='o', label='A位置缺失率')
    plt.xlabel('归一化序列位置 (%)')
    plt.ylabel('错误率')
    plt.title('A位置的突变与缺失对比')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # T的突变与缺失
    plt.subplot(2, 2, 2)
    plt.plot(norm_pos_rates.index, norm_pos_rates['Tto_MutRate'], 'g-', marker='o', label='T位置突变率')
    plt.plot(norm_pos_rates.index, norm_pos_rates['T_DelRate'], 'g--', marker='o', label='T位置缺失率')
    plt.xlabel('归一化序列位置 (%)')
    plt.ylabel('错误率')
    plt.title('T位置的突变与缺失对比')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # C的突变与缺失
    plt.subplot(2, 2, 3)
    plt.plot(norm_pos_rates.index, norm_pos_rates['Cto_MutRate'], 'b-', marker='o', label='C位置突变率')
    plt.plot(norm_pos_rates.index, norm_pos_rates['C_DelRate'], 'b--', marker='o', label='C位置缺失率')
    plt.xlabel('归一化序列位置 (%)')
    plt.ylabel('错误率')
    plt.title('C位置的突变与缺失对比')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # G的突变与缺失
    plt.subplot(2, 2, 4)
    plt.plot(norm_pos_rates.index, norm_pos_rates['Gto_MutRate'], 'y-', marker='o', label='G位置突变率')
    plt.plot(norm_pos_rates.index, norm_pos_rates['G_DelRate'], 'y--', marker='o', label='G位置缺失率')
    plt.xlabel('归一化序列位置 (%)')
    plt.ylabel('错误率')
    plt.title('G位置的突变与缺失对比')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(dirs['images_dir'], 'normalized_position_base_mut_vs_del_comparison.png'))
    
    # 4. 总突变+缺失率
    plt.figure(figsize=(15, 6))
    plt.plot(norm_pos_rates.index, norm_pos_rates['TotalMutDelRate'], 'g-', marker='o')
    plt.xlabel('归一化序列位置 (%)')
    plt.ylabel('总突变+缺失率')
    plt.title('归一化位置的总突变+缺失率分布')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(dirs['images_dir'], 'normalized_position_total_rate.png'))
    
    # 5. 组合对比图
    plt.figure(figsize=(15, 8))
    plt.plot(norm_pos_rates.index, norm_pos_rates['MutationRate'], 'b-', marker='o', label='总突变率(ATCG)')
    plt.plot(norm_pos_rates.index, norm_pos_rates['DelRate'], 'r-', marker='o', label='缺失率')
    plt.plot(norm_pos_rates.index, norm_pos_rates['TotalMutDelRate'], 'g-', marker='o', label='总突变+缺失率')
    plt.xlabel('归一化序列位置 (%)')
    plt.ylabel('错误率')
    plt.title('归一化位置的各类错误率对比')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(dirs['images_dir'], 'normalized_position_error_types_comparison.png'))
    
    # 归一化位置的突变率分布
    data_file = os.path.join(dirs['data_dir'], 'normalized_position_mutation_rates.csv')
    norm_pos_rates.to_csv(data_file)
    
    # 生成HTML报告
    report_path = generate_html_report(file_path, dirs)
    
    print("分析完成，结果已保存到目录：")
    print(f"  基础目录: {dirs['base_dir']}")
    print(f"  图像目录: {dirs['images_dir']}")
    print(f"  数据目录: {dirs['data_dir']}")
    print(f"  报告目录: {dirs['reports_dir']}")
    print(f"HTML报告已生成：{report_path}")
    
    # 提示打开报告的方法
    print(f"您可以手动打开HTML报告并导出PDF：")
    print(f"  1. 使用浏览器打开: {report_path}")
    print(f"  2. 点击页面顶部的“打印报告为PDF”按钮")
    print(f"  3. 建议使用文件名: {dirs['file_name_no_ext']}_report.pdf")
    print(f"  4. 或直接使用命令: start {report_path}")

def generate_html_report(input_file_path, dirs, output_path=None):
    """
    生成包含所有分析结果的HTML报告
    
    参数:
    input_file_path: 输入数据文件路径
    dirs: 目录路径字典
    output_path: 输出报告路径，默认为dirs['reports_dir']目录下的'{文件名}_report.html'
    
    返回:
    生成的HTML报告路径
    """
    if output_path is None:
        report_filename = f"{dirs['file_name_no_ext']}_report.html"
        output_path = os.path.join(dirs['reports_dir'], report_filename)
    
    # 获取输入文件名
    file_name = os.path.basename(input_file_path)
    
    # HTML头部
    html = '''
    <!DOCTYPE html>
    <html>
    <head>
        <title>突变分析报告</title>
        <meta charset="utf-8">
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            .container { max-width: 1200px; margin: 0 auto; }
            h1 { text-align: center; color: #333; }
            .figure-container { margin: 20px 0; text-align: center; }
            .figure-container img { max-width: 100%; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }
            .figure-title { font-weight: bold; margin: 10px 0; }
            .section { margin: 30px 0; }
            table { width: 100%; border-collapse: collapse; }
            table, th, td { border: 1px solid #ddd; }
            th, td { padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; }
            .report-header { display: flex; justify-content: space-between; margin-bottom: 20px; }
            .logo { text-align: center; font-size: 24px; font-weight: bold; }
            .info { font-size: 14px; color: #666; }
            @media print {
                .page-break { page-break-before: always; }
                .no-print { display: none; }
            }
            .print-button { 
                display: inline-block; 
                padding: 10px 15px; 
                background: #4CAF50; 
                color: white; 
                border-radius: 4px; 
                cursor: pointer; 
                margin-bottom: 20px; 
                font-size: 16px;
                font-weight: bold;
            }
        </style>
    </head>
    <body>
        <div class="container">
            <div class="report-header">
                <div class="logo">突变分析报告</div>
                <div class="info">
                    分析日期：CURRENT_DATE<br>
                    分析文件：FILENAME<br>
                </div>
            </div>
            
            <button class="print-button no-print" onclick="window.print()">打印报告为PDF</button>
            
            <div class="section">
                <h2>1. 总体突变分析</h2>
                
                <div class="figure-container">
                    <img src="../images/normalized_position_mutation_rate.png" alt="总突变率分布">
                    <div class="figure-title">图1: 归一化位置的总突变率分布</div>
                </div>
                
                <div class="figure-container">
                    <img src="../images/normalized_position_total_rate.png" alt="总突变+缺失率分布">
                    <div class="figure-title">图2: 归一化位置的总突变+缺失率分布</div>
                </div>
            </div>
            
            <div class="section page-break">
                <h2>2. 目标碱基分析 (toA/T/C/G)</h2>
                
                <div class="figure-container">
                    <img src="../images/normalized_position_to_base_specific_rates.png" alt="目标碱基突变率">
                    <div class="figure-title">图3: 突变至特定碱基的比率分布 (toA/toT/toC/toG)</div>
                </div>
            </div>
            
            <div class="section page-break">
                <h2>3. 模板碱基分析 (Ato/Tto/Cto/Gto)</h2>
                
                <div class="figure-container">
                    <img src="../images/normalized_position_from_base_specific_rates.png" alt="模板碱基突变率">
                    <div class="figure-title">图4: 不同位置碱基的突变率分布 (Ato/Tto/Cto/Gto)</div>
                </div>
            </div>
            
            <div class="section page-break">
                <h2>4. 碱基特异性突变与缺失对比</h2>
                
                <div class="figure-container">
                    <img src="../images/normalized_position_base_mut_vs_del_comparison.png" alt="突变与缺失对比">
                    <div class="figure-title">图5: 各碱基位置的突变与缺失率对比</div>
                </div>
                
                <div class="figure-container">
                    <img src="../images/normalized_position_base_specific_deletion_rates.png" alt="碱基特异性缺失率">
                    <div class="figure-title">图6: 不同碱基位置的缺失率分布</div>
                </div>
            </div>
            
            <div class="section">
                <h2>5. 数据表格</h2>
                <p>归一化位置突变率数据可在 <a href="../data/normalized_position_mutation_rates.csv">normalized_position_mutation_rates.csv</a> 文件中找到。位置1-100 为 5p-3p</p>
            </div>
            
            <div class="section">
                <p style="text-align: center; color: #666; font-style: italic;">
                    此报告由突变分析脚本自动生成
                </p>
            </div>
        </div>
        
        <script>
            document.querySelector('.print-button').addEventListener('click', function() {
                window.print();
            });
        </script>
    </body>
    </html>
    '''
    
    # 替换特定标记
    html = html.replace('CURRENT_DATE', datetime.datetime.now().strftime('%Y-%m-%d'))
    html = html.replace('FILENAME', file_name)
    
    # 写入HTML文件
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    
    # 不自动打开HTML报告，只返回路径
    
    return output_path

def open_html_in_browser(html_path):
    """
    在浏览器中打开HTML报告，供用户直接另存为PDF
    
    参数:
    html_path: HTML报告路径
    
    返回:
    None
    """
    try:
        # 使用系统默认浏览器打开HTML文件
        html_path_url = 'file:///' + os.path.abspath(html_path).replace('\\', '/')
        print(f"打开HTML报告命令: start {html_path}")
        print("打开后请点击网页中的“打印报告为PDF”按钮，或使用浏览器的Ctrl+P功能导出PDF")
        # 不自动打开浏览器
        # webbrowser.open(html_path_url)
        return True
    except Exception as e:
        print(f"准备HTML报告时出错: {str(e)}")
        return False

def select_folder():
    """
    Open a folder selection dialog for the user to select a folder to analyze
    
    Returns:
    Selected folder path, or None if user cancels
    """
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    folder_path = filedialog.askdirectory(title="Select folder containing .mpileup.cns.filter.xls files")
    return folder_path if folder_path else None

def batch_process_files(folder_path):
    """
    Batch process all .mpileup.cns.filter.xls files in the specified folder
    
    Parameters:
    folder_path: Path to the folder to process
    
    Returns:
    Number of files processed
    """
    # Find all matching files
    pattern = os.path.join(folder_path, "*.mpileup.cns.filter.xls")
    files = glob.glob(pattern)
    
    if not files:
        print(f"ERROR: No .mpileup.cns.filter.xls files found in folder {folder_path}")
        return 0
    
    print(f"Found {len(files)} files to process...")
    
    # Process each file
    for i, file_path in enumerate(files):
        print(f"\n[{i+1}/{len(files)}] Processing: {os.path.basename(file_path)}")
        try:
            analyze_mutation_distribution(file_path)
        except Exception as e:
            print(f"Error processing file {file_path}: {str(e)}")
    
    print(f"\nBatch processing complete! Processed {len(files)} files")
    return len(files)

if __name__ == "__main__":
    print("Mutation Analysis Tool - Batch Processing Mode")
    print("Please select a folder containing .mpileup.cns.filter.xls files...")
    
    # Let user select a folder
    folder_path = select_folder()
    
    if folder_path:
        print(f"Selected folder: {folder_path}")
        # Batch process files
        batch_process_files(folder_path)
    else:
        print("User canceled operation or no folder was selected")
