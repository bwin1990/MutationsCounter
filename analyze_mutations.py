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

BASES = ['A', 'T', 'C', 'G']

def _require_columns(df, required_cols):
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

def _safe_divide(num, den):
    try:
        if den == 0 or pd.isna(den):
            return np.nan
        return num / den
    except Exception:
        return np.nan

def compute_global_summary(df):
    """
    Compute weighted (sum-based) global substitution/deletion rates and deletion-by-template.
    """
    _require_columns(df, ['Depths', 'Muts', 'delcount', 'Template'])

    total_depths = float(df['Depths'].sum())
    total_muts = float(df['Muts'].sum())
    total_del = float(df['delcount'].sum())

    summary = {
        'total_depths': total_depths,
        'total_substitution_muts': total_muts,
        'total_deletions': total_del,
        'substitution_total_rate': _safe_divide(total_muts, total_depths),
        'deletion_total_rate': _safe_divide(total_del, total_depths),
        'total_error_rate': _safe_divide(total_muts + total_del, total_depths),
    }

    del_by_template = (
        df.groupby('Template')[['delcount', 'Depths']]
          .sum()
          .reindex(BASES)
    )
    del_by_template['DelRate'] = del_by_template['delcount'] / del_by_template['Depths']
    del_by_template.index.name = 'Template'

    return summary, del_by_template

def compute_substitution_matrix(df):
    """
    Compute from->to substitution count matrix and percent matrix (off-diagonal sums to 100%).
    """
    _require_columns(df, ['Template'])
    count_cols = [f"{b}count" for b in BASES]
    _require_columns(df, count_cols)

    by_from = df.groupby('Template')[count_cols].sum().reindex(BASES, fill_value=0)
    counts = pd.DataFrame(0.0, index=BASES, columns=BASES)
    for b in BASES:
        counts[b] = by_from[f"{b}count"].astype(float).values

    # Do not report diagonal as "substitution to self"
    np.fill_diagonal(counts.values, 0.0)

    total_sub_counts = float(np.nansum(counts.values))
    if total_sub_counts > 0:
        percent = counts / total_sub_counts * 100.0
    else:
        percent = counts.copy()

    return counts, percent, total_sub_counts

def compute_position_base_specific_deletion(df):
    """
    For each Pos, compute weighted deletion rate by template base (A/T/C/G).
    """
    _require_columns(df, ['Pos', 'Template', 'delcount', 'Depths'])
    grouped = df.groupby(['Pos', 'Template'])[['delcount', 'Depths']].sum()
    grouped['DelRate'] = grouped['delcount'] / grouped['Depths']
    rates = grouped['DelRate'].unstack('Template').reindex(columns=BASES)
    rates.index.name = 'Pos'
    return rates

def plot_substitution_matrix_heatmap(percent_matrix, summary, out_path):
    data = percent_matrix.values.astype(float)
    mask = np.eye(len(BASES), dtype=bool)
    data_masked = np.ma.array(data, mask=mask)

    plt.figure(figsize=(7, 6))
    cmap = plt.cm.YlOrRd.copy()
    cmap.set_bad(color='#eeeeee')
    vmax = np.nanmax(data) if np.isfinite(data).any() else 1
    if vmax == 0:
        vmax = 1
    im = plt.imshow(data_masked, cmap=cmap, vmin=0, vmax=vmax)
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Percent of substitutions (%)')
    plt.xticks(range(len(BASES)), BASES)
    plt.yticks(range(len(BASES)), BASES)
    plt.xlabel('To (mutated base)')
    plt.ylabel('From (template base)')
    sub_rate = summary.get('substitution_total_rate', np.nan)
    plt.title(f"Substitution Spectrum (Total substitution rate: {sub_rate:.6g})")

    for i, fb in enumerate(BASES):
        for j, tb in enumerate(BASES):
            if i == j:
                plt.text(j, i, '—', ha='center', va='center', color='#888888', fontsize=12)
                continue
            v = percent_matrix.loc[fb, tb]
            if pd.isna(v):
                txt = 'NA'
            else:
                txt = f"{v:.1f}%"
            color = 'white' if (not pd.isna(v) and v >= 15) else 'black'
            plt.text(j, i, txt, ha='center', va='center', color=color, fontsize=10)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

def plot_deletion_by_template_bar(del_by_template, summary, out_path):
    plt.figure(figsize=(7, 5))
    rates = del_by_template['DelRate'].reindex(BASES)
    plt.bar(BASES, rates.values, color=['#4C78A8', '#F58518', '#54A24B', '#E45756'])
    plt.xlabel('Template base')
    plt.ylabel('Deletion rate (weighted)')
    del_rate = summary.get('deletion_total_rate', np.nan)
    plt.title(f"Deletion Rate by Template Base (Total deletion rate: {del_rate:.6g})")
    for i, v in enumerate(rates.values):
        if pd.isna(v):
            continue
        plt.text(i, v, f"{v:.2e}", ha='center', va='bottom', fontsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

def plot_position_base_specific_deletion_lines(pos_base_del_rates, out_path):
    plt.figure(figsize=(15, 6))
    x = pos_base_del_rates.index.values
    colors = {'A': '#4C78A8', 'T': '#F58518', 'C': '#54A24B', 'G': '#E45756'}
    for b in BASES:
        if b not in pos_base_del_rates.columns:
            continue
        plt.plot(x, pos_base_del_rates[b].values, label=f'{b} deletion rate', linewidth=1.2, color=colors.get(b))
    plt.xlabel('Position (Pos)')
    plt.ylabel('Deletion rate (weighted)')
    plt.title('Per-position Deletion Rate by Template Base')
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend(ncol=4, fontsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

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
    for col in ['Pos', 'Muts', 'Depths', 'delcount', 'Acount', 'Tcount', 'Ccount', 'Gcount']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
    
    # 计算每个序列的长度
    sequence_lengths = df.groupby('Chrom')['Pos'].max().to_dict()
    print(f"序列长度范围: {min(sequence_lengths.values())} - {max(sequence_lengths.values())}")

    # ===== New: Global summary + substitution matrix + per-position base-specific deletions =====
    summary, deletion_by_template = compute_global_summary(df)
    sub_counts, sub_percent, total_sub_counts = compute_substitution_matrix(df)
    pos_base_del_rates = compute_position_base_specific_deletion(df)

    # Save new CSV outputs
    pd.DataFrame([summary]).to_csv(os.path.join(dirs['data_dir'], 'global_summary.csv'), index=False)
    deletion_by_template.to_csv(os.path.join(dirs['data_dir'], 'deletion_rate_by_template.csv'))
    sub_counts.to_csv(os.path.join(dirs['data_dir'], 'substitution_matrix_counts.csv'))
    sub_percent.to_csv(os.path.join(dirs['data_dir'], 'substitution_matrix_percent.csv'))
    pos_base_del_rates.to_csv(os.path.join(dirs['data_dir'], 'position_base_specific_deletion_rates.csv'))

    # New visualizations
    plot_substitution_matrix_heatmap(
        sub_percent,
        summary,
        os.path.join(dirs['images_dir'], 'substitution_matrix_heatmap.png'),
    )
    plot_deletion_by_template_bar(
        deletion_by_template,
        summary,
        os.path.join(dirs['images_dir'], 'deletion_rate_by_template_bar.png'),
    )
    plot_position_base_specific_deletion_lines(
        pos_base_del_rates,
        os.path.join(dirs['images_dir'], 'position_base_specific_deletion_rates.png'),
    )
    
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
    position_rates.to_csv(os.path.join(dirs['data_dir'], 'position_mutation_rates.csv'))
    
    # 可视化各种突变率沿序列位置的分布
    # 1. 总突变率分布
    plt.figure(figsize=(15, 6))
    plt.bar(position_rates.index, position_rates['MutationRate'])
    plt.xlabel('序列位置')
    plt.ylabel('平均突变率')
    plt.title('总突变率(ATCG)在不同序列位置的分布')
    plt.savefig(os.path.join(dirs['images_dir'], 'mutation_distribution.png'))
    plt.close()
    
    # 2. 各碱基突变率分布图 (toA/T/C/G)
    plt.figure(figsize=(15, 8))
    plt.plot(position_rates.index, position_rates['ARate'], 'r-', label='突变为A(toA)')
    plt.plot(position_rates.index, position_rates['TRate'], 'g-', label='突变为T(toT)')
    plt.plot(position_rates.index, position_rates['CRate'], 'b-', label='突变为C(toC)')
    plt.plot(position_rates.index, position_rates['GRate'], 'y-', label='突变为G(toG)')
    plt.xlabel('序列位置')
    plt.ylabel('碱基特异性突变率')
    plt.title('不同碱基在各序列位置的突变率(toA/T/C/G)')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(dirs['images_dir'], 'base_specific_mutation_rates.png'))
    plt.close()
    
    # 2.1 模板碱基突变率分布图 (Ato/Tto/Cto/Gto)
    plt.figure(figsize=(15, 8))
    plt.plot(position_rates.index, position_rates['Ato_MutRate'], 'r-', label='A位置突变率(Ato)')
    plt.plot(position_rates.index, position_rates['Tto_MutRate'], 'g-', label='T位置突变率(Tto)')
    plt.plot(position_rates.index, position_rates['Cto_MutRate'], 'b-', label='C位置突变率(Cto)')
    plt.plot(position_rates.index, position_rates['Gto_MutRate'], 'y-', label='G位置突变率(Gto)')
    plt.xlabel('序列位置')
    plt.ylabel('碱基特异性突变率')
    plt.title('不同模板碱基位置的突变率(Ato/Tto/Cto/Gto)')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(dirs['images_dir'], 'template_base_specific_mutation_rates.png'))
    plt.close()
    
    # 3. 总缺失率分布图
    plt.figure(figsize=(15, 6))
    plt.bar(position_rates.index, position_rates['DelRate'])
    plt.xlabel('序列位置')
    plt.ylabel('平均总缺失率')
    plt.title('总缺失率在不同序列位置的分布')
    plt.savefig(os.path.join(dirs['images_dir'], 'deletion_rate_distribution.png'))
    plt.close()
    
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
    plt.close()
    
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
    plt.close()
    
    # 3.3 模板碱基突变率与缺失率的对比
    plt.figure(figsize=(15, 12))
    # A的突变与缺失
    plt.subplot(2, 2, 1)
    plt.plot(position_rates.index, position_rates['Ato_MutRate'], 'r-', label='A位置突变率')
    plt.plot(position_rates.index, position_rates['A_DelRate'], 'r--', label='A位置缺失率')
    plt.xlabel('序列位置')
    plt.ylabel('错误率')
    plt.title('A位置的突变与缺失对比')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # T的突变与缺失
    plt.subplot(2, 2, 2)
    plt.plot(position_rates.index, position_rates['Tto_MutRate'], 'g-', label='T位置突变率')
    plt.plot(position_rates.index, position_rates['T_DelRate'], 'g--', label='T位置缺失率')
    plt.xlabel('序列位置')
    plt.ylabel('错误率')
    plt.title('T位置的突变与缺失对比')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # C的突变与缺失
    plt.subplot(2, 2, 3)
    plt.plot(position_rates.index, position_rates['Cto_MutRate'], 'b-', label='C位置突变率')
    plt.plot(position_rates.index, position_rates['C_DelRate'], 'b--', label='C位置缺失率')
    plt.xlabel('序列位置')
    plt.ylabel('错误率')
    plt.title('C位置的突变与缺失对比')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # G的突变与缺失
    plt.subplot(2, 2, 4)
    plt.plot(position_rates.index, position_rates['Gto_MutRate'], 'y-', label='G位置突变率')
    plt.plot(position_rates.index, position_rates['G_DelRate'], 'y--', label='G位置缺失率')
    plt.xlabel('序列位置')
    plt.ylabel('错误率')
    plt.title('G位置的突变与缺失对比')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(dirs['images_dir'], 'base_mut_vs_del_comparison.png'))
    plt.close()
    
    # 4. 突变+缺失总率分布图
    plt.figure(figsize=(15, 6))
    plt.bar(position_rates.index, position_rates['TotalMutDelRate'])
    plt.xlabel('序列位置')
    plt.ylabel('突变+缺失总率')
    plt.title('突变+缺失总率在不同序列位置的分布')
    plt.savefig(os.path.join(dirs['images_dir'], 'total_mutation_deletion_rate.png'))
    plt.close()
    
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
    plt.close()
    
    # 打印额外的分析信息
    print("\n位置特异性突变率统计已保存。")
    print("各种突变类型的位置分布图已生成。")
    
    '''
    # 位置归一化分析（已注释）
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
    
    # 归一化位置的突变率分布（已注释）
    # data_file = os.path.join(dirs['data_dir'], 'normalized_position_mutation_rates.csv')
    # norm_pos_rates.to_csv(data_file)
    '''
    # 生成HTML报告
    report_data = {
        'summary': summary,
        'deletion_by_template': deletion_by_template,
        'substitution_matrix_counts': sub_counts,
        'substitution_matrix_percent': sub_percent,
    }
    report_path = generate_html_report(file_path, dirs, report_data=report_data)
    
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

def generate_html_report(input_file_path, dirs, report_data=None, output_path=None):
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

    file_name = os.path.basename(input_file_path)
    today = datetime.datetime.now().strftime('%Y-%m-%d')

    # Optional report data (new global summary + matrices)
    summary = (report_data or {}).get('summary', {})
    deletion_by_template = (report_data or {}).get('deletion_by_template', None)
    sub_counts = (report_data or {}).get('substitution_matrix_counts', None)
    sub_percent = (report_data or {}).get('substitution_matrix_percent', None)

    def fmt(v, digits=6):
        try:
            if v is None or (isinstance(v, float) and np.isnan(v)):
                return ''
            return f"{float(v):.{digits}g}"
        except Exception:
            return str(v)

    summary_table = f"""
        <table>
          <tr><th>Metric</th><th>Value</th></tr>
          <tr><td>Total depths</td><td>{fmt(summary.get('total_depths'))}</td></tr>
          <tr><td>Total substitution muts (sum Muts)</td><td>{fmt(summary.get('total_substitution_muts'))}</td></tr>
          <tr><td>Total deletions (sum delcount)</td><td>{fmt(summary.get('total_deletions'))}</td></tr>
          <tr><td>Total substitution rate</td><td>{fmt(summary.get('substitution_total_rate'))}</td></tr>
          <tr><td>Total deletion rate</td><td>{fmt(summary.get('deletion_total_rate'))}</td></tr>
          <tr><td>Total error rate (sub+del)</td><td>{fmt(summary.get('total_error_rate'))}</td></tr>
        </table>
    """

    del_table_html = ""
    if isinstance(deletion_by_template, pd.DataFrame):
        del_df = deletion_by_template.copy()
        for c in ['delcount', 'Depths', 'DelRate']:
            if c in del_df.columns:
                del_df[c] = del_df[c].map(lambda x: '' if pd.isna(x) else f"{float(x):.6g}")
        del_table_html = del_df.to_html(classes='dataframe', border=0)

    sub_percent_table_html = ""
    if isinstance(sub_percent, pd.DataFrame):
        sp = sub_percent.copy()
        sp = sp.reindex(index=BASES, columns=BASES)
        for i in BASES:
            for j in BASES:
                if i == j:
                    sp.loc[i, j] = '—'
                else:
                    v = sp.loc[i, j]
                    sp.loc[i, j] = '' if pd.isna(v) else f"{float(v):.1f}%"
        sub_percent_table_html = sp.to_html(classes='dataframe', border=0, escape=False)

    html = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>突变分析报告 / Mutation Report</title>
  <style>
    :root {{
      --text: #1f2328;
      --muted: #57606a;
      --border: #d0d7de;
      --panel: #f6f8fa;
    }}
    body {{ font-family: Arial, sans-serif; margin: 18px; color: var(--text); }}
    .container {{ max-width: 1200px; margin: 0 auto; }}
    .header {{ display: flex; justify-content: space-between; align-items: baseline; gap: 16px; }}
    .title {{ font-size: 22px; font-weight: 700; }}
    .meta {{ font-size: 13px; color: var(--muted); line-height: 1.4; }}
    .toolbar {{ margin: 14px 0 6px; }}
    .print-button {{
      display: inline-block; padding: 10px 14px; background: #2da44e; color: white;
      border-radius: 6px; border: none; cursor: pointer; font-weight: 700;
    }}
    .section {{ margin: 22px 0; padding: 14px; background: var(--panel); border: 1px solid var(--border); border-radius: 10px; }}
    .section h2 {{ margin: 0 0 10px; font-size: 18px; }}
    details {{ margin: 12px 0; background: white; border: 1px solid var(--border); border-radius: 10px; overflow: hidden; }}
    summary {{
      list-style: none;
      cursor: pointer;
      padding: 10px 12px;
      font-weight: 700;
      background: linear-gradient(90deg, #fff8c5, #ffffff);
      border-bottom: 1px solid var(--border);
    }}
    summary::-webkit-details-marker {{ display: none; }}
    .details-body {{ padding: 12px; }}
    .figure {{ margin: 10px 0 2px; text-align: center; }}
    .figure img {{ max-width: 100%; border: 1px solid var(--border); border-radius: 10px; background: white; }}
    .caption {{ margin-top: 8px; font-size: 13px; color: var(--muted); }}
    table {{ width: 100%; border-collapse: collapse; background: white; }}
    th, td {{ border: 1px solid var(--border); padding: 8px 10px; font-size: 13px; }}
    th {{ background: #fff8c5; text-align: left; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    @media print {{
      .no-print {{ display: none; }}
      .section {{ break-inside: avoid; }}
      body {{ margin: 0; }}
    }}
  </style>
</head>
<body>
  <div class="container">
    <div class="header">
      <div class="title">突变分析报告 / Mutation Report</div>
      <div class="meta">
        Date: {today}<br>
        Input: {file_name}
      </div>
    </div>

    <div class="toolbar no-print">
      <button class="print-button" onclick="window.print()">Print to PDF</button>
    </div>

    <div class="section">
      <h2>1) Global Summary (weighted by depths)</h2>
      {summary_table}
      <p class="caption">
        CSV: <a href="../data/global_summary.csv">global_summary.csv</a>,
        <a href="../data/deletion_rate_by_template.csv">deletion_rate_by_template.csv</a>
      </p>
    </div>

    <details open>
      <summary>2) Global Substitution Spectrum</summary>
      <div class="details-body">
        <div class="figure">
          <img src="../images/substitution_matrix_heatmap.png" alt="substitution matrix heatmap">
          <div class="caption">Off-diagonal cells sum to 100%. Diagonal is not applicable.</div>
        </div>
        {sub_percent_table_html}
        <p class="caption">
          CSV: <a href="../data/substitution_matrix_counts.csv">substitution_matrix_counts.csv</a>,
          <a href="../data/substitution_matrix_percent.csv">substitution_matrix_percent.csv</a>
        </p>
      </div>
    </details>

    <details open>
      <summary>3) Global Deletion (by template base)</summary>
      <div class="details-body">
        <div class="figure">
          <img src="../images/deletion_rate_by_template_bar.png" alt="deletion rate by template base">
        </div>
        {del_table_html}
      </div>
    </details>

    <details open>
      <summary>4) Per-position Deletion by Template Base (no length normalization)</summary>
      <div class="details-body">
        <div class="figure">
          <img src="../images/position_base_specific_deletion_rates.png" alt="per-position base-specific deletion rates">
          <div class="caption">Each line is a weighted deletion rate among sites where Template==A/T/C/G at that Pos.</div>
        </div>
        <p class="caption">CSV: <a href="../data/position_base_specific_deletion_rates.csv">position_base_specific_deletion_rates.csv</a></p>
      </div>
    </details>

    <details>
      <summary>5) Additional Position-wise Plots (legacy)</summary>
      <div class="details-body">
        <div class="figure">
          <img src="../images/mutation_distribution.png" alt="mutation distribution">
          <div class="caption">Average mutation rate by Pos.</div>
        </div>
        <div class="figure">
          <img src="../images/base_specific_mutation_rates.png" alt="to-base mutation rates">
        </div>
        <div class="figure">
          <img src="../images/template_base_specific_mutation_rates.png" alt="from-base mutation rates">
        </div>
        <div class="figure">
          <img src="../images/deletion_rate_distribution.png" alt="deletion distribution">
        </div>
        <p class="caption">CSV: <a href="../data/position_mutation_rates.csv">position_mutation_rates.csv</a></p>
      </div>
    </details>

    <p class="caption" style="text-align:center;">Generated by MutationsCounter</p>
  </div>
</body>
</html>
"""

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)

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
    # Find all matching files (support a few common naming variants)
    patterns = [
        os.path.join(folder_path, "*.mpileup.cns.filter.xls"),
        os.path.join(folder_path, "*.mpileup*.cns*.xls"),
    ]
    files = []
    for p in patterns:
        files.extend(glob.glob(p))
    files = sorted(set(files))
    
    if not files:
        print(f"ERROR: No mpileup CNS .xls files found in folder {folder_path}")
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
