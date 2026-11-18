import os
import argparse
import pandas as pd
import gzip

def remove_duplicate_barcodes(folder_path, output_folder):
    """
    删除文件夹下所有文件中，在指定列中重复的行（包括跨文件和文件内的重复），
    并将结果保存到指定输出文件夹中，保持原文件顺序。
    
    参数：
        folder_path (str): 输入文件夹路径
        output_folder (str): 输出文件夹路径
    """
    # 检查并创建输出文件夹
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    unique_barcodes = set()  # 全局唯一的元素集合
    files = [f for f in os.listdir(folder_path) if f.endswith(".txt.gz")]
    column_name = 'cell_bc'

    for file in files:
        file_path = os.path.join(folder_path, file)
        
        # 读取文件
        with gzip.open(file_path, 'rt') as gz_file:
            df = pd.read_csv(gz_file, sep="\t")
        
        # 确保列名正确
        if column_name not in df.columns:
            print(f"文件 {file} 中未找到列 {column_name}，跳过处理。")
            continue

        df = df.drop_duplicates(subset=["cell_bc", "x_pos", "y_pos"], keep="first").copy()
        
        # 找到当前文件中重复的行
        file_duplicates = df[column_name].duplicated(keep=False)
        
        # 找到跨文件的重复值
        cross_file_duplicates = df[column_name].isin(unique_barcodes)
        
        # 合并条件：文件内重复或跨文件重复
        to_remove = file_duplicates | cross_file_duplicates
        
        # 删除重复的行
        df_cleaned = df[~to_remove].copy()
        
        # 更新全局唯一值集合
        unique_barcodes.update(df_cleaned[column_name])
        
        # 保留清理后仍存在的行的顺序
        df_cleaned = df_cleaned.sort_index()

        # 保存到输出文件夹
        output_path = os.path.join(output_folder, file)
        with gzip.open(output_path, 'wt') as gz_out:
            df_cleaned.to_csv(gz_out, sep="\t", index=False)
        
        print(f"文件 {file} 处理完成，结果保存为 {output_path}。")

if __name__ == "__main__":
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description="Remove duplicate elements across and within files.")
    parser.add_argument("input_folder", type=str, help="Path to the input folder containing .txt.gz files")
    parser.add_argument("output_folder", type=str, help="Path to the output folder for cleaned files")
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 调用主函数
    remove_duplicate_barcodes(args.input_folder, args.output_folder)
