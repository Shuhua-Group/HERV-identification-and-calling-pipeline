def calculate_overlap_percentage(ref_start, ref_end, sv_start, sv_end):
    overlap_start = max(ref_start, sv_start)
    overlap_end = min(ref_end, sv_end)
    overlap_length = max(0, overlap_end - overlap_start)
    
    ref_length = ref_end - ref_start
    overlap_percentage = (overlap_length / ref_length) * 100
    
    return overlap_percentage
# 输入文件和输出文件的路径
from sys import argv
wd=argv[1]
input_file_path = wd+ '/merge_filter/reference_interval'
output_file_path = wd+ '/merge_filter/reference_interval_filtered_results.txt'


# 打开输入文件和输出文件
with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
    # 写入输出文件的表头，假设第一行是表头
    outfile.write(infile.readline())
    
    # 遍历文件的每一行
    for line in infile:
        fields = line.strip().split('\t')
        
        # 假设文件格式为：col1 col2 col3 col4 other_col1 other_col2 ...
        ref_start, ref_end, sv_start, sv_end = map(int, fields[:4])
        
        # 计算重合部分长度占 reference 区间长度的比例
        overlap_percentage = calculate_overlap_percentage(ref_start, ref_end, sv_start, sv_end)
        
        # 如果满足条件，写入输出文件
        if overlap_percentage > 80:
            outfile.write(line)

print(f"Filtered results have been written to {output_file_path}")

