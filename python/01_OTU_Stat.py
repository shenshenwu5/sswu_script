# 假设矩阵文件名为"abundance_matrix.csv"，分类文件名为"sample_groups.csv"
# 矩阵文件中，行代表细菌，列代表样本，每个值代表该样本中该细菌的丰度
# 分类文件中，有两列，第一列是样本名，第二列是分组信息

import pandas as pd

# 这里需要根据您实际文件存储的路径来修改
abundance_matrix_path = '02_otutab.csv'
grouping_path = '01_16S_soil_202_group.csv'
output_path = '03_output_results.txt'

# 加载细菌丰度矩阵，假设文件中包含行名和列名
abundance_matrix = pd.read_csv(abundance_matrix_path, index_col=0)

# 加载样本分组信息，假设文件中包含行名和列名
grouping = pd.read_csv(grouping_path)
grouping.columns = ['Sample', 'Group']

# 创建一个字典来存储结果
results = {}

# 对每个分组进行循环
for group, group_data in grouping.groupby('Group'):
    # 获取当前分组的样本列表
    group_samples = group_data['Sample'].tolist()

    # 获取当前分组的样本中每个细菌的丰度
    group_abundances = abundance_matrix[group_samples]

    # 统计每个细菌在超过80%的样本中丰度大于0.1的情况
    count_bacteria = group_abundances.apply(lambda x: (x > 0.1).sum(), axis=1)
    selected_bacteria = count_bacteria[count_bacteria >= 0.8 * len(group_samples)]

    # 保存结果
    results[group] = selected_bacteria

# 将结果输出到文件
with open(output_path, 'w') as f:
    for group, bacteria_counts in results.items():
        f.write(f'Group: {group}\n')
        f.write('Bacteria with abundance > 0.1 in >= 80% of samples:\n')
        for bacterium, count in bacteria_counts.iteritems():
            f.write(f'{bacterium}: {count}\n')
        f.write('\n')
