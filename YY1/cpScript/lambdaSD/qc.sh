cd results
python /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY1/bin/lambdaSD/calculate_index_hopping_v2.py ../mapfile .
python /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY1/bin/lambdaSD/lambdaSD_matrix_modify.py misassigned_counts_matrix.xls
