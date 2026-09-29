perl /mnt/gpfs1/Users/wangning/project/qican/mNGS/QC_stat.pl CN002341-S01-D01-L02-UDI1 >>stat.txt
perl /mnt/gpfs1/Users/wangning/project/qican/mNGS/QC_stat.pl CN002386-S01-D01-L01-UDI2 >>stat.txt 
perl /mnt/gpfs1/Users/wangning/project/qican/mNGS/QC_stat.pl CN002387-S01-D01-L01-UDI3 >>stat.txt
perl /mnt/gpfs1/Users/wangning/project/qican/mNGS/QC_stat.pl CN002388-S01-D01-L01-UDI4 >>stat.txt

perl /mnt/gpfs1/Users/wangning/project/qican/mNGS/get_zymo_result.pl CN002341-S01-D01-L02-UDI1 > CN002341-S01-D01-L02-UDI1.txt
perl /mnt/gpfs1/Users/wangning/project/qican/mNGS/classfiy.pl CN002386-S01-D01-L01-UDI2 > CN002386-S01-D01-L01-UDI2.txt
perl /mnt/gpfs1/Users/wangning/project/qican/mNGS/classfiy.pl CN002387-S01-D01-L01-UDI3 > CN002387-S01-D01-L01-UDI3.txt
perl /mnt/gpfs1/Users/wangning/project/qican/mNGS/classfiy.pl CN002388-S01-D01-L01-UDI4 > CN002388-S01-D01-L01-UDI4.txt

/mnt/gpfs/Users/fanlei/software/STRetch_failed/tools/miniconda/envs/STR/bin/python3 /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY1/bin/mNGS/export_feishu_abundance_history.py
/mnt/gpfs/Users/wangning/software/Miniconda/envs/qiime2-2019.7/bin/Rscript /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY1/bin/mNGS/prepare_abundance_boxplot_inputs.R
/mnt/gpfs/Users/wangning/software/Miniconda/envs/qiime2-2019.7/bin//Rscript /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY1/bin/mNGS/plot_abundance_boxplot_generic.R
