# source /mnt/gpfs1/Users/yangjinxurong/software/miniconda3/bin/activate /mnt/gpfs1/Users/yangjinxurong/software/miniconda3/envs/nextflow/
source /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/activate base
nextflow -q -log log/log run /mnt/gpfs1/Users/yangjinxurong/projects/qican/YY1/bin/TB2/TB.nf -params-file param.yaml -c /mnt/gpfs1/Users/yangjinxurong/projects/qican/YY1/bin/TB2/config -bg -with-trace -resume
