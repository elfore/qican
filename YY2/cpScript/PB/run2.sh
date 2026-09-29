source /mnt/gpfs1/Users/yangjinxurong/software/miniconda3/bin/activate /mnt/gpfs1/Users/yangjinxurong/software/miniconda3/envs/nextflow/
nextflow -q -log log/log run /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY2/cpScript/PB/antipsychotic_drugs_simple_v1.1.PE_dev_result2.nf -params-file param2.yaml -c /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY2/cpScript/PB/config -bg -with-trace -resume
