#!/bin/bash

# 1. 定义源目录（你 tree 命令看到的那个路径）
SRC="/mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/PE150/SKII15266"

# 2. 定义目标目录（当前所在目录）
DEST=$(pwd)

echo "正在从源路径同步脚本和配置文件..."
echo "源路径: $SRC"
echo "目标路径: $DEST"
echo "--------------------------------------"

# 3. 使用 find 命令筛选并拷贝
# 逻辑：
# - maxdepth 2: 只匹配到 tree 显示的两层深度
# - ! -path "*/CN*" : 排除掉以 CN 开头的样本文件夹
# - ! -path "*/work*" : 排除 work 目录
# - ! -path "*/result*" : 排除结果目录
# - \( -name "*.sh" -o -name "*.yaml" ... \): 只匹配特定的脚本和配置文件

find "$SRC" -maxdepth 2 \
    ! -path "*/CN*" \
    ! -path "*/.nextflow*" \
    ! -path "*/shell/*" \
    ! -path "*/work*" \
    ! -path "*/log*" \
    ! -path "*/result*" \
    ! -path "*/QC*" \
    ! -path "*/FASTQ*" \
    ! -name "trace-*" \
    \( \
        -name "*.sh" \
        -o -name "*.yaml" \
        -o -name "*.nf" \
        -o -name "*.config" \
        -o -name "*.info" \
        -o -type d \
    \) -print0 | while IFS= read -r -d '' item; do
    
    # 计算相对路径
    rel_path=${item#$SRC/}
    
    # 如果是目录且不是源根目录本身，则创建
    if [ -d "$item" ]; then
        if [ "$item" != "$SRC" ]; then
            mkdir -p "$DEST/$rel_path"
        fi
    # 如果是文件，则拷贝
    elif [ -f "$item" ]; then
        cp "$item" "$DEST/$rel_path"
        echo "已拷贝: $rel_path"
    fi
done

echo "--------------------------------------"
echo "同步完成！"

# fix this: create path skID_SE_1

sed -i "s#260508150833_B162_SKII15275-YY2BRS2hot-260508150833#batchID#g" param.yaml
sed -i "s#SKII15275#skID#g" param.yaml 

sed -i "s#260429170625_B143_SKII15266-YY2-23-DN-260429170627#batchID#g" param.yaml
sed -i "s#SKII15266#skID#g" param.yaml 

sed -i "s#260429170625_B143_SKII15266-YY2-23-DN-260429170627#batchID#g" sample.info
sed -i "s#SKII15266#skID#g" sample.info
sed -i "s#SKII15266#skID#g" genemid.sh

sed -i "s#260429170714_B141_SKII15265-YY223-CTRL-260429170714#batchID#g" param.yaml
sed -i "s#SKII15265#skID#g" param.yaml 