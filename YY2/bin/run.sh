#!/bin/bash

# ==============================================================================
# 脚本名称: generate_run.sh
# 功能描述: 根据 batchID 自动提取 skID，创建运行目录，替换参数并校验数据完整性
# 特点: 支持在任意路径执行，模版路径已硬编码
# ==============================================================================


# 1. 参数校验
if [ -z "$1" ]; then
    echo "使用方法: sh $0 <batchID>"
    echo "示例: sh $0 260508151426_B172_SKII15277-YY2-WBhot-260508151422"
    exit 1
fi

BATCH_ID=$1

# 2. 自动识别 skID
# 逻辑: 260508151426_B172_SKII15277-YY2-WBhot-260508151422 -> 第三个下划线后，第一个横杠前
SK_ID=$(echo "$BATCH_ID" | cut -d'_' -f3 | cut -d'-' -f1)

if [ -z "$SK_ID" ]; then
    echo "错误: 无法从 batchID ($BATCH_ID) 中解析出 skID。"
    exit 1
fi

sh /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY2/bin/generate_run.sh $BATCH_ID
sh /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY2/bin/batch_execute.sh $SK_ID

# 已经创建 crontab -e
# */10 * * * * sh /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY2/bin/robot.sh >> /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY2/monitor.log 2>&1