#!/bin/bash

# 输入文件名称（请确保文件名正确）
INPUT_FILE="param.yaml"

if [ ! -f "$INPUT_FILE" ]; then
    echo "错误: 找不到输入文件 $INPUT_FILE"
    exit 1
fi

echo "开始检查文件中的 FASTQ 路径..."
echo "--------------------------------------"

# 初始化统计变量
total_count=0
missing_count=0

# 1. 使用 grep 提取所有包含 .fq.gz 的行
# 2. 使用 tr 将逗号替换为空格，以便处理一行中有 R1,R2 的情况
# 3. 循环检查每一个路径
for fq_path in $(grep ".fq.gz" "$INPUT_FILE" | sed 's/.*: //' | tr ',' ' '); do
    ((total_count++))
    
    if [ -f "$fq_path" ]; then
        # 如果文件存在，可以选择不输出，或者输出 OK
        # echo "[OK] $fq_path"
        :
    else
        echo "[缺失] $fq_path"
        ((missing_count++))
    fi
done

echo "--------------------------------------"
echo "检查完成！"
echo "总计检查路径数: $total_count"

if [ "$missing_count" -eq 0 ]; then
    echo "结果: 所有文件均存在。"
else
    echo "结果: 共有 $missing_count 个文件不存在，请检查路径！"
fi