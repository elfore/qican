#!/bin/bash

# 检查是否输入了 yaml 文件
if [ -z "$1" ]; then
    echo "使用方法: sh check_yaml_samples.sh <你的参数文件>"
    exit 1
fi

YAML_FILE=$1

missing_count=0
total_count=0

# 使用 < <(...) 结构，把文本处理的结果作为输入流喂给 while 循环
while read -r file; do
    [ -z "$file" ] && continue
    ((total_count++))
    if [ ! -f "$file" ]; then
        ((missing_count++))
        echo -e "\e[31m[ 缺失 ]\e[0m $file"
    fi
done < <(grep -E "\.fastq\.gz" "$YAML_FILE" | tr ';' '\n' | sed -e 's/^[ \t]*//' -e 's/.*: //')

echo "--------------------------------------------------"
echo "检查完毕！共检测了 $total_count 个文件。"
if [ "$missing_count" -eq 0 ]; then
    echo -e "\e[32m🎉 恭喜，所有测序文件完整存在！\e[0m"
else
    echo -e "\e[31m⚠️ 警告：发现 $missing_count 个文件缺失，请检查！\e[0m"
fi
echo "--------------------------------------------------"