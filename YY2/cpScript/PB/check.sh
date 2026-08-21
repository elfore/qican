#!/bin/bash

# 输入文件名称（请确保文件名正确）
INPUT_FILE="param.yaml"

# 设定触发报错的最小文件大小阈值（单位：字节 Byte）
# 10240 Byte = 10 KB。如果一个 fq.gz 小于 10KB，基本可以判定没有有效数据
MIN_SIZE=10240 

if [ ! -f "$INPUT_FILE" ]; then
    echo "错误: 找不到输入文件 $INPUT_FILE"
    exit 1
fi

echo "开始检查文件中的 FASTQ 路径及状态..."
echo "--------------------------------------"

# 初始化统计变量
total_count=0
missing_count=0
empty_count=0
toosmall_count=0

# 提取 yaml 中的 fq.gz 路径
for fq_path in $(grep ".fq.gz" "$INPUT_FILE" | sed 's/.*: //' | tr ',' ' '); do
    ((total_count++))

    # 1. 检查文件是否存在
    if [ ! -f "$fq_path" ]; then
        echo "[缺失] ❌ 文件不存在: $fq_path"
        ((missing_count++))
        continue
    fi

    # 获取文件实际大小（字节）
    # 使用 stat 兼容 Linux 和 macOS
    if [[ "$OSTYPE" == "darwin"* ]]; then
        file_size=$(stat -f%z "$fq_path")
    else
        file_size=$(stat -c%s "$fq_path")
    fi

    # 2. 检查文件是否有内容（对于 .gz 来说，至少要大于 20 字节才算有压缩内容）
    if [ "$file_size" -le 20 ]; then
        echo "[空文件] ❌ 文件无实际内容: $fq_path (大小: ${file_size}B)"
        ((empty_count++))
        continue
    fi

    # 3. 检查文件是否太小
    if [ "$file_size" -lt "$MIN_SIZE" ]; then
        echo "[过小] ⚠️ 文件太小可能会导致分析报错: $fq_path (大小: $((file_size / 1024)) KB)"
        ((toosmall_count++))
        continue
    fi

    # 如果通过所有检查，可以选择静默或输出 OK
    # echo "[正常] 🟢 $fq_path"
done

echo "--------------------------------------"
echo "检查完成！"
echo "总计检查路径数: $total_count"
echo "  - 缺失文件数: $missing_count"
echo "  - 空 文件 数: $empty_count"
echo "  - 过小文件数: $toosmall_count"

# 最终结论判断
if [ "$missing_count" -eq 0 ] && [ "$empty_count" -eq 0 ] && [ "$toosmall_count" -eq 0 ]; then
    echo "结果: 所有文件均正常，可以放心进行下游分析！"
    exit 0
else
    echo "结果: 发现异常文件，请在跑流程前及时排查！"
    exit 1
fi