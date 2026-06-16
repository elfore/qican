#!/bin/bash

# ==============================================================================
# 脚本名称: generate_run.sh
# 功能描述: 根据 batchID 自动提取 skID，创建运行目录，替换参数并校验数据完整性
# 特点: 支持在任意路径执行，模版路径已硬编码
# 使用方法: ./generate_run.sh <batchid>
# 示例: ./generate_run.sh 260508151426_B172_SKII15277-YY2-WBhot-260508151422
# ==============================================================================

# 0. 硬编码源模版路径
SOURCE_PE_DIR="/mnt/gpfs1/Users/yangjinxurong/projects/qican/YY2/cpScript"
SOURCE_SE_1_DIR="/mnt/gpfs1/Users/yangjinxurong/projects/qican/YY2/cpScriptSE/SE_1/"
SOURCE_SE_2_DIR="/mnt/gpfs1/Users/yangjinxurong/projects/qican/YY2/cpScriptSE/SE_2/"
OUTPUT_DIR="/mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/PE150/"
# OUTPUT_DIR="/mnt/gpfs1/Users/yangjinxurong/projects/qican/YY2"

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

echo "---------------------------------------"
echo "源模版路径: $SOURCE_PE_DIR"
echo "解析信息:"
echo "  BatchID: $BATCH_ID"
echo "  skID:    $SK_ID"
echo "---------------------------------------"

# 3. 在目的目录下创建运行目录
OUTPUT_DIR_PE="$OUTPUT_DIR/$SK_ID"
OUTPUT_DIR_SE_1="$OUTPUT_DIR/${SK_ID}_SE_1"
OUTPUT_DIR_SE_2="$OUTPUT_DIR/${SK_ID}_SE_2"

if [ -d "$OUTPUT_DIR_PE" ]; then
    echo "警告: 目录 $OUTPUT_DIR_PE 已存在，程序退出。"
    exit 1
else
    mkdir -p "$OUTPUT_DIR_PE"
    echo "已创建运行目录: $OUTPUT_DIR_PE"
fi

if [ -d "$OUTPUT_DIR_SE_1" ]; then
    echo "警告: 目录 $OUTPUT_DIR_SE_1 已存在，程序退出。"
    exit 1
else
    mkdir -p "$OUTPUT_DIR_SE_1"
    echo "已创建运行目录: $OUTPUT_DIR_SE_1"
fi

if [ -d "$OUTPUT_DIR_SE_2" ]; then
    echo "警告: 目录 $OUTPUT_DIR_SE_2 已存在，程序退出。"
    exit 1
else
    mkdir -p "$OUTPUT_DIR_SE_2"
    echo "已创建运行目录: $OUTPUT_DIR_SE_2"
fi

# 4. 整理文件结构 (从 SOURCE_PE_DIR 复制模版)
TEMPLATES=("newNational" "PA" "PB" "PD" "tumor_panel" "WGS")
echo "正在从源路径同步PE模版..."
for item in "${TEMPLATES[@]}"; do
    if [ -d "$SOURCE_PE_DIR/$item" ]; then
        cp -r "$SOURCE_PE_DIR/$item" "$OUTPUT_DIR_PE/"
        echo "  [已同步] $item"
    else
        echo "  [错误] 模版目录不存在: $SOURCE_PE_DIR/$item"
        exit 1
    fi
done

echo "正在从源路径同步SE_1模版..."
for item in "${TEMPLATES[@]}"; do
    if [ -d "$SOURCE_SE_1_DIR/$item" ]; then
        cp -r "$SOURCE_SE_1_DIR/$item" "$OUTPUT_DIR_SE_1/"
        echo "  [已同步] $item"
    else
        echo "  [错误] 模版目录不存在: $SOURCE_SE_1_DIR/$item"
    fi
done

echo "正在从源路径同步SE_2模版..."
for item in "${TEMPLATES[@]}"; do
    if [ -d "$SOURCE_SE_2_DIR/$item" ]; then
        cp -r "$SOURCE_SE_2_DIR/$item" "$OUTPUT_DIR_SE_2/"
        echo "  [已同步] $item"
    else
        echo "  [错误] 模版目录不存在: $SOURCE_SE_2_DIR/$item"
    fi
done


# 5. 执行内容替换
echo "正在执行参数替换 (skID -> $SK_ID, batchID -> $BATCH_ID)..."
# 递归查找新PE目录下所有文件进行替换 (排除 .nf 文件)
find "$OUTPUT_DIR_PE" -type f ! -name "*.nf" | while read -r file; do
    # 仅处理文本文件
    if file "$file" | grep -q "text"; then
        sed -i "s/batchID/$BATCH_ID/g" "$file"
        sed -i "s/skID/$SK_ID/g" "$file"

    fi
done

# 递归查找新SE_1目录下所有文件进行替换 (排除 .nf 文件)
find "$OUTPUT_DIR_SE_1" -type f ! -name "*.nf" | while read -r file; do
    # 仅处理文本文件
    if file "$file" | grep -q "text"; then
        sed -i "s/batchID/$BATCH_ID/g" "$file"
        sed -i "s/skID/$SK_ID/g" "$file"

    fi
done

# 递归查找新SE_2目录下所有文件进行替换 (排除 .nf 文件)
find "$OUTPUT_DIR_SE_2" -type f ! -name "*.nf" | while read -r file; do
    # 仅处理文本文件
    if file "$file" | grep -q "text"; then
        sed -i "s/batchID/$BATCH_ID/g" "$file"
        sed -i "s/skID/$SK_ID/g" "$file"

    fi
done

echo "参数替换完成（已跳过所有 .nf 文件）。"

# 6. 检查 fq.gz 文件是否存在
echo "---------------------------------------"
echo "正在检查 fq.gz 数据文件完整性..."

MISSING_COUNT=0
TOTAL_FQ=0

# 从替换后的PE配置文件中提取所有以 .fq.gz 结尾的绝对路径
FQ_PATHS=$(grep -roh "/[^[:space:],;]*\.fq\.gz" "$OUTPUT_DIR_PE" | sort -u)

if [ -z "$FQ_PATHS" ]; then
    echo "提示: 未在配置文件中检测到任何 .fq.gz 路径。"
    exit 1
    rm -rf $OUTPUT_DIR_PE
else
    while read -r path; do
        [ -z "$path" ] && continue
        TOTAL_FQ=$((TOTAL_FQ + 1))
        if [ ! -f "$path" ]; then
            # echo "[缺失] $path"
            MISSING_COUNT=$((MISSING_COUNT + 1))
        fi
    done <<< "$FQ_PATHS"

    if [ "$MISSING_COUNT" -eq 0 ]; then
        echo "检查通过: 成功验证 $TOTAL_FQ 个 fq.gz 文件，全部存在。"
    else
        echo "检查结果: 发现 $MISSING_COUNT 个文件缺失 (总计 $TOTAL_FQ 个路径)。"
        echo "请检查原始数据路径是否正确。"
        rm -rf $OUTPUT_DIR_PE
    fi
fi

# 从替换后的SE_1配置文件中提取所有以 .fq.gz 结尾的绝对路径
FQ_PATHS=$(grep -roh "/[^[:space:],;]*\.fq\.gz" "$OUTPUT_DIR_SE_1" | sort -u)

if [ -z "$FQ_PATHS" ]; then
    echo "提示: 未在配置文件中检测到任何 .fq.gz 路径。"
    exit 1
    rm -rf $OUTPUT_DIR_SE_1
else
    while read -r path; do
        [ -z "$path" ] && continue
        TOTAL_FQ=$((TOTAL_FQ + 1))
        if [ ! -f "$path" ]; then
            # echo "[缺失] $path"
            MISSING_COUNT=$((MISSING_COUNT + 1))
        fi
    done <<< "$FQ_PATHS"

    if [ "$MISSING_COUNT" -eq 0 ]; then
        echo "检查通过: 成功验证 $TOTAL_FQ 个 fq.gz 文件，全部存在。"
    else
        echo "检查结果: 发现 $MISSING_COUNT 个文件缺失 (总计 $TOTAL_FQ 个路径)。"
        echo "请检查原始数据路径是否正确。"
        rm -rf $OUTPUT_DIR_SE_1
    fi
fi

# 从替换后的SE_2配置文件中提取所有以 .fq.gz 结尾的绝对路径
FQ_PATHS=$(grep -roh "/[^[:space:],;]*\.fq\.gz" "$OUTPUT_DIR_SE_2" | sort -u)

if [ -z "$FQ_PATHS" ]; then
    echo "提示: 未在配置文件中检测到任何 .fq.gz 路径。"
    exit 1
    rm -rf $OUTPUT_DIR_SE_2
else
    while read -r path; do
        [ -z "$path" ] && continue
        TOTAL_FQ=$((TOTAL_FQ + 1))
        if [ ! -f "$path" ]; then
            # echo "[缺失] $path"
            MISSING_COUNT=$((MISSING_COUNT + 1))
        fi
    done <<< "$FQ_PATHS"

    if [ "$MISSING_COUNT" -eq 0 ]; then
        echo "检查通过: 成功验证 $TOTAL_FQ 个 fq.gz 文件，全部存在。"
    else
        echo "检查结果: 发现 $MISSING_COUNT 个文件缺失 (总计 $TOTAL_FQ 个路径)。"
        echo "请检查原始数据路径是否正确。"
        rm -rf $OUTPUT_DIR_SE_2
    fi
fi


echo "---------------------------------------"
echo "全部流程处理完毕！"
echo "结果目录: $OUTPUT_DIR_PE"
echo "---------------------------------------"
