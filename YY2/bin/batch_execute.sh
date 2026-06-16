#!/bin/bash

# ==============================================================================
# 脚本名称: batch_execute.sh
# 功能描述: 遍历指定目录下的子任务，校验 fq.gz/fa.gz 数据完整性并执行运行脚本
# 使用方法: ./batch_execute.sh <skid>
# 示例: ./batch_execute.sh SKII15277
# ==============================================================================

# 1. 参数校验
if [ -z "$1" ]; then
    echo "使用方法: sh $0 <目录名>"
    echo "示例: sh $0 SKII15277"
    exit 1
fi

# 2. 定义工作根目录和目标目录
# 基础路径保持不变
# BASE_PATH="/mnt/gpfs1/Users/yangjinxurong/projects/qican/YY2"
BASE_PATH="/mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/PE150/"
TARGET_DIR="$1"
WORK_DIR="$BASE_PATH/$TARGET_DIR"

if [ ! -d "$WORK_DIR" ]; then
    echo "目标目录 $WORK_DIR 不存在，不执行脚本"
    exit 1
fi

echo "开始处理工作目录: $WORK_DIR"
echo "========================================================================"

# 3. 遍历子任务文件夹
for sub_dir in "$WORK_DIR"/*; do
    # 仅处理目录
    if [ ! -d "$sub_dir" ]; then
        continue
    fi

    task_name=$(basename "$sub_dir")
    echo "[任务: $task_name] 正在检查数据状态..."

    # 提取该子任务所有配置文件中引用的数据文件路径 (.fq.gz 或 .fa.gz)
    DATA_FILES=$(grep -roh "/[^[:space:],;]*\.f[aq]\.gz" "$sub_dir" | sort -u)

    if [ -z "$DATA_FILES" ]; then
        echo "  - 提示: 未发现数据路径引用，默认视为数据就绪。"
        DATA_READY=true
    else
        DATA_READY=true
        while read -r file_path; do
            if [ -z "$file_path" ]; then continue; fi
            if [ ! -f "$file_path" ]; then
                echo "  - [缺失数据] $file_path"
                DATA_READY=false
            fi
        done <<< "$DATA_FILES"
    fi

    # 4. 如果数据完整，则尝试执行脚本
    if [ "$DATA_READY" = true ]; then
        echo "  - [状态] 数据完整，准备执行。"
        
        EXEC_SCRIPT=""
        # for s_name in "run2.sh" "do.sh"; do
        for s_name in "run.sh" "genemid.sh"; do
            if [ -f "$sub_dir/$s_name" ]; then
                EXEC_SCRIPT="$s_name"
                break
            fi
        done

        if [ -n "$EXEC_SCRIPT" ]; then
            echo "  - [启动] 正在进入目录执行 $EXEC_SCRIPT ..."
            (
                cd "$sub_dir" || exit
                nohup sh "$EXEC_SCRIPT" > "${task_name}_execution.log" 2>&1 &
                echo "  - [提交] 任务已在后台启动。日志: ${sub_dir}/${task_name}_execution.log"
            )
        else
            echo "  - [跳过] 未在目录下找到可执行脚本 (run.sh/do.sh)。"
        fi
    else
        echo "  - [跳过] 因关键数据文件缺失，不执行此任务。"
    fi
    echo "------------------------------------------------------------------------"
done

echo "所有任务预检与提交流程执行完毕。"
echo "可以使用 'ps -ef | grep .sh' 查看进程状态。"
