#!/usr/bin/env bash
# SGA (Simple GEO Analyzer) — 一键安装脚本
# 用法: bash setup.sh
set -euo pipefail

SGA_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "========================================"
echo "  SGA (Simple GEO Analyzer) 安装程序"
echo "========================================"
echo ""

# 1. 检查 conda
if ! command -v conda &> /dev/null; then
    echo "[错误] 未找到 conda，请先安装 Miniconda 或 Anaconda。"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi
echo "[✓] 检测到 conda: $(conda --version)"

# 2. 创建 conda 环境
ENV_NAME="SGA"
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "[!] conda 环境 '${ENV_NAME}' 已存在，跳过创建。"
else
    echo "[→] 正在创建 conda 环境 '${ENV_NAME}'..."
    conda env create -f "${SGA_DIR}/environment.yml"
    echo "[✓] conda 环境创建完成。"
fi

# 3. 生成配置文件
CONFIG_FILE="${SGA_DIR}/conf/config.yaml"
CONFIG_TEMPLATE="${SGA_DIR}/conf/config.yaml.template"
if [ -f "${CONFIG_FILE}" ]; then
    echo "[!] 配置文件已存在，跳过生成。"
else
    echo "[→] 正在从模板生成配置文件..."
    cp "${CONFIG_TEMPLATE}" "${CONFIG_FILE}"
    echo "[✓] 配置文件已生成: conf/config.yaml"
    echo ""
    echo "  ┌──────────────────────────────────────────┐"
    echo "  │  请编辑 conf/config.yaml 设置分析参数：   │"
    echo "  │  - tar_gene: 目标基因                     │"
    echo "  │  - gse_id: GEO 数据集 ID                  │"
    echo "  │  - analysis_mode: 分析模式                │"
    echo "  │  执行 SGA help 查看完整配置说明           │"
    echo "  └──────────────────────────────────────────┘"
    echo ""
fi

# 4. 恢复卸载时保留的内容
BACKUP_DIR="${HOME}/SGA_backup"
if [ -d "${BACKUP_DIR}" ]; then
    echo "[→] 检测到 SGA 备份目录: ${BACKUP_DIR}"
    if [ -t 0 ]; then
        for SUBDIR in conf data res; do
            SRC="${BACKUP_DIR}/${SUBDIR}"
            if [ -d "${SRC}" ]; then
                read -p "  是否恢复 ${SUBDIR}/ 目录？[y/N]: " -r RESTORE
                if [[ "${RESTORE}" =~ ^[Yy]$ ]]; then
                    cp -r "${SRC}" "${SGA_DIR}/${SUBDIR}"
                    echo "    [✓] 已恢复 ${SUBDIR}/"
                fi
            fi
        done
        echo "[✓] 备份恢复完成。"
    else
        echo "[!] 非交互模式，跳过恢复。手动恢复: cp -r ${BACKUP_DIR}/* ${SGA_DIR}/"
    fi
fi

# 5. 注册 SGA 命令到 ~/.local/bin
BIN_DIR="${HOME}/.local/bin"
SGA_BIN="${BIN_DIR}/SGA"
if [ -x "${SGA_BIN}" ]; then
    echo "[!] SGA 命令已注册 (${SGA_BIN})，跳过。"
else
    echo "[→] 正在注册 SGA 命令到 ~/.local/bin..."
    mkdir -p "${BIN_DIR}"
    cat > "${SGA_BIN}" << EOF
#!/usr/bin/env bash
# SGA (Simple GEO Analyzer) command
conda run -n SGA python ${SGA_DIR}/main.py "\$@"
EOF
    chmod +x "${SGA_BIN}"
    echo "[✓] SGA 可执行脚本已写入: ${SGA_BIN}"

    # 确保 ~/.local/bin 在 PATH 中
    MARKER="# >>> SGA bin path >>>"
    if ! echo "${PATH}" | tr ':' '\n' | grep -qFx "${BIN_DIR}"; then
        # 优先写 ~/.bash_profile（如果存在），否则写 ~/.profile
        if [ -f "${HOME}/.bash_profile" ]; then
            PROFILE="${HOME}/.bash_profile"
        else
            PROFILE="${HOME}/.profile"
        fi
        if ! grep -q "${MARKER}" "${PROFILE}" 2>/dev/null; then
            cat >> "${PROFILE}" << EOP

${MARKER}
export PATH="${BIN_DIR}:\${PATH}"
# <<< SGA bin path <<<
EOP
            echo "[✓] PATH 已追加到 ${PROFILE}"
        fi
    else
        echo "[✓] ~/.local/bin 已在 PATH 中，无需修改。"
    fi
fi

# 6. 安装 TumorDecon（免疫浸润分析依赖）
echo "[→] 安装 TumorDecon（免疫浸润分析依赖）..."
conda run -n SGA pip install TumorDecon 2>/dev/null || {
    echo "[!] TumorDecon 安装失败（可能需要手动安装），忽略。"
}

echo ""
echo "========================================"
echo "  安装完成！"
echo "========================================"
echo ""
echo "  如果 ~/.local/bin 新加入 PATH，请执行："
echo "    source ~/.profile"
echo "  （或重新打开终端）"
echo ""
echo "  使用前请先编辑配置文件："
echo "    vim conf/config.yaml"
echo ""
echo "  然后即可使用："
echo "    SGA help"
echo "    SGA gse_id=GSE143318 analysis_mode=diff"
echo ""
