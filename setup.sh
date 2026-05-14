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
echo "[√] 检测到 conda: $(conda --version)"

# 2. 创建/更新 conda 环境
ENV_NAME="SGA"
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "[!] conda 环境 '${ENV_NAME}' 已存在，更新依赖..."
    conda env update -f "${SGA_DIR}/environment.yml" --prune
else
    echo "[→] 正在创建 conda 环境 '${ENV_NAME}'..."
    conda env create -f "${SGA_DIR}/environment.yml"
fi
echo "[√] conda 环境就绪。"

# 3. 生成配置文件
CONFIG_FILE="${SGA_DIR}/conf/config.yaml"
CONFIG_TEMPLATE="${SGA_DIR}/conf/config.yaml.template"
if [ -f "${CONFIG_FILE}" ]; then
    echo "[!] 配置文件已存在，跳过生成。"
else
    echo "[→] 正在从模板生成配置文件..."
    cp "${CONFIG_TEMPLATE}" "${CONFIG_FILE}"
    echo "[√] 配置文件已生成: conf/config.yaml"
    echo ""
    echo "  ┌──────────────────────────────────────────┐"
    echo "  │  请编辑 conf/config.yaml 设置分析参数：   │"
    echo "  │  - tar_gene: 目标基因                     │"
    echo "  │  - gse_id: GEO 数据集 ID                  │"
    echo "  │  - analysis_mode: 分析模式                │"
    echo "  │  执行 sga help 查看完整配置说明           │"
    echo "  └──────────────────────────────────────────┘"
    echo ""
fi

# 4. 恢复卸载时保留的内容
BACKUP_DIR="${HOME}/SGA_backup"
if [ -d "${BACKUP_DIR}" ]; then
    echo "[→] 检测到 SGA 备份目录: ${BACKUP_DIR}"
    if [ -t 0 ]; then
        read -p "  是否恢复备份？[y/N]: " -r RESTORE
        if [[ "${RESTORE}" =~ ^[Yy]$ ]]; then
            for SUBDIR in conf data res; do
                SRC="${BACKUP_DIR}/${SUBDIR}"
                [ -d "${SRC}" ] && cp -r "${SRC}" "${SGA_DIR}/"
            done
            echo "  [√] 备份已恢复。"
        fi
    else
        echo "[!] 非交互模式，跳过恢复。手动: cp -r ${BACKUP_DIR}/* ${SGA_DIR}/"
    fi
fi

# 5. 生成帮助文本（help.yaml → help.txt）
echo "[→] 正在生成帮助文本..."
conda run -n SGA python "${SGA_DIR}/conf/generate_help.py" 2>/dev/null || {
    echo "[!] help.txt 生成失败，跳过。"
}
echo "[√] 帮助文本已生成: conf/help.txt"

# 6. 注册 sga 命令到 ~/.local/bin
BIN_DIR="${HOME}/.local/bin"
SGA_BIN="${BIN_DIR}/sga"
echo "[→] 正在注册 sga 命令到 ~/.local/bin..."
mkdir -p "${BIN_DIR}"

CORE_KEYS="gse_id tar_gene multi_gene analysis_mode immune_method data_dir debug storage process p_threshold strict_filter log2fc_threshold enrichment_source_mode enrichment_gene_sets control_label exp_label exp_type group_select_col group_memory_enabled group_memory_use plot_data_warning overwrite_figures organism"

cat > "${SGA_BIN}" << 'SGAEOF'
#!/usr/bin/env bash
# sga (Simple GEO Analyzer) — version / help / config 由 bash 直接处理
set -euo pipefail

SGA_DIR="SGA_DIR_PLACEHOLDER"
CONF="${SGA_DIR}/conf/config.yaml"
HELP_TXT="${SGA_DIR}/conf/help.txt"
CORE_KEYS="CORE_KEYS_PLACEHOLDER"

_sga_version() {
    git -C "$SGA_DIR" describe --tags --abbrev=0 2>/dev/null || echo "未知"
}

case "${1:-}" in
    version|--version|-v)
        _sga_version
        exit 0
        ;;
    help|--help|-h)
        [ -f "$HELP_TXT" ] && cat "$HELP_TXT" || echo "help.txt 未找到，请重新运行 setup.sh。"
        exit 0
        ;;
    config)
        shift
        if [ $# -eq 0 ]; then
            printf "  %-28s = %s\n" "version" "$(_sga_version)"
            for key in $CORE_KEYS; do
                val=$(grep "^${key}:" "$CONF" 2>/dev/null | head -1 | sed "s/^${key}: *//")
                printf "  %-28s = %s\n" "$key" "${val:-<未设置>}"
            done
            echo ""
            echo "  用法: sga config [--all | <key>]"
            echo "  sga config --all    查看完整配置"
            echo "  sga config gse_id   查看单个配置项"
        elif [ "$1" = "--all" ]; then
            cat "$CONF"
        else
            if [ "$1" = "version" ]; then
                _sga_version
            else
                val=$(grep "^${1}:" "$CONF" 2>/dev/null | head -1 | sed "s/^${1}: *//")
                [ -n "${val:-}" ] && echo "$val" || { echo "配置键 \"${1}\" 不存在"; exit 1; }
            fi
        fi
        exit 0
        ;;
    *)
        conda run -n SGA python "${SGA_DIR}/main.py" "$@"
        ;;
esac
SGAEOF

sed -i "s|SGA_DIR_PLACEHOLDER|${SGA_DIR}|" "${SGA_BIN}"
sed -i "s|CORE_KEYS_PLACEHOLDER|${CORE_KEYS}|" "${SGA_BIN}"
chmod +x "${SGA_BIN}"
echo "[√] sga 可执行脚本已写入: ${SGA_BIN}"

# 7. 确保 ~/.local/bin 在 PATH 中
MARKER="# >>> SGA bin path >>>"
if ! echo "${PATH}" | tr ':' '\n' | grep -qFx "${BIN_DIR}"; then
    PROFILE="${HOME}/.bash_profile"
    [ -f "$PROFILE" ] || PROFILE="${HOME}/.profile"
    if ! grep -q "${MARKER}" "${PROFILE}" 2>/dev/null; then
        cat >> "${PROFILE}" << EOP

${MARKER}
export PATH="${BIN_DIR}:\${PATH}"
# <<< SGA bin path <<<
EOP
        echo "[√] PATH 已追加到 ${PROFILE}"
    fi
else
    echo "[√] ~/.local/bin 已在 PATH 中，无需修改。"
fi

echo ""
echo "========================================"
echo "  安装完成！"
echo "========================================"
echo ""
echo "  首次使用请先编辑配置文件："
echo "    vim conf/config.yaml"
echo ""
echo "  然后即可使用："
echo "    sga help"
echo "    sga gse_id=GSE143318 analysis_mode=diff"
echo ""
