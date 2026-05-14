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

# 5. 生成帮助文本（help.yaml → help.txt）
echo "[→] 正在生成帮助文本..."
conda run -n SGA python -c "
import yaml, os, sys
help_path = os.path.join('${SGA_DIR}', 'conf', 'help.yaml')
with open(help_path) as f:
    data = yaml.safe_load(f)
width = 72
lines = []
lines.append('=' * width)
lines.append('SGA (Simple GEO Analyzer) — 配置项参考')
lines.append('=' * width)
lines.append('')
lines.append(data.get('usage', '').strip())
lines.append('')
for section in data.get('sections', []):
    lines.append(f'[{section[\"title\"]}]')
    lines.append('-' * 48)
    for field in section.get('fields', []):
        key = field.get('key', '')
        ftype = field.get('type', '')
        default = field.get('default', '')
        desc = field.get('desc', '').strip()
        choices = field.get('choices')
        note = field.get('note')
        lines.append(f'  {key}')
        lines.append(f'    类型: {ftype}')
        lines.append(f'    默认: {default}')
        if choices:
            lines.append(f'    可选: {\", \".join(choices)}')
        for line in desc.split(chr(10)):
            lines.append(f'    {line.strip()}')
        if note:
            lines.append(f'    注意: {note}')
        lines.append('')
with open(os.path.join('${SGA_DIR}', 'conf', 'help.txt'), 'w') as f:
    f.write(chr(10).join(lines))
print('help.txt 已生成')
" 2>/dev/null || echo "[!] help.txt 生成失败（可通过 sga help 查看在线帮助）"
echo "[✓] 帮助文本已生成: conf/help.txt"

# 6. 注册 SGA 命令到 ~/.local/bin
BIN_DIR="${HOME}/.local/bin"
SGA_BIN="${BIN_DIR}/sga"
echo "[→] 正在注册 sga 命令到 ~/.local/bin..."
mkdir -p "${BIN_DIR}"
cat > "${SGA_BIN}" << 'SGAEOF'
#!/usr/bin/env bash
# SGA (Simple GEO Analyzer) — version / help / config 由 bash 直接处理
set -euo pipefail

SGA_DIR="SGA_DIR_PLACEHOLDER"
CONF="${SGA_DIR}/conf/config.yaml"
HELP_TXT="${SGA_DIR}/conf/help.txt"

_sga_version() {
    git -C "$SGA_DIR" describe --tags --abbrev=0 2>/dev/null || echo "未知"
}

case "${1:-}" in
    version|--version|-v)
        _sga_version
        exit 0
        ;;
    help|--help|-h)
        if [ -f "$HELP_TXT" ]; then
            cat "$HELP_TXT"
        else
            echo "help.txt 未找到，请重新运行 setup.sh。"
        fi
        exit 0
        ;;
    config)
        shift
        if [ $# -eq 0 ]; then
            printf "  %-28s = %s\n" "version" "$(_sga_version)"
            for key in gse_id tar_gene multi_gene analysis_mode immune_method data_dir debug storage process p_threshold strict_filter log2fc_threshold enrichment_source_mode enrichment_gene_sets control_label exp_label exp_type group_select_col group_memory_enabled group_memory_use plot_data_warning overwrite_figures organism; do
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
                if [ -n "${val:-}" ]; then
                    echo "$val"
                else
                    echo "配置键 \"${1}\" 不存在"
                    exit 1
                fi
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
chmod +x "${SGA_BIN}"
echo "[✓] sga 可执行脚本已写入: ${SGA_BIN}"

# 确保 ~/.local/bin 在 PATH 中
MARKER="# >>> SGA bin path >>>"
if ! echo "${PATH}" | tr ':' '\n' | grep -qFx "${BIN_DIR}"; then
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

# 7. 安装 TumorDecon（免疫浸润分析依赖）
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
echo "    sga help"
echo "    sga gse_id=GSE143318 analysis_mode=diff"
echo ""
