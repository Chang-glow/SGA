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

# 3. 注册 SGA 命令到 ~/.bashrc
BASHRC="${HOME}/.bashrc"
MARKER="# >>> SGA command >>>"
if grep -q "${MARKER}" "${BASHRC}" 2>/dev/null; then
    echo "[!] SGA 命令已注册，跳过。"
else
    echo "[→] 正在注册 SGA 命令到 ~/.bashrc..."
    cat << 'EOF' >> "${BASHRC}"

# >>> SGA command >>>
SGA() {
    conda run -n SGA python SGA_DIR_PLACEHOLDER/main.py "$@"
}
# <<< SGA command <<<
EOF
    # 替换路径占位符为实际路径
    sed -i "s|SGA_DIR_PLACEHOLDER|${SGA_DIR}|g" "${BASHRC}"
    echo "[✓] SGA 命令已注册。"
fi

# 4. 安装 TumorDecon（免疫浸润分析依赖）
echo "[→] 安装 TumorDecon（免疫浸润分析依赖）..."
conda run -n SGA pip install TumorDecon 2>/dev/null || {
    echo "[!] TumorDecon 安装失败（可能需要手动安装），忽略。"
}

echo ""
echo "========================================"
echo "  安装完成！"
echo "========================================"
echo ""
echo "  请执行以下命令使 SGA 命令生效："
echo "    source ~/.bashrc"
echo ""
echo "  然后即可使用："
echo "    SGA help"
echo "    SGA gse_id=GSE143318 analysis_mode=diff"
echo ""
