#!/usr/bin/env bash
# SGA (Simple GEO Analyzer) — 卸载脚本
# 用法: bash uninstall.sh
set -euo pipefail

SGA_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "========================================"
echo "  SGA (Simple GEO Analyzer) 卸载程序"
echo "========================================"
echo "  项目目录: ${SGA_DIR}"
echo ""

# 1. 移除 ~/.local/bin/SGA
SGA_BIN="${HOME}/.local/bin/SGA"
if [ -f "${SGA_BIN}" ]; then
    rm -f "${SGA_BIN}"
    echo "[✓] 已移除: ${SGA_BIN}"
else
    echo "[!] 未找到 ${SGA_BIN}，跳过。"
fi

# 2. 清理 profile 中的 PATH 标记块
MARKER="# >>> SGA bin path >>>"
for PROFILE in "${HOME}/.bash_profile" "${HOME}/.profile"; do
    if [ -f "${PROFILE}" ] && grep -q "${MARKER}" "${PROFILE}" 2>/dev/null; then
        sed -i "/^${MARKER}/,/^# <<< SGA bin path <<</d" "${PROFILE}"
        echo "[✓] 已从 ${PROFILE} 移除 SGA PATH 条目。"
    fi
done

# 3. 移除项目目录（可选择性保留 res / data / conf）
echo ""
if [ -t 0 ]; then
    read -p "是否移除整个项目目录？[y/N]: " -r REPLY
    if [[ "${REPLY}" =~ ^[Yy]$ ]]; then
        BACKUP_DIR="${HOME}/SGA_backup"
        KEPT_ANY=false

        for SUBDIR in conf data res; do
            SRC="${SGA_DIR}/${SUBDIR}"
            if [ -d "${SRC}" ]; then
                read -p "  保留 ${SUBDIR}/ 目录？[y/N]: " -r KEEP
                if [[ "${KEEP}" =~ ^[Yy]$ ]]; then
                    DST="${BACKUP_DIR}/${SUBDIR}"
                    mkdir -p "$(dirname "${DST}")"
                    cp -r "${SRC}" "${DST}"
                    echo "    [✓] 已备份到 ${DST}"
                    KEPT_ANY=true
                fi
            fi
        done

        cd "${HOME}"
        rm -rf "${SGA_DIR}"
        echo "[✓] 已移除项目目录: ${SGA_DIR}"

        if [ "${KEPT_ANY}" = true ]; then
            echo "[!] 保留内容已备份到: ${BACKUP_DIR}"
        fi
    else
        echo "[!] 保留项目目录: ${SGA_DIR}"
    fi
else
    echo "[!] 非交互模式，保留项目目录。如需移除: rm -rf ${SGA_DIR}"
fi

# 4. 可选移除 conda 环境
ENV_NAME="SGA"
if conda env list 2>/dev/null | grep -q "^${ENV_NAME} "; then
    if [ -t 0 ]; then
        echo ""
        read -p "是否移除 conda 环境 '${ENV_NAME}'？[y/N]: " -r REPLY
        if [[ "${REPLY}" =~ ^[Yy]$ ]]; then
            conda env remove -n "${ENV_NAME}" -y
            echo "[✓] conda 环境 '${ENV_NAME}' 已移除。"
        else
            echo "[!] 保留 conda 环境 '${ENV_NAME}'。如需手动移除: conda env remove -n ${ENV_NAME}"
        fi
    else
        echo "[!] 非交互模式，保留 conda 环境 '${ENV_NAME}'。如需移除: conda env remove -n ${ENV_NAME}"
    fi
fi

echo ""
echo "========================================"
echo "  卸载完成"
echo "========================================"
