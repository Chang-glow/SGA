#!/usr/bin/env python3
"""从 help.yaml 生成 help.txt，供 sga help 命令使用。"""
import yaml
import os
import sys

def main():
    help_path = os.path.join(os.path.dirname(__file__), 'help.yaml')
    out_path = os.path.join(os.path.dirname(__file__), 'help.txt')

    with open(help_path) as f:
        data = yaml.safe_load(f)

    lines = []
    lines.append('=' * 72)
    lines.append('SGA (Simple GEO Analyzer) — 配置项参考')
    lines.append('=' * 72)
    lines.append('')
    lines.append(data.get('usage', '').strip())
    lines.append('')

    for section in data.get('sections', []):
        lines.append(f'[{section["title"]}]')
        lines.append('-' * 48)
        for field in section.get('fields', []):
            key = field.get('key', '')
            lines.append(f'  {key}')
            lines.append(f'    类型: {field.get("type", "")}')
            lines.append(f'    默认: {field.get("default", "")}')
            choices = field.get('choices')
            if choices:
                lines.append(f'    可选: {", ".join(choices)}')
            for line in field.get('desc', '').strip().split('\n'):
                lines.append(f'    {line.strip()}')
            note = field.get('note')
            if note:
                lines.append(f'    注意: {note}')
            lines.append('')

    with open(out_path, 'w') as f:
        f.write('\n'.join(lines))
    print('help.txt 已生成')

if __name__ == '__main__':
    main()
