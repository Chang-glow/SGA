#!/usr/bin/env python3
"""从 help.yaml 生成 help.txt，供 sga help 命令使用。"""
import yaml
import os
import sys

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    help_path = os.path.join(script_dir, 'help.yaml')
    out_path = os.path.join(project_dir, 'docs', 'help.md')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(help_path) as f:
        data = yaml.safe_load(f)

    lines = []
    lines.append('# SGA (Simple GEO Analyzer) — 配置项参考')
    lines.append('')
    lines.append(data.get('usage', '').strip())
    lines.append('')

    for section in data.get('sections', []):
        lines.append(f'## {section["title"]}')
        lines.append('')
        lines.append('| 配置项 | 类型 | 默认值 | 说明 |')
        lines.append('|--------|------|--------|------|')
        for field in section.get('fields', []):
            key = field.get('key', '')
            ftype = field.get('type', '')
            fdefault = str(field.get('default', '')).replace('|', '\\|')
            desc_lines = field.get('desc', '').strip().split('\n')
            desc = ' '.join(s.strip() for s in desc_lines)
            choices = field.get('choices')
            if choices:
                desc += f' 可选: {", ".join(choices)}'
            note = field.get('note')
            if note:
                desc += f' 注意: {note}'
            lines.append(f'| `{key}` | {ftype} | `{fdefault}` | {desc} |')
        lines.append('')

    with open(out_path, 'w') as f:
        f.write('\n'.join(lines))
    print('help.txt 已生成')

if __name__ == '__main__':
    main()
