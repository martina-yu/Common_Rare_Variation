"""Build a file and task catalog without importing analyses or reading data rows.

Run from any directory: python scripts/maintenance/organize_project.py
Only metadata and documentation under docs/project_catalog are generated.
"""

from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.functions.paths import PROJECT_ROOT

import ast
import csv
import hashlib
import io
import json
import re

ROOT = PROJECT_ROOT
OUT = ROOT / 'docs/project_catalog'
SKILL_URL = 'https://github.com/fmschulz/omics-skills/tree/main/skills/bioinformatics-project'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save_csv(name, fields, rows):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    (OUT / name).write_text(stream.getvalue(), encoding='utf-8')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sources = sorted(
        p
        for p in ROOT.glob('scripts/**/*.py')
        if not any(part in {'checkpoints', '.ipynb_checkpoints', '__pycache__'} for part in p.parts)
    )
    tasks, source_rows, mentions = [], [], []
    for path in sources:
        relative = path.relative_to(ROOT).as_posix()
        tree = ast.parse(path.read_text(encoding='utf-8'))
        doc = ast.get_docstring(tree) or ''
        role = (
            'helper'
            if '/functions/' in relative or path.name in {'function.py', '__init__.py'}
            else (
                'catalog'
                if path.name == 'organize_project.py'
                else (
                    'optional'
                    if path.name.startswith(('10_', '90_'))
                    else 'analysis' if re.match(r'0[1-9]_', path.name) else 'audit_or_manuscript'
                )
            )
        )
        source_rows.append(
            dict(
                path=relative,
                role=role,
                description=doc.splitlines()[0] if doc else '',
                sha256=digest(path),
            )
        )
        defaults = set()
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == 'run_local_step'
            ):
                for keyword in node.keywords:
                    if keyword.arg == 'default_tasks':
                        defaults.update(ast.literal_eval(keyword.value))
        # Read the explicit sequential sections without importing the analysis.
        literal_tasks = next(
            (
                ast.literal_eval(node.value)
                for node in tree.body
                if isinstance(node, ast.Assign)
                and any(
                    isinstance(target, ast.Name) and target.id == 'TASKS' for target in node.targets
                )
            ),
            (),
        )
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == 'parse_tasks'
                and len(node.args) >= 3
            ):
                defaults.update(ast.literal_eval(node.args[2]))
                if not literal_tasks and isinstance(node.args[1], (ast.List, ast.Tuple)):
                    literal_tasks = ast.literal_eval(node.args[1])
        for task in literal_tasks:
            section = next(
                (
                    node.lineno
                    for node in ast.walk(tree)
                    if isinstance(node, ast.If)
                    and isinstance(node.test, ast.Compare)
                    and isinstance(node.test.left, ast.Name)
                    and node.test.left.id == 'selected_task'
                    and any(
                        isinstance(value, ast.Constant) and value.value == task
                        for value in node.test.comparators
                    )
                ),
                '',
            )
            if not section:
                # Single-section stages execute directly under the main guard.
                section = max((
                    node.lineno for node in tree.body
                    if isinstance(node, ast.If) and '__name__' in ast.unparse(node.test)
                ), default='')
            tasks.append(
                dict(
                    script=relative,
                    task=task,
                    function='sequential section',
                    line=section,
                    default=str(task in defaults).lower(),
                    command=f'python {relative} --task {task}',
                )
            )
        functions = {
            n.name: n.lineno
            for n in tree.body
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
        }
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'computations' for t in node.targets)
                and isinstance(node.value, ast.Dict)
            ):
                for key, value in zip(node.value.keys, node.value.values):
                    if isinstance(key, ast.Constant) and isinstance(value, ast.Name):
                        task, function = key.value, value.id
                        tasks.append(
                            dict(
                                script=relative,
                                task=task,
                                function=function,
                                line=functions.get(function, ''),
                                default=str(task in defaults).lower(),
                                command=f'python {relative} --task {task}',
                            )
                        )
        for node in ast.walk(tree):
            if (
                path.name != 'organize_project.py'
                and isinstance(node, ast.Constant)
                and isinstance(node.value, str)
            ):
                value = node.value
                if (
                    len(value) < 500
                    and '\n' not in value
                    and re.search(r'\.(csv|tsv|json|svg|png|jpg|gz|npz|docx|md|yaml|yml)$', value)
                ):
                    mentions.append((relative, node.lineno, value))

    workflows = []
    for path in sorted(ROOT.glob('workflow/*.smk')):
        content = path.read_text()
        for match in re.finditer(r'^rule\s+(\w+):\s*\n(.*?)(?=^rule\s|\Z)', content, re.M | re.S):
            rule, body = match.groups()
            sections = dict(
                re.findall(
                    r'^    (input|output|shell):\s*\n(.*?)(?=^    \w+:|\Z)', body, re.M | re.S
                )
            )
            command = ' '.join(re.findall(r"['\"]([^'\"\n]+)['\"]", sections.get('shell', '')))
            for kind in ('input', 'output'):
                for value in re.findall(
                    r"['\"]((?:scripts|config|data|results|manuscript)/[^'\"]+)['\"]",
                    sections.get(kind, ''),
                ):
                    workflows.append(
                        dict(
                            workflow=path.relative_to(ROOT).as_posix(),
                            rule=rule,
                            direction=kind,
                            path=value,
                            exists=(
                                str((ROOT / value).exists()).lower()
                                if '{' not in value
                                else 'pattern'
                            ),
                            command=command,
                        )
                    )

    result_rows, unresolved = [], []
    for path in sorted((ROOT / 'results').rglob('*')):
        if (
            not path.is_file()
            or path.is_symlink()
            or any(part.startswith('.') for part in path.relative_to(ROOT).parts)
        ):
            continue
        relative = path.relative_to(ROOT).as_posix()
        declared = [r for r in workflows if r['direction'] == 'output' and r['path'] == relative]
        refs = [
            f'{source}:{line}'
            for source, line, literal in mentions
            if literal == relative or literal == path.name or literal.endswith('/' + path.name)
        ]
        status = 'workflow_output' if declared else 'source_reference' if refs else 'unresolved'
        category = path.relative_to(ROOT / 'results').parts[0]
        if path.parent == ROOT / 'results':
            category = 'documentation'
            status = 'organization_documentation'
        result_rows.append(
            dict(
                path=relative,
                category=category,
                bytes=path.stat().st_size,
                sha256=digest(path),
                provenance_status=status,
                workflow_rules=';'.join(f"{r['workflow']}:{r['rule']}" for r in declared),
                commands=';'.join(sorted({r['command'] for r in declared})),
                source_references=';'.join(sorted(set(refs))),
            )
        )
        if status == 'unresolved':
            unresolved.append(relative)

    assert len({r['path'] for r in result_rows}) == len(result_rows)
    assert len({(r['script'], r['task']) for r in tasks}) == len(tasks)
    assert tasks, 'No analysis tasks found'
    assert all((ROOT / r['path']).is_file() for r in result_rows)
    assert all((ROOT / r['script']).is_file() and r['line'] for r in tasks)
    save_csv('scripts.csv', ['path', 'role', 'description', 'sha256'], source_rows)
    save_csv('tasks.csv', ['script', 'task', 'function', 'line', 'default', 'command'], tasks)
    save_csv(
        'results.csv',
        [
            'path',
            'category',
            'bytes',
            'sha256',
            'provenance_status',
            'workflow_rules',
            'commands',
            'source_references',
        ],
        result_rows,
    )
    save_csv(
        'workflow_paths.csv',
        ['workflow', 'rule', 'direction', 'path', 'exists', 'command'],
        workflows,
    )
    unresolved_text = '# 尚未定位来源的结果\n\n静态扫描没有找到明确的工作流输出声明或文件名引用。这不表示文件无用、失效或可以删除；动态拼接路径和历史归档可能补充来源。\n\n'
    unresolved_text += '\n'.join(f'- [{p}](../../{p})' for p in unresolved) + '\n'
    (OUT / 'UNRESOLVED.md').write_text(unresolved_text)
    counts = {
        s: sum(r['provenance_status'] == s for r in result_rows)
        for s in ('workflow_output', 'source_reference', 'unresolved', 'organization_documentation')
    }
    summary = dict(
        source_files=len(source_rows),
        tasks=len(tasks),
        results=len(result_rows),
        result_provenance=counts,
        missing_literal_workflow_inputs=[
            r for r in workflows if r['direction'] == 'input' and r['exists'] == 'false'
        ],
        method='AST task discovery; literal source references; explicit Snakemake output declarations',
        real_analysis_rerun=False,
    )
    (OUT / 'audit.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
    text = f'''# 脚本与结果导航

整理日期：2026-10-13。采用 [omics-skills / bioinformatics-project]({SKILL_URL}) 的项目布局与来源追溯规则，适配现有 AoU Workbench、Snakemake 和 YAML 环境。

## 从哪里开始

- [项目现状](../../SUMMARY.md)：分析分支、阅读顺序与复现限制。
- [脚本说明](../../scripts/README.md)和[函数索引](../../scripts/FUNCTION_INDEX.md)：分析入口及计算位置。
- [tasks.csv](tasks.csv)：从当前源码提取的 {len(tasks)} 个任务、函数行号和单任务命令。
- [scripts.csv](scripts.csv)：{len(source_rows)} 个活动脚本和辅助模块的用途与校验和。
- [results.csv](results.csv)：{len(result_rows)} 个结果文件的位置、类别、大小、SHA-256 和来源证据。
- [workflow_paths.csv](workflow_paths.csv)：Snakemake 规则中的显式输入输出和命令。
- [UNRESOLVED.md](UNRESOLVED.md)：尚未定位来源的文件。
- [audit.json](audit.json)：清单统计与当前缺失的显式工作流输入。

## 怎样判断来源

`workflow_output`（{counts['workflow_output']} 个）表示工作流明确声明该文件是输出；命令是工作流提供的复现入口，仍需满足输入与环境要求。

`source_reference`（{counts['source_reference']} 个）只表示活动源码中出现完整路径或文件名，可作为生成者或消费者的定位线索；不能证明该脚本生成了当前文件，也不能证明它的运行版本。

`unresolved`（{counts['unresolved']} 个）表示本次静态扫描未找到上述证据。动态路径、历史脚本、归档和手工转换仍需人工核对。不要根据编号猜测生成脚本，不要根据该状态删除结果。

`organization_documentation`（{counts['organization_documentation']} 个）是本次新增的结果目录说明和组织日志，单独记录，不计入分析结果的来源缺口。

SHA-256 描述整理时磁盘文件的内容，不证明科学结论正确或当前脚本已经重跑。目录中的结果编号通常来自合并前步骤，不能直接等同当前 `01`–`09` 编号。

## 更新索引

```bash
python scripts/maintenance/organize_project.py
```

此命令只读取脚本、工作流和结果文件并生成本目录的索引；不导入分析模块，不读取参与者数据表，也不执行数据提取、建模、作图或稿件装配。文件散列逐块读取。

原 `data/`、`results/table/`、`results/figures/`、`results/report/`、`manuscript/` 与 `archive/` 继续承担现有职责。项目规定优先于技能的通用布局：参与者级中间数据保留在本地 `data/`，不迁移到公开结果目录；已有环境文件保留，精确环境锁缺口记录在方法说明中。
'''
    (OUT / 'README.md').write_text(text)
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k != 'missing_literal_workflow_inputs'},
            ensure_ascii=False,
        )
    )


if __name__ == '__main__':
    main()
