# Documentation: scripts/README.md.
import argparse
from contextlib import contextmanager
import inspect
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]


@contextmanager
def component_arguments(name, arguments):
    previous = sys.argv
    sys.argv = [name, *arguments]
    try:
        yield
    finally:
        sys.argv = previous


def run_local_step(description, computations, default_tasks):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--task', action='append', choices=list(computations))
    parser.add_argument(
        '--plan', action='store_true', help='Show computations without executing them.'
    )
    parser.add_argument('--list-tasks', action='store_true')
    parser.add_argument('--where', nargs='?', const='*', metavar='TASK')
    parser.add_argument('component_args', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    names = args.task if args.task is not None else default_tasks
    if args.list_tasks or args.where is not None:
        names = list(computations) if args.where in (None, '*') else [args.where]
        if any(name not in computations for name in names):
            parser.error('Unknown computation: ' + str(args.where))
    forwarded = (
        args.component_args[1:] if args.component_args[:1] == ['--'] else args.component_args
    )
    if forwarded and len(names) != 1:
        parser.error('Component arguments require exactly one --task.')
    if not names:
        parser.error('Choose an explicit --task, or use --list-tasks.')
    for name in names:
        operation = computations[name]
        if args.plan or args.list_tasks or args.where is not None:
            path = Path(inspect.getsourcefile(operation)).resolve().relative_to(ROOT)
            print(
                f'{name}: {path}:{inspect.getsourcelines(operation)[1]} -> {operation.__name__}()'
            )
        else:
            print(f'Running {name}: {operation.__name__}()', flush=True)
            with component_arguments(operation.__name__, forwarded):
                operation()


def parse_tasks(description, tasks, defaults):
    """Select explicit sequential sections; metadata commands never run analyses."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--task', action='append', choices=tasks)
    parser.add_argument('--plan', action='store_true')
    parser.add_argument('--list-tasks', action='store_true')
    parser.add_argument('--where', nargs='?', const='*')
    parser.add_argument('component_args', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    selected = args.task if args.task is not None else list(defaults)
    if args.list_tasks or args.where is not None:
        selected = list(tasks) if args.where in (None, '*') else [args.where]
        if any(name not in tasks for name in selected):
            parser.error('Unknown task')
    forwarded = (
        args.component_args[1:] if args.component_args[:1] == ['--'] else args.component_args
    )
    if forwarded and len(selected) != 1:
        parser.error('Component arguments require exactly one --task.')
    if not selected:
        parser.error('Choose an explicit --task, or use --list-tasks.')
    if args.plan or args.list_tasks or args.where is not None:
        print('\n'.join(selected))
        raise SystemExit(0)
    sys.argv = [sys.argv[0], *forwarded]
    return selected
