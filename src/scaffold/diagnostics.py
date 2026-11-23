"""Read-only checks for a local conversion installation."""
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import sys


def doctor(moldir):
    packages = {}
    for name in ('scaffold', 'anarcii', 'biopython', 'gemmi', 'PyYAML', 'torch', 'boltzgen'):
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = None
    from .schemes import PRESETS
    path = Path(moldir).resolve()
    try:
        from boltzgen.data.parse.schema import YamlDesignParser
        parser_error = None
    except Exception as exc:
        parser_error = str(exc)
    return {'passed': all(packages.values()) and path.is_file() and parser_error is None,
            'python': sys.version, 'python_executable': sys.executable, 'packages': packages,
            'moldir': str(path), 'moldir_exists': path.is_file(), 'official_parser_import_error': parser_error,
            'schemes': PRESETS, 'note': 'Read-only import/file checks; does not run numbering or parse a scaffold.'}
