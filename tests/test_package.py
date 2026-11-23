"""Installed package and both terminal entry points must use the public name."""
from importlib import metadata
from pathlib import Path
import os
import subprocess
import sys
import sysconfig
import tomllib

import scaffold
from scaffold.cli import main


ROOT = Path(__file__).resolve().parents[1]


def test_installed_distribution_version_and_console_entry_point():
    project = tomllib.loads((ROOT / 'pyproject.toml').read_text(encoding='utf-8'))['project']
    distribution = metadata.distribution('scaffold')
    assert project['name'] == distribution.metadata['Name'] == 'scaffold'
    assert scaffold.__version__ == project['version'] == distribution.version
    assert project['scripts'] == {'scaffold': 'scaffold.cli:main'}
    entries = [entry for entry in distribution.entry_points if entry.group == 'console_scripts']
    assert [(entry.name, entry.value) for entry in entries] == [('scaffold', 'scaffold.cli:main')]
    assert entries[0].load() is main


def _assert_version_command(command, tmp_path):
    # Run outside the checkout and do not inherit a source-tree PYTHONPATH.
    # A stale editable install or a missing Windows launcher must fail here.
    env = {key: value for key, value in os.environ.items() if key.upper() != 'PYTHONPATH'}
    result = subprocess.run(command, cwd=tmp_path, env=env, capture_output=True,
                            text=True, encoding='utf-8', timeout=30)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == f'scaffold {scaffold.__version__}'


def test_python_module_entry_point_works_outside_checkout(tmp_path):
    _assert_version_command([sys.executable, '-m', 'scaffold', '--version'], tmp_path)


def test_installed_console_script_works_outside_checkout(tmp_path):
    filename = 'scaffold.exe' if os.name == 'nt' else 'scaffold'
    executable = Path(sysconfig.get_path('scripts')) / filename
    assert executable.is_file(), f'Console script was not installed: {executable}'
    _assert_version_command([str(executable), '--version'], tmp_path)
