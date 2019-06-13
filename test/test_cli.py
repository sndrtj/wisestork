from click.testing import CliRunner
import pytest

from wiseguy.wiseguy import count_cli, gcc_cli, newref_cli, zscore_cli


@pytest.fixture()
def runner():
    return CliRunner()


def test_cli_count_help(runner):
    result = runner.invoke(count_cli, "--help")
    assert result.exit_code == 0


def test_cli_gc_correct_help(runner):
    result = runner.invoke(gcc_cli, "--help")
    assert result.exit_code == 0


def test_cli_newref_help(runner):
    result = runner.invoke(newref_cli, "--help")
    assert result.exit_code == 0


def test_cli_zscore_help(runner):
    result = runner.invoke(zscore_cli, "--help")
    assert result.exit_code == 0
