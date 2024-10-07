from typer.testing import CliRunner

from climateset.cli import app

runner = CliRunner()


def test_app():
    result = runner.invoke(app, ["download", "--help"])
    assert result.exit_code == 0
    assert "Download ClimateSet data via configuration file." in result.stdout
