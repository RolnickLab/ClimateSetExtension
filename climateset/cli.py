import typer
from click import Context
from typer.core import TyperGroup

from climateset import CONFIGS
from climateset.download.downloader import download_from_config_file

MINIMAL_DATASET = CONFIGS / "minimal_dataset.yaml"


class OrderCommands(TyperGroup):
    def list_commands(self, ctx: Context):
        return list(self.commands)


app = typer.Typer(
    no_args_is_help=True,
    cls=OrderCommands,
)


@app.command(name="info", help="General info.")
def info():
    typer.echo("\nWelcome to the ClimateSet CLI!\n")
    typer.echo("More information to come!\n")


@app.command(name="download", help="Download ClimateSet data via configuration file.")
def download_command(config: str = MINIMAL_DATASET):
    download_from_config_file(config)


if __name__ == "__main__":
    app()
