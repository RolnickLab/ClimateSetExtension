import typer

from climateset import CONFIGS
from climateset.download import download_from_config_file, downloader_config
from climateset.download.cmip6_downloader import CMIP6Downloader
from climateset.download.input4mips_downloader import Input4MipsDownloader

app = typer.Typer(no_args_is_help=True)

CONFIG_PATH = CONFIGS / "minimal_dataset.yaml"


@app.command(
    name="download-basic",
    help="Download ClimateSet data by building the config objects. See function content for more details.",
)
def basic_download():
    """
    By default, will download to the DATA_DIR folder. You can override this behavior modifying the config objects or by
    adding the `data_dir` key in the config file under each project.

    ex.
    CMIP6:
      models: [ "NorESM2-LM" ]
      variables: [ "tas" ]
      experiments: [ "historical", "ssp126" ]
      data_dir: "PATH/TO/DATA_DIR"
    input4MIPs:
      variables: [ "CO2", "CH4" ]
      experiments: [ "historical","ssp126" ]
      data_dir: "PATH/TO/DATA_DIR"
    """
    input4mips_config = downloader_config.create_input4mips_downloader_config_from_file(CONFIG_PATH)
    cmip6_config = downloader_config.create_cmip6_downloader_config_from_file(CONFIG_PATH)

    # If you want to specify where data will be downloaded, change the following:
    # input4mips_config.data_dir = "PATH_TO_DATA_DIR"
    # cmip6_config.data_dir = "PATH_TO_DATA_DIR"

    input4mips_downloader = Input4MipsDownloader(input4mips_config)
    input4mips_downloader.download()

    cmip6_downloader = CMIP6Downloader(cmip6_config)
    cmip6_downloader.download()


@app.command(
    name="download-from-config",
    help="Download ClimateSet data via download_from_config_file() function. See function content for more details.",
)
def alternative_approach():
    """
    By default, will download to the DATA_DIR folder. You can override this behavior by adding the `data_dir` key in the
    config file under each project.

    ex.
    CMIP6:
      models: [ "NorESM2-LM" ]
      variables: [ "tas" ]
      experiments: [ "historical", "ssp126" ]
      data_dir: "PATH/TO/DATA_DIR"
    input4MIPs:
      variables: [ "CO2", "CH4" ]
      experiments: [ "historical","ssp126" ]
      data_dir: "PATH/TO/DATA_DIR"
    """
    download_from_config_file(CONFIG_PATH)


if __name__ == "__main__":
    app()
