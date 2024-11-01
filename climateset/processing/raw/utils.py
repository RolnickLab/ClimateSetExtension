from pathlib import Path


def create_generic_output_path(output_dir: Path, path: str, file: str) -> Path:
    """Creates an output path and the necessary parent directories.
    'Generic' in the sense that file names are not adapted as in other
    create_output_path functions used in the resolution processers.
    Args:
        output_dir (Path): First part of the path for the output dir
        path (str): Can e.g. stem from os.walk - path of the current file.
        file (str): Can e.g. stem from os.walk - current file.
    Returns:
        Path: New path pointing where the output file can be stored.
    """
    topic_dir = file.split("_")[0]
    out_path = Path(output_dir / topic_dir / Path(path.split(f"{topic_dir}/")[-1]) / file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    return out_path
