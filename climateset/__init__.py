from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = APP_ROOT.parent
CONFIGS = PROJECT_ROOT / "configs"
DATA_DIR = PROJECT_ROOT / "data"
RAW_DATA = DATA_DIR / "raw"
PROCESSED_DATA = DATA_DIR / "processed"
LOAD_DATA = DATA_DIR / "load"
META_DATA = DATA_DIR / "meta"
SCRIPT_DIR = PROJECT_ROOT / "scripts"
