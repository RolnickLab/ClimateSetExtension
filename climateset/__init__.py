from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = APP_ROOT.parent
CONFIG = PROJECT_ROOT / "config"
DATA_DIR = PROJECT_ROOT / "data"
SCRIPT_DIR = PROJECT_ROOT / "scripts"
