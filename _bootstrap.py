"""Make ``src/`` importable when scripts are run from a fresh clone (no pip install needed)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
