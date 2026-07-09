import sys
from pathlib import Path

# Add resources/home/ubuntu to sys.path so run_qc_flags can be imported directly
sys.path.insert(0, str(Path(__file__).parent.parent / 'resources' / 'home' / 'ubuntu'))
