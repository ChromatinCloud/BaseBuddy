#!/usr/bin/env python3
"""
Convenience launcher for BaseBuddy GUI
The actual GUI script has been moved to scripts/run_gui.py
"""

import subprocess
import sys
from pathlib import Path

script_path = Path(__file__).parent / "scripts" / "run_gui.py"
subprocess.run([sys.executable, str(script_path)] + sys.argv[1:])