#!/usr/bin/env python3
"""
Launch BaseBuddy GUI

This script sets up the Python path and launches the BaseBuddy GUI.
"""

import sys
from pathlib import Path

# Add src directory to Python path
src_dir = Path(__file__).parent / "src"
sys.path.insert(0, str(src_dir))

# Import and run the GUI
from basebuddy.gui.main_app import run_gui

if __name__ == "__main__":
    print("Launching BaseBuddy GUI...")
    print("Note: This requires customtkinter to be installed.")
    print("If not installed, run: pip install customtkinter")
    print("")
    
    try:
        run_gui()
    except ImportError as e:
        print(f"Error: Missing dependency - {e}")
        print("\nTo install required GUI dependencies:")
        print("  pip install customtkinter")
    except Exception as e:
        print(f"Error launching GUI: {e}")
        import traceback
        traceback.print_exc()