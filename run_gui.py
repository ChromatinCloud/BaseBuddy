#!/usr/bin/env python3
"""
BaseBuddy GUI Launcher
Handles both development and bundled app environments
"""

import os
import sys
from pathlib import Path

def setup_bundled_environment():
    """Set up environment when running as bundled app"""
    if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
        # Running as bundled app
        bundle_dir = Path(sys._MEIPASS)
        tools_dir = bundle_dir / 'tools'
        
        # Add tools directory to PATH
        os.environ['PATH'] = f"{tools_dir}{os.pathsep}{os.environ.get('PATH', '')}"
        
        # Set environment variables
        os.environ['BASEBUDDY_TOOLS_DIR'] = str(tools_dir)
        os.environ['BASEBUDDY_BUNDLED'] = '1'
        
        # Ensure tools are executable
        for tool in tools_dir.glob('*'):
            if tool.is_file():
                try:
                    tool.chmod(0o755)
                except:
                    pass
                    
        return True
    return False

def main():
    """Main entry point"""
    # Set up environment if bundled
    is_bundled = setup_bundled_environment()
    
    # Add src to path for imports
    if not is_bundled:
        src_path = Path(__file__).parent / 'src'
        if src_path.exists():
            sys.path.insert(0, str(src_path))
    
    # Import and run the GUI
    try:
        from basebuddy.gui.main_app import BaseBuddyGUI
        
        # Create and run the app
        app = BaseBuddyGUI()
        app.mainloop()
        
    except ImportError as e:
        print(f"Error importing BaseBuddy GUI: {e}")
        print("Make sure you're in the correct environment or the app is properly built.")
        sys.exit(1)
    except Exception as e:
        print(f"Error running BaseBuddy GUI: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()