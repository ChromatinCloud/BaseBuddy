#!/usr/bin/env python3
# This script is a stub. Please use 'basebuddy' CLI or 'python -m basebuddy.cli'.
import sys
print("Warning: basebuddy.py is deprecated. Use 'basebuddy' CLI or 'python -m basebuddy.cli'.", file=sys.stderr)
try:
    from basebuddy.cli import main
    if __name__ == "__main__":
        main()
except ImportError:
    print("Error: Could not import basebuddy.cli. Ensure BaseBuddy is installed correctly.", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"An unexpected error occurred: {e}", file=sys.stderr)
    sys.exit(1)
