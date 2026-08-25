#!/usr/bin/env python3

from __future__ import annotations

import os
import sys
from pathlib import Path


def main() -> None:
    script = Path(__file__).with_name("checkCollWin.py")
    argv = [sys.executable, str(script), "--no-geometry", *sys.argv[1:]]
    os.execv(sys.executable, argv)


if __name__ == "__main__":
    main()
