
#!/usr/bin/env python

"""Simple test that the imports can be done."""

import os
from pathlib import Path
import subprocess

scripts_pth=Path(Path(__file__).parent.parent.resolve() / "scripts").resolve()

print(scripts_pth)


subprocess.run(["python", scripts_pth / "ucs_batch_run.py", "--help"])
    
