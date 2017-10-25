#!/usr/bin/env python


import sys
import subprocess
import os
import glob
import shutil
import zipfile

sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../'))
from ModFEMtest import TestDiff

arg1=sys.argv[1]
os.chdir(arg1)

msg = TestDiff()
print(msg)

