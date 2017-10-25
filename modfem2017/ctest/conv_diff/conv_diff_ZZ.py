#!/usr/bin/env python

import sys
import os

sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../'))
from ModFEMtest import TestZZ

arg1=sys.argv[1]
os.chdir(arg1)

msg = TestZZ()
print(msg)

