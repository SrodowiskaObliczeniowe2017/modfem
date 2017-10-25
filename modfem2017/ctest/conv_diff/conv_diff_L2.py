#!/usr/bin/env python

import sys
import os

sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../'))
from ModFEMtest import TestL2

arg1=sys.argv[1]
os.chdir(arg1)

msg = TestL2()
print(msg)

