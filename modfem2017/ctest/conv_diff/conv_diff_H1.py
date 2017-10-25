#!/usr/bin/env python

import sys
import os

sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../'))
from ModFEMtest import TestH1

arg1=sys.argv[1]
os.chdir(arg1)

msg = TestH1()
print(msg)

