#!/usr/bin/env python

import sys
import os

sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../'))
from ModFEMtest import ModFEMtestRunner

bin = sys.argv[1]
arg1=sys.argv[2]
ppath=sys.argv[3]
arg3=sys.argv[4]
n_proc=0
config_name=""
if len(sys.argv) >= 6 :
    n_proc=sys.argv[5]
if len(sys.argv) == 7 :
    config_name=sys.argv[6]

runner = ModFEMtestRunner(bin,arg1,ppath,arg3,n_proc,config_name)
print(runner.Run())

