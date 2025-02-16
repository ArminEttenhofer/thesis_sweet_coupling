#! /usr/bin/env python3
import sys
from mule.sdc.testing import checkConvergence

checkConvergence(checkOrder=True, generatePlot=False)

if len(sys.argv) <= 1:
    print("*"*80)
    print("Convergence tests successful")
    print("*"*80)