#!/usr/bin/env python3
"""
This module runs any/all of the 'broad' test for SOFT.
"""

import broadutil
import subprocess
import sys

# Tests
import distribution.runtest as distribution
import orbits.classification.runtest as orbitclassification

TESTS = [
    {'name': 'distribution', 'desc': 'Verify that distribution functions are properly applied', 'run': distribution.run},
    {'name': 'orbit-classification', 'desc': 'Verify that orbits are properly classified', 'run': orbitclassification.run}
]

SOFTPATH = '../../build/src/soft'

def broad_help():
    global TESTS
    print("runbroad.py [test1 [test2 [...]]]\n")

    print("Available tests:\n")

    for test in TESTS:
        print("  {0:<20}  {1}".format(test['name'], test['desc']))

    print("\nInvoke with 'all' to run all tests.")

def main(argv):
    broadutil.init()

    if len(argv) == 0:
        broad_help()

    for i in range(0, len(argv)):
        run_test(argv[i])

def run_test(testname):
    global TESTS

    if testname == 'all':
        for test in TESTS:
            test['run'](SOFTPATH)
    else:
        for test in TESTS:
            if testname == test['name']:
                test['run'](SOFTPATH)
                break

if __name__ == '__main__':
    main(sys.argv[1:])

