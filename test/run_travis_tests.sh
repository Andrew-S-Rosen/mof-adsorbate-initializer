#!/bin/sh

# run_travis_tests.sh: Run the tests for Travis CI

python examples/add_O.py &&
python examples/add_O_OMS_omd.py &&
python examples/add_O_OMS_zeo.py &&
python examples/add_CH4_PEG_ASCII.py &&
python test/analyze_test.py