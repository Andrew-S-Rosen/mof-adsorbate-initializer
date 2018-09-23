#!/bin/sh

# run_travis_tests.sh: Run the tests for Travis CI
cd examples
python add_O.py &&
python add_O_OMS_omd.py &&
python add_O_OMS_zeo.py &&
python add_CH4_PEG_ASCII.py &&
cd ../
python test/analyze_test.py