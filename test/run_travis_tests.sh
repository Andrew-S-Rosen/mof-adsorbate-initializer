#!/bin/sh

# run_travis_tests.sh: Run the tests for Travis CI

python ../examples/add_O.py &&
python ../examples/add_H.py &&
python ../examples/add_CH4.py &&
python analyze_test.py