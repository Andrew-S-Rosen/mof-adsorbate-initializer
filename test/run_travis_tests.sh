#!/bin/sh

# run_travis_tests.sh: Run the tests for Travis CI

python adsorbate_addition/examples/add_O.py &&
python adsorbate_addition/examples/add_H.py &&
python adsorbate_addition/examples/add_CH4.py &&
python analyze_test.py