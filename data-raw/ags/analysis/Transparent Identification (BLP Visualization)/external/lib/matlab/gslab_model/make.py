#! /usr/bin/env python
#****************************************************
# GET LIBRARY
#****************************************************
import os
import gslab_make

from gslab_make.get_externals import *
from gslab_make.make_log import *
from gslab_make.run_program import *
from gslab_make.dir_mod import *

#****************************************************
# MAKE.PY STARTS
#****************************************************

# SET DEFAULT OPTIONS
set_option(makelog = 'log/make.log', output_dir = './log', temp_dir = '')

clear_dirs('log/', 'external/', 'depend/')
start_make_logging()

# GET EXTERNALS AND DEPENDS
get_externals('externals.txt', './external')
get_externals('depends.txt', './depend')

# RUN ALL TESTS
run_matlab(program = 'test/run_all_tests.m')

end_make_logging()

raw_input('\n Press <Enter> to exit.')
