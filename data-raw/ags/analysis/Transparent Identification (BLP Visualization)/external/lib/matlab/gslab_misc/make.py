#! /usr/bin/env python
#****************************************************
# GET LIBRARY
#****************************************************
import os
from gslab_make.get_externals import *
from gslab_make.make_log import *
from gslab_make.run_program import *
from gslab_make.dir_mod import *

#****************************************************
# MAKE.PY STARTS
#****************************************************
set_option(makelog = 'log/make.log', output_dir = './log', temp_dir = '')

clear_dirs('./log', './external')
start_make_logging()

# GET_EXTERNALS
get_externals('externals.txt', './external')

# ANALYSIS
run_matlab(program = 'test/run_all_tests', changedir = True)

end_make_logging()

raw_input('\n Press <Enter> to exit.')
