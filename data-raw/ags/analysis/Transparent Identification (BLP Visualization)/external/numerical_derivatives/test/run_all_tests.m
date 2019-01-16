function run_all_tests

echo on
diary '../log/test.log'

addpath(genpath('../external/matlab_xunit'))
addpath(genpath('../external/gslab_misc'))
addpath(genpath('../m'))

runxunit ../test

exit
