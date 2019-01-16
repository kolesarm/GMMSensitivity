diary off
diary '../log/test.log'
clear
clear global

addpath(genpath('../external/matlab_xunit/src/'))
addpath(genpath('../m'))

runxunit ../test -verbose

diary off
exit
