function testBlpGetMisspec
%
% Testing the BlpGetMisspec function using fake data
%
addpath(genpath(fullfile(fileparts(pwd), 'external')))
addpath(genpath(fullfile(fileparts(pwd), 'depend')))
addpath(genpath(fullfile(fileparts(pwd), 'm')))

delta = 1;
avg_price = 10;
avg_mc 	  = 5;
k_xi 	  = 20;

M_same  = [1 1 1 2 2 3 3 3]';
M_rival = [7 7 7 6 6 5 5 5]';

test_loop = @(demand_supply, same_rival)BlpGetMisspec(demand_supply, same_rival, delta,...
                                           avg_price, avg_mc, k_xi, ...
                                           M_same, M_rival);

misspec_perturb = test_loop('supply', 'same')
expected = [0 0 0 0 0 0 0 0 -2 -2 -2 -4 -4 -6 -6 -6]';
assertEqual(misspec_perturb, expected);

misspec_perturb = test_loop('supply', 'rival')
expected = [0 0 0 0 0 0 0 0 -14 -14 -14 -12 -12 -10 -10 -10]';
assertEqual(misspec_perturb, expected);

misspec_perturb = test_loop('demand', 'same')
expected = [0.5 0.5 0.5 1 1 1.5 1.5 1.5 0 0 0 0 0 0 0 0]';
assertEqual(misspec_perturb, expected);

misspec_perturb = test_loop('demand', 'rival')
expected = [3.5 3.5 3.5 3 3 2.5 2.5 2.5 0 0 0 0 0 0 0 0]';
assertEqual(misspec_perturb, expected);

% Check that code handles case where dimensions don't match
M_rival = [7 7 7]';
test_loop = @(demand_supply, same_rival)BlpGetMisspec(demand_supply, same_rival, delta,...
                                           avg_price, avg_mc, k_xi, ...
                                           M_same, M_rival);
assertExceptionThrown(@()test_loop('demand', 'same'), '');


end