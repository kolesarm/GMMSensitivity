function test_compute_ahat
     
    addpath(genpath('../m'));
    addpath(genpath('../external'));
    addpath(genpath('../depend'));

    xTol        = 10^(-4); 
    x0          = [1; 0.1];
    nparam      = length(x0);
    
    g = testfunc(x0);
    nmom = length(g);
    
    weight      = eye(nmom);
    jacob_func  = @(param) NumJacob(@testfunc, param, xTol);
    
    ahat_numeric = compute_ahat(jacob_func, x0, xTol, weight, g);

    % The value returned by ahat_numeric has the correct dimensions
    assertEqual(size(ahat_numeric), [nparam, nparam]);
    
    % compute_ahat returns approximately the same result as analytic
    % derivation and basic matlab arithmetic.
    G1_analytic = [2*x0(2),            2*x0(1); ...
                   2*x0(2)/x0(1)^3,    -1/x0(1)^2; ...
                   exp(x0(1) + x0(2)), exp(x0(1) + x0(2))];
                   
    G2_analytic = [2*x0(1),            0; ...
                   -1/x0(1)^2,         0; ...
                   exp(x0(1) + x0(2)), exp(x0(1) + x0(2))];

    ahat_analytic = [G1_analytic' * weight * g, G2_analytic' * weight * g];
    
    
    assertElementsAlmostEqual(ahat_numeric, ahat_analytic, 'absolute', 10^-2);
    
    % The result of compute_ahat is symmetric
    assertElementsAlmostEqual(ahat_numeric, ahat_numeric', 'absolute', 10^-4);
end

function f = testfunc(param)
    f_1 = param(1)*param(1)*param(2);
    f_2 = param(2) / param(1);
    f_3 = exp(param(1) + param(2));
    
    f = [f_1; f_2; f_3];
end


    
    
    
    