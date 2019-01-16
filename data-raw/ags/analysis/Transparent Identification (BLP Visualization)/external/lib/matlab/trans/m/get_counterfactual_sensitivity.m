%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GET_COUNTERFACTUAL_SENSITIVITY
%
% Takes as input:
% - Sensitivity matrix lambda (P by K)
% - Jacobian of transformation to apply to parameters (C by P)
%
% Optionally takes:
% - VCov of parameters (P by P)
% - VCov of moments (K by K)
%
% Returns as output:
% - Sensitivity of counterfactual to moments (C by K)
%
% Optionally returns if VCovs are provided:
% - Standardized sensitivity of counterfactual to moments via delta method (C by K)
% - VCov of counterfactual via delta method (C by C)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CfacLambda, CfacLambdaTilde, CfacVcov] = ...
get_counterfactual_sensitivity(Lambda, cfac_jacobian, param_vcov, moment_vcov)

    CfacLambda = cfac_jacobian * Lambda;

    if nargin == 4
        CfacVcov = cfac_jacobian * param_vcov * cfac_jacobian';
        se_cfac = sqrt(diag(CfacVcov));
        se_mom = sqrt(diag(moment_vcov));
        CfacLambdaTilde = get_standardized_sensitivity(CfacLambda, se_cfac, se_mom);
    else
        nargoutchk(0, 1);
    end
end
