function [func, ynew] = generateSplinedCurve(params)
%% generateSplinedCurve

%% implementor area
% TODO: implement params.coeff case
import common.smoothing.*;
n = length(params.xtrain);
func = cell(1, n-1);
if ~isfield(params, 'coeff')
    coeff = calculateSplineCoeff(params.xtrain, params.ytrain);
else
    coeff = params.coeff;
end

for i = 1:n-1
    % for each interval
    func{i} = @(x) coeff(i, 3) .* x.^3 + coeff(i, 2) .* x.^2 + coeff(i, 1) .* x + params.ytrain(i);
end
if ~isfield(params, 'xnew')
    % if params.xnew does not exist, return spline functions as cell, func{}
    params.xnew = [];
    ynew = [];
else
    % else, return ynew as spline functions' output
    ynew = zeros(size(params.xnew));
    len = length(params.xnew);
    leno= length(params.xtrain);
    for i = 1:len
        for j = 1:leno-1
            if params.xnew(i) - params.xtrain(j) >= 0 && params.xnew(i) - params.xtrain(j+1) <= 0
                ynew(i) = func{j}(params.xnew(i) - params.xtrain(j));
                break;
            end
        end
    end
end