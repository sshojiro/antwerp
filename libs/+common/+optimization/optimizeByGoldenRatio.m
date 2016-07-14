function optimum = optimizeByGoldenRatio(lb, ub, func, option)
%% optimizeByGoldenRatio
% 
% @desc: 
% optimization by golden ratio,
% which is one-variable optimization problem
% 
% @input:
% lb: lower bound
% ub: upper bound
% func: objective fujnction
% 
% @example:
% lb = -100;
% ub = 100;
% objfun = @(x) x.^3 + x;
% opt = common.optimization.optimizeByGoldenRatio(lb, ub, objfun);
% %=> returns -100
% 
% lb = -1;
% ub = 1;
% objfun = @(x) x.^2;
% opt = common.optimization.optimizeByGoldenRatio(lb, ub, objfun);
% %=> returns zero

%% set constants and defaults
TAU.PLUS = (1+sqrt(5))/2;
if nargin < 4, option = struct(); end

%% initialize
tmp.lb = lb;
tmp.ub = ub;
notconverge = true;
if ~isfield(option, 'epsilon')
    epsilon = 0.000000001;% default
else
    epsilon = option.epsilon;
end

%% optimization process
while notconverge
    new.lb = tmp.lb + ((TAU.PLUS - 1)/TAU.PLUS) * (tmp.ub - tmp.lb);
    new.ub = tmp.lb + (1/TAU.PLUS) * (tmp.ub - tmp.lb);
    %--- update lb, ub ---%
    fl = func(new.lb);
    fu = func(new.ub);
    if fl >= fu,
        tmp.lb = new.lb;
    else
        tmp.ub = new.ub;
    end
    if isfield(option, 'disp'), disp(tmp); end
    if (tmp.ub - tmp.lb) < epsilon, notconverge = false; end
end
optimum.input = mean([tmp.ub tmp.lb]);
optimum.output = func(optimum.input);
end