function fig = drawSGParamOptResult(window_num, polynomial_idx, RMSE_Xmix, label)
% @desc:
% show the result figure of Savitzky-Golay parameter plot
% the scatter gram has polynomial's maximum order as x axis
% and RMSE(RMSE_R and RMSE_Xmix are possibly used) as y axis. 
%
% @input:
% window_num: window number candidates
% polynomial_idx: polynomial's maximum order candidates
% RMSE_Xmix: RMSE as y plot
% label: label cell, {'title', 'xlabel', 'ylabel'}
%
% @output:
% fig: figure object which is generated
fig = figure;
color = {'m', 'r', 'g', 'b', 'k'};% one line which implies with the same w
shape = {'+', 'o', '*', '.', 'x'};% one line which implies with the same p
Legends = cell(numel(window_num), 1);
Series = zeros(numel(window_num), 1);
for I=1:numel(window_num)
    for J=1:numel(polynomial_idx)
        idx = (I-1)*numel(window_num)+J;
        colorspec = [color{I} shape{J}];
        f = plot(polynomial_idx(J), RMSE_Xmix(idx), colorspec);
        set(f, 'LineWidth', 3);
        set(f, 'MarkerSize', 15);
        if J==1
            Legends{I} = sprintf('w=%d', I);
            Series(I) = f(1, 1);
        end
        hold on;
    end
end
fs = 30;% FontSize
if length(label) < 3 || nargin < 4, label = {'RMSE X_{mix} vs. polynomial order', 'polynomial order', 'RMSE X_{mix}'}; end;
title(label{1}, 'FontSize', fs);
xlabel(label{2}, 'FontSize', fs);
ylabel(label{3}, 'FontSize', fs);
legend(Series, Legends, 'FontSize', fs);
set(gca, 'FontSize', fs);

end