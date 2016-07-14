function [fig, x, y, coeff] = plotSplinedCurve(xoriginal,yoriginal,k)
%% spline plot
% https://github.com/entzel/Splines
% EXAMPLE:
% x = [4.2,1.8,0,1.5,2, 1.5, 1.8,3.59,3,3]
% y = [2.2,0,1.1,4, 3.7,4, 5.8, 5.9 ,4.6, 5.9]
% t = [0 1 2 3 4 5 6 7 8 9 10]
% [fig,u,v] = splineplot(t, x, 10)
% [fig,u,g] = splineplot(t, y, 10)
% plot(v,g)
%
% coeff = [a1 a2 a3]
% y = ystart + a1 x + a2 x^2 + a3 x^3
% at every interval

%% implementation
n = length(xoriginal);
len = n + (k-1) * (n-1);
coeff = common.smoothing.calculateSplineCoeff(xoriginal, yoriginal);
x = zeros(len, 1);
y = zeros(len, 1);
for i=1:n-1% i is the index of interval
    xs = linspace(xoriginal(i), xoriginal(i+1), k+1);
    dx = xs - xoriginal(i);
    % 3rd-order term
    ys = coeff(i,3) * dx;
    % 2nd-order term
    ys = (ys + coeff(i,2)) .* dx;
    % 1st-order term and start point
    ys = (ys + coeff(i,1)) .* dx + yoriginal(i);
    for j=1:k
        x((i-1)*k+j) = xs(j);
        y((i-1)*k+j) = ys(j);
    end
end
x(len) = xoriginal(end);
y(len) = yoriginal(end);
figure;
fig = plot(xoriginal, yoriginal, 'o', x, y);

