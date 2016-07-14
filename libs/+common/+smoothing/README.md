Splines
=======
use splinecoeff.m to get a vector of spline coefficients, it takes a vector of x values and a vector of y values
splineplot uses splinecoeff to plot a spline given a set of values. To get a cursive E plotted with cubic splines, enter
the following into the matlab command prompt:
x= [4.2,1.8,0,1.5,2, 1.5, 1.8,3.59,3,3]
y = [2.2,0,1.1,4, 3.7,4, 5.8, 5.9 ,4.6, 5.9]
t = [0 1 2 3 4 5 6 7 8 9 10]
[u,v] = splineplot(t,x,10)
[u,g] = splineplot(t, y,10)
plot(v,g)
