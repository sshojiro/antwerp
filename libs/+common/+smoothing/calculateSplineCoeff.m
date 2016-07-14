%% calculateSplineCoeff
% Calculates coefficients of cubic spline
% Input: x,y vectors of data points ie ([1,2,4,7],[3,4,1,2])
%   plus two optional extra data v1, vn
% Output: matrix of coefficients b1,c1,d1;b2,c2,d2;...
function coeff = calculateSplineCoeff(x,y)
n=length(x);
A=zeros(n,n);            % matrix A is nxn determined by number of x vars
b=zeros(n,1);            %b in Ax =b  
for i=1:n-1              % define the deltas for S2
    dx(i) = x(i+1)-x(i);
    dy(i)=y(i+1)-y(i);
end
for i=2:n-1              % load the A matrix
    A(i,i-1:i+1)=[dx(i-1) 2*(dx(i-1)+dx(i)) dx(i)];
    b(i)=3*(dy(i)/dx(i) - dy(i-1)/dx(i-1)); % right-hand side
end
% Set endpoint conditions

A(1,1) = 1;              % natural spline conditions
A(n,n) = 1;
coeff=zeros(n,3);
coeff(:,2)=A\b;          % solve for c coefficients

for i=1:n-1              % solve for b and d
    coeff(i,3)=(coeff(i+1,2)-coeff(i,2))/(3*dx(i));
    coeff(i,1)=dy(i)/dx(i)-dx(i)*(2*coeff(i,2)+coeff(i+1,2))/3;
end
coeff=coeff(1:n-1,1:3);
