function Xnonnegative = cutNegativeToZero( X )
%% Xnonnegative = cutNegativeToZero( X )
% @desc: cut off values below zero
% return +1 -> +1
%        -1 -> 0
%
% ex)
% X
% X =
%     0.8252   -1.0582   -0.2725
%     1.3790   -0.4686    1.0984
% common.cutNegativeToZero(X)
% ans =
%     0.8252         0         0
%     1.3790         0    1.0984

    negative = find(X<0);
    X(negative) = zeros(size(negative));
    Xnonnegative = X;
end

