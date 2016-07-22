%% generator.data
% 
% data generation


%% init
clear; close all; rehash;
pathhandle;
load_constants;

%% input

data_used = data.xlab.data20160713;

N = 100;
D = 33;

X = randn(N, D);
y = randn(N, 1) .^2;

mkdir([pwd env.dir.data data_used.dir.data]);

save([pwd env.dir.data data_used.file.data], 'X', 'y');