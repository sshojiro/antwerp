%% sample_code
% sample

%% init
clear; rehash; close all;
pathhandle;
load_constants;
DirName = util.setResultDir(mfilename(env.builtin.fullpath));

%% input
fold = 5;
saveFiles = false;

data_used = data.xlab.data20160713;
load([pwd env.dir.data data_used.file.data], 'X', 'y');


%% model construction

[Zx, xmean, xstd] = zscore(X);
[Zy, ymean, ystd] = zscore(y);

b = inv(Zx' * Zx) * Zx' * Zy;
predy = @(x) x * b * ystd + ymean; 

fig = FG.plotWithDiagnalLine({y}, {predy(Zx)}, {'' 'y_{obs}' 'y_{pred}'}, ...
    {}, {'.'}, true, [-5 5]);
util.saveJpg(fig, DirName, 'yyplot', 0, true, saveFiles);

if saveFiles
    save([DirName env.division util.addPrefixTime('result')]);
end