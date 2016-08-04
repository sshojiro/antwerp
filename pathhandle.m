function pathhandle(mode)
%% pathhandle(mode)
% 
% @input
% mode:
%     'add': add directories
%     'rm': remove directories
if nargin == 0, mode = 'add'; end
DIRS.ON_PATH = {'libs/'
    'experiments/'
    'scripts/'}';
if strcmp(mode, 'add')
    for targetdir = DIRS.ON_PATH
        addpath(genpath(targetdir{1}));
    end
elseif strcmp(mode, 'rm')
    for targetdir = DIRS.ON_PATH
        rmpath(genpath(targetdir{1}));
    end
end