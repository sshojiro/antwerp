classdef Tracker < handle
    properties(Constant)
        PROP = struct('R2', 'r2',...
            'RMSE', 'rmse',...
            'Y', 'y');
        PROPTYPE = struct('CALC', 'calc',...
            'CV', 'cv',...
            'TEST', 'test',...
            'VALID', 'valid');
    end
    properties(GetAccess=public, SetAccess=private)
        e;% matrix
        times; % matrix
        r2;% struct
        rmse;% struct
        y;% struct
    end
    methods
        %--- public methods below ---%
        function self = Tracker(estimator, options)
            %% Tracker
            % 
            % @input:
            % estimator: estimator to adjust(supposed to be a vector)
            % options: option
            %
            % Ex) The class is able to change the y-value length
            % opt.y.TEST = 15;
            % t = common.Tracker(1:10, opt); t.y
            % ans = 
            %      calc: [0 0 0 0 0 0 0 0 0 0]
            %        cv: [0 0 0 0 0 0 0 0 0 0]
            %      test: [15x10 double]
            %     valid: [0 0 0 0 0 0 0 0 0 0]

            if nargin < 2, options = struct(); end
            self.e = estimator;
            self.times = 1:length(self.e);
            for names = struct2cell(self.PROP)'
                name = upper(names{1});
                if ~isfield(options, self.PROP.(name))
                   self.init(self.PROP.(name));
                else
                   self.init(self.PROP.(name), options);
                end
            end
        end
        
        function setVal(self, prop, proptype, idx, val)
            %% setVal
            %
            % @desc:
            % set {val} into {idx}-th item in {proptype}
            prop_tmp = self.(prop).(proptype);
            prop_tmp(:, idx) = val;
            self.(prop).(proptype) = prop_tmp;
        end
    end
    methods(Access=private)
        %--- private methods below ---%
        function init(self, prop, opt)
            if nargin < 3, opt = struct(); end
            for names = struct2cell(self.PROPTYPE)'
                name = upper(names{1});
                if ~isfield(opt, (prop))
                    self.(prop).(self.PROPTYPE.(name)) = zeros(1, length(self.e));
                else
                    if isfield(opt.(prop), (name))
                        self.(prop).(self.PROPTYPE.(name)) = zeros(opt.(prop).(name), length(self.e));
                    else
                        self.(prop).(self.PROPTYPE.(name)) = zeros(1, length(self.e));
                    end
                end
            end
        end
    end
end