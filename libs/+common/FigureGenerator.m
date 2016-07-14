classdef FigureGenerator
	%% DESCRIPTION:
	%    Module for drawing figures
	%  EXAMPLE:(usually used with common.Utility module)
	% 	 util = common.Utility;
	%    FG = common.FigureGenerator;
	%    fig = FG.plotSeveralBandedSpectra(xaxis, Region, Spectra, {'title', 'xlabel', 'ylabel'}, 1:2, {'spectrum1', 'spectrum2'}, CustomLineStyle);
	% 	 util.saveJpg(fig, DirName, 'pure_spectra', 0, true, saveFiles);
	properties(Constant)
		Line = {'r-', 'm-', 'b-', 'g-', 'y-', 'c-', 'w-', 'k-'};
		util = common.Utility;
		DashedLine = {'--m', '--r', '--b', '--g', '--y'};
		FontSize = 'FontSize';
		FontSizeVal  = 30;% fontsize of letters on figure
		LineWidth= 'LineWidth';
		LineWidthVal = 2;% line width value
		MarkerLineWidth = 'LineWidth';
		MarkerSize = 'MarkerSize';
		MarkerLineWidthVal = 3;% marker line width value
		MarkerSizeVal = 15;% marker size value
		PlotStyle = {'xm', '*r', '+b'};
		IdxName = 'x';
	end
	methods
		function FG = FigureGenerator
			%% INITIALIZER
			if nargin ~= 0,
				sprintf('initialized');
			end
        end
        function fig = addFigToSubplot(FG, fig_target, coordinate, target, fig)
            %% fig = addFigToSubplot(FG, fig_target, coordinate, target, fig)
            %
            % @input:
            % fig: a figure object to add
            % sub: sub
            % 
            % @output:
            % figadded
            %
            % @reference:
            % http://www.mathworks.com/matlabcentral/answers/101806-how-can-i-insert-my-matlab-figure-fig-files-into-multiple-subplots
            %
            % @example:
            % load_constants;
            % ord = [2 1];
            % fig1 = figure;
            % plot(1:10, sin(1:10), '-');
            % fig2 = figure;
            % plot(1:0.1:10, cos(1:0.1:10), '-');
            % fig = figure;% this fig must be defined just before running the method.
            % fig = FG.addFigToSubplot( fig1, ord, 1, fig );
            % fig = FG.addFigToSubplot( fig2, ord, 2, fig );
            if nargin < 5, fig = figure; end
            sub = subplot(coordinate(1), coordinate(2), target);
            ax = gca(fig_target);
            ffig = get(ax, 'children');
            copyobj(ffig, sub);
        end
        
        function fig = plot3DScattergram(FG, Data, Labels, option, num_grids, lb, ub)
            %% plot3DScattergram(FG, Data, Labels, option, num_grid, lb, ub)
            % 
            % FG.plot3DScattergram(Data, Labels, option, Grid, lb, ub)
			%
            % @input:
            % Data: data that scatters
            % num_grid: number of grids
            % Labels: labels of the scattergram
            % option
            % option.tickinterval: tick label interval
            % option.hasgrid: show grid or no
            % option.legends
            % lb: lower bound
            % ub: upper bound
            % 
            % @example:
            % load_constants;
            % num_sample = 100;
            % num_dim = 3;
            % Data = {
            %     randn(num_sample, num_dim)
            %     randn(num_sample, num_dim)
            %     randn(num_sample, num_dim)
            %     randn(num_sample, num_dim)
            %     };
            % Labels = {'x' 'y' 'z'};
            % option.hasgrid = true;
            % option.legends = {
            %     'cluster1'
            %     'cluster2'
            %     'cluster3'
            %     'cluster4'
            %     };
            % fig = FG.plot3DScattergram(Data, Labels, option, 3);
            
            num_dim = 3;% 3d plot
            if nargin < 7,
                Dcat = cell2mat(Data);
                ub = max(Dcat);
            end
            if nargin < 7, lb = min(Dcat); end
            if nargin < 6,
                Grid = cell(1, num_dim);
                for i = 1:num_dim
                    Grid{i} = linspace(lb(i), ub(i), num_grids)';
                end
            end
            if nargin < 4, option = struct(); end
            fig = figure;
            co = get(gca,'ColorOrder');
            hold on;
            if isfield(option, 'tickinterval')
                tickinterval = option.tickinterval;
            else
                tickinterval = 2;
            end
            opt.isnew = true;
            p4legends = zeros(1, length(Data));
            for i = 1:length(Data)
                if isfield(option, 'color')
                    color_one_marker = option.color{i};
                else
                    color_one_marker = common.Utility.generateColorSpec(i, 'Marker');
                    % co(mod(i,size(co,1))+1,:);
                end
                if ~isempty(Data{i})
                    if isfield(option, 'markersize')
                        ptmp = plot3(Data{i}(:,1),Data{i}(:,2),Data{i}(:,3), ...
                            color_one_marker, ...
                            'MarkerSize', option.markersize);
                    else
                        ptmp = plot3(Data{i}(:,1),Data{i}(:,2),Data{i}(:,3), ...
                            color_one_marker,...
                            FG.MarkerSize, FG.MarkerSizeVal);
                    end
                    p4legends(i) = ptmp(1, 1);
                else
                    if isfield(option, 'legends')
                        option.legends(i,:) = [];
                    end
                end
            end
            p4legends(p4legends==0) = [];
            if isfield(option, 'hasgrid') && option.hasgrid
                num_grid = length(Grid{1});
                num_line = num_grid^2;
                x = repmat(Grid{1}', 2, num_grid);
                y = [repmat(lb(2), 1, num_line); repmat(ub(2), 1, num_line)];
                z = repmat(Grid{3}(common.cyclic(1:num_grid))', 2, 1);
                plot3(x,y,z,'-k');
                x = repmat(Grid{1}(common.cyclic(1:num_grid))', 2, 1);
                y = repmat(Grid{2}', 2, num_grid);
                z = [repmat(lb(3), 1, num_line); repmat(ub(3), 1, num_line)];
                plot3(x,y,z,'-k');
                x = [repmat(lb(1), 1, num_line); repmat(ub(1), 1, num_line)];
                y = repmat(Grid{2}(common.cyclic(1:num_grid), :)', 2, 1);
                z = repmat(Grid{3}', 2, num_grid);
                plot3(x,y,z,'-k');
            end
            zlim([lb(3) ub(3)]);
            if length(Labels) >= 4
                zlabel(Labels{4}, FG.FontSize, FG.FontSizeVal);
            end
            FG.addTickToFigure('z', Grid{3}(1:tickinterval:end)', opt);
            
            if isfield(option, 'legends')
                if iscell(option.legends)
                    legend(p4legends, option.legends, FG.FontSize, FG.FontSizeVal);
                else
                    error('   option.legends must be a cell that contains strings');
                end
            end
            xlim([lb(1) ub(1)]);
            ylim([lb(2) ub(2)]);
            if length(Labels) >= 3
                title(Labels{1}, FG.FontSize, FG.FontSizeVal);
                xlabel(Labels{2}, FG.FontSize, FG.FontSizeVal);
                ylabel(Labels{3}, FG.FontSize, FG.FontSizeVal);
            elseif length(Labels) == 2
                xlabel(Labels{1}, FG.FontSize, FG.FontSizeVal);
                ylabel(Labels{2}, FG.FontSize, FG.FontSizeVal);
            end
            FG.addTickToFigure('x', Grid{1}(1:tickinterval:end)', opt);
            FG.addTickToFigure('y', Grid{2}(1:tickinterval:end)', opt);
            box on; axis square;
            az = 18;
            el = 10;
            view(az, el);
		end

        function fig = generateColorMap(FG, Object, MapT, DivT, GridT, MinLO, MaxLO, Labels, option)
            %% generateColorMap(Object, MapT, DivT, GridT, MinLO, MaxLO, Labels, option)
            % fig = FG.generateColorMap(Object, MapT, DivT, GridT, MinLO, MaxLO, Labels, option)
            %
			% @desc: plot of chip densities for all samples
            %       generate heat map
			% @input:
			% Object: ChipDensAll
			% MapT:
			% DivT:
			% GridT:
			% @output:
			% fig: figure object to save
            if nargin < 9, option = struct(); end
            fig = figure;
            const = layout.const().ja;
            hold on;
            imagesc(MapT{1},MapT{2},reshape(Object,DivT{1},DivT{2})');
            colormap(hot)% make the map gray-scale
            if ~isfield(option, 'hascolorbar'), colorbar; end% set colorbar
            xlim([GridT{1}(1) GridT{1}(end)]);
            ylim([GridT{2}(1) GridT{2}(end)]);
            if length(Labels)==3
                title(Labels{1}, FG.FontSize, FG.FontSizeVal);
                xlabel(Labels{2}, FG.FontSize, FG.FontSizeVal);
                ylabel(Labels{3}, FG.FontSize, FG.FontSizeVal);
            elseif length(Labels)==2
                xlabel(Labels{1}, FG.FontSize, FG.FontSizeVal);
                ylabel(Labels{2}, FG.FontSize, FG.FontSizeVal);
            else
                xlab = [const.ordinal{1} const.principlecomp ' -' const.percent];%sprintf('‘æ1Žå¬•ª -“');
                ylab = [const.ordinal{2} const.principlecomp ' -' const.percent];%sprintf('‘æ2Žå¬•ª -“');
                xlabel(xlab, FG.FontSize, FG.FontSizeVal);
                ylabel(ylab, FG.FontSize, FG.FontSizeVal);
            end
            if isfield(option, 'tickinterval')
                tickinterval = option.tickinterval;
            else
                tickinterval = 2;
            end
            opt.isnew = true;
            FG.addTickToFigure('x', GridT{1}(1:tickinterval:end)', opt);
            FG.addTickToFigure('y', GridT{2}(1:tickinterval:end)', opt);
            box on;axis square;
            % draw white lines like grids
            if isfield(option, 'hasgrid') && option.hasgrid
                plot([MinLO(1);MaxLO(1)],[GridT{2}';GridT{2}'],'-w');
                plot([GridT{1}';GridT{1}'],[MinLO(2);MaxLO(2)],'-w');
            end
            Cmax = max(Object(~isinf(Object)));
            Cmin = min(Object(~isinf(Object)));
            Crange = max(abs([Cmax,Cmin]));
            if isfield(option, 'hassymmetricbar') && option.hassymmetricbar
                if ~isempty(Crange)
                    caxis([common.Utility.cutup(-Crange) Crange]);
                end
            else
                caxis([0 0.7]);% colorbar range
            end
            if isfield(option, 'colorrange')
                if length(option.colorrange) == 2
                    caxis([common.Utility.cutup(-option.colorrange(1)) option.colorrange(2)]);
                else
                    caxis([common.Utility.cutup(-option.colorrange(1)) option.colorrange(1)]);
                end
            end

            if isfield(option, 'mapcolor')
                colormap(option.mapcolor);
            end
        end

        function addTickToFigure(FG, direction, tk, options )
            %% addTickToFigure
            %ADDXTICK --- Add new tick to x-axis
            %   Ex.    common.FigureGenerator.addTickToFigure('x', [1:2:21]);
            %   Ex2.
            % opt.isnew = true; -> override XTickLabels All
            % opt.tklbl = 2/pi; -> put the defined XTickLabel below the
            % xaxis
            % common.FigureGenerator.addTickToFigure('x', [1:2:21], opt);
            % see http://ichiro-maruta.blogspot.jp/2009/01/matlab.html
            switch direction
                case 'x'
                    axisdirection = 'XTick';
                    ticklabel = 'XTickLabel';
                case 'X'
                    axisdirection = 'XTick';
                    ticklabel = 'XTickLabel';
                case 'y'
                    axisdirection = 'YTick';
                    ticklabel = 'YTickLabel';
                case 'Y'
                    axisdirection = 'YTick';
                    ticklabel = 'YTickLabel';
                case 'z'
                    axisdirection = 'ZTick';
                    ticklabel = 'ZTickLabel';
                case 'Z'
                    axisdirection = 'ZTick';
                    ticklabel = 'ZTickLabel';
                otherwise
                    axisdirection = 'XTick';
            end
            if size(tk, 1) == length(tk)
                tmp = tk';
                tk = tmp;
            end
            if nargin==2
                set(gca,axisdirection,unique(sort([get(gca,axisdirection), tk ] )), FG.FontSize, FG.FontSizeVal);
            else
                if isfield(options, 'isnew') && options.isnew
                    % clear and insert again
                    origin = tk;
                    [sorted, idx_sort] = sort(origin);
                    if isfield(options, 'tklbl'), ticklabels = options.tklbl; end
                else
                    origin = [get(gca,axisdirection), tk ];
                    [sorted, idx_sort] = sort(origin);
                    if isfield(options, 'tklbl'), ticklabels = [cellstr(get(gca, ticklabel))' options.tklbl]; end
                end
                [coord, idx_unique] = unique(sorted);
                
                if ~isfield(options, 'tklbl')
                    set(gca, axisdirection, coord, ...
                    FG.FontSize, FG.FontSizeVal);
                else
                    idx_sort = idx_sort(idx_unique);
                    set(gca, axisdirection, coord, FG.FontSize, FG.FontSizeVal);
                    if isfield(options, 'tklbl')
                        set(gca, ticklabel, ticklabels(idx_unique(idx_sort)), ...
                            FG.FontSize, FG.FontSizeVal);
                    end
                end
            end% end if
        end% end function

		function fig = plotFrequency(FG, Data, Labels, BinHandle)
			% @formula: fig = FG.plotFrequency(Data, Labels, BinHandle)
			% @desc:  draw frequency plot(histogram)
			% @input: (nxd matrix)Data: unsupervised data
			% 		  (nxd matrix)Labels: labels of the graph
			% 		  (boolean)BinHandle: binary handler for histogram
			% @output:(figure)fig: figure object which is a histogram
			fig = figure;
			if nargin < 4
				Xvalues = [];
				interval=(max(Data) - min(Data))/size(Data,2);
				hist(Data, min(Data):interval:max(Data));
			else
				hist(Data, BinHandle);
			end
			if nargin < 3, Labels = {}; end
			if size(Labels, 2) < 1,
				Title = 'Data frequency';
			else
				Title = Labels{1};
			end
			if size(Labels, 2) < 2,
				XLabel = 'Value';
			else
				XLabel = Labels{2};
			end
			if size(Labels, 2) < 3,
				YLabel = 'Frequency';
			else
				YLabel = Labels{3};
			end
			title(Title, FG.FontSize, FG.FontSizeVal);
			xlabel(XLabel, FG.FontSize, FG.FontSizeVal);
			ylabel(YLabel, FG.FontSize, FG.FontSizeVal);
			set(gca, FG.FontSize, FG.FontSizeVal);
		end

		function fig = plotRMSEMoleFracAndSpectra(FG, Robs, Rcalc, Xobs, Xcalc)
			% @formula: FigureGenerator.plotRMSEMoleFracAndSpectra(Robs, Rcalc, Xobs, Xcalc)
			% @example:
			%       FG = FigureGenerator;
			%       fig = FG.plotRMSEMoleFracAndSpectra(Robs, Rcalc, Xobs, Xcalc);
			% 		(Robs n Rcalc are 1x3 matrices, while Xobs n Xcalc are 1153x3 matrices, for example.)
			% @desc:  plot data to show original(correct answer) data and calculated data
			% @input: (cxn matrix)Robs: R observed
			%		  (cxn matrix)Rcalc: R calculated
			%		  (dxn matrix)Xobs: Xmix observed
			%		  (dxn matrix)Robs: Xmix calculated
			% @output:(Figure) return figure object to save outside of this function
			eval = common.Evaluation;
			RMSEMoleFrac = eval.computeRMSE(Robs, Rcalc);% supposed to be cx1 matrix
			RMSESpectra  = eval.computeRMSE(Xobs, Xcalc);% supposed to be cx1 matrix
			StructM = FG.generateInputStructure(RMSEMoleFrac);
			StructS = FG.generateInputStructure(RMSESpectra);
			Labels = {'RMSE plot', 'RMSE of MoleFraction', 'RMSE of Spectra'};
			Color = {'ro','bx','m*'};
			fig = FG.plotWithDiagnalLine(StructM, StructS, Labels, {'RMSE R-X plot'}, Color);
		end

		function fig = plotSeveralBandedSpectra(FG, xaxis, Region, Spectra, Labels, NumSpectra, Legends, CustomLineStyle, option)
			% @formula: FigureGenerator.plotSeveralBandedSpectra(xaxis, Region, Spectra, Labels, NumSpectra, Legends, CustomLineStyle)
			% @example:
			%       ra = [563:670 874:906];
			%       Region = struct('x1', ra, 'x2', ra, 'x3', ra);
			%       Spectra= struct('x1', RXpred_lin, 'x2', RXpred_non, 'x3', Xmix);
			%       FG.plotSeveralBandedSpectra(xaxis, Region,  Spectra,...
			%           {'Lambert-Beers law Spectra (Linear IOT vs. Nonlinear IOT)', 'cm-1', '(Xmix - f(Rpred*Xpure)).^2'}, 1, {'Linear', 'NonLinear', 'Original'});
			% @desc:  plot spectra data easily
			% @input: (1xd matrix)xaxis: xaxis numbers of the plot, which mean wave numbers
			%         (Struct)Region: structure containing regional data
			%           (1xd integer matrix)Region.x*: incremental index for plotting figures
			%         (Struct)Spectra: spectral data to plot
			%           (nxd matrix)Spectra.x*: x* has the same length with xaxis.x*
			%                                   sometimes Spectra.x* are expected to have n components' data.
			%         (StrCell)Labels: {title, xlabel, ylabel}. Is added in this order.
			%         (scholar)NumSpectra: how many spectra will be shown
			%         (StrCell)Legends: the permutation of legend names
			%         (StrCell, optional)CustomLineStyle: the default line style is the FG.DashedLine.
			%                                             To replace this, insert StrCell of '-m' etc. into CustomLineStyle
            %         option.linewidth: line width
			% @output:(Figure) return figure object to save outside of this function
			fig = figure;
			util = common.Utility;
			cnt = 1;
			Plot = struct();
            if nargin < 9, option = struct(); end
			if nargin < 8,
				LineStyle = FG.DashedLine;
				CustomLineStyle = {};
			else
				LineStyle = CustomLineStyle;
            end
            if iscell(Region)
                Regiontmp = struct();
                for i = 1:length(Region)
                    Regiontmp.(['x' num2str(i)]) = Region{i};
                end
                Region = Regiontmp;
            end
            if iscell(Spectra)
                Spectratmp = struct();
                for i = 1:length(Spectra)
                    Spectratmp.(['x' num2str(i)]) = Spectra{i};
                end
                Spectra = Spectratmp;
            end
			for pair = fieldnames(Region)',
				SpectraUtil = util.splitVector(Region.(pair{1}));
				for idx = fieldnames(SpectraUtil.index)',
					% plot one series in one color
					if size(Spectra.(pair{1}), 1) > 1,
						if size(Spectra.(pair{1}), 1) < size(NumSpectra, 2),
							IndexSpectrum = 1:size(Spectra.(pair{1}), 1);
						else
							IndexSpectrum = NumSpectra;% sometimes here is 1:8, and so on.
						end
					else
					    IndexSpectrum = 1;
					end
					if size(Region.(pair{1}), 2) ~= size(Spectra.(pair{1}), 2),
					    Absorption = Spectra.(pair{1})(IndexSpectrum, SpectraUtil.value.(idx{1}));
					else
					    Absorption = Spectra.(pair{1})(IndexSpectrum, SpectraUtil.index.(idx{1}));
					end
					% for normal spectra
                    Plot.(pair{1}) = plot(xaxis(SpectraUtil.value.(idx{1})),...
                        Absorption, LineStyle{cnt});
                    if isfield(option, 'linewidth')
                        if isvector(option.linewidth)
                           words = strsplit(pair{1}, 'x');
                           set(Plot.(pair{1}), FG.LineWidth, option.linewidth(str2num(words{2})));
                        else
                           set(Plot.(pair{1}), FG.LineWidth, option.linewidth);
                        end
                    else
                        set(Plot.(pair{1}), FG.LineWidth, FG.LineWidthVal);
                    end

					% for Lambda
					% Plot.(pair{1}) = plot(xaxis(SpectraUtil.index.(idx{1})),...
					%     Absorption, LineStyle{cnt});
					hold on;
				end
				cnt = cnt + 1;
			end
			Series = zeros(1, cnt-1);
			for I = 1:cnt-1,
				if size(Plot.([FG.IdxName num2str(I)]), 1) > 1,
					% if the plot has several curves, indicate one scholar to get a legend on it.
					Series(I) = Plot.([FG.IdxName num2str(I)])(1, 1);
				else
					Series(I) = Plot.([FG.IdxName num2str(I)]);
				end
			end
			title(Labels{1}, FG.FontSize, FG.FontSizeVal);
			if size(Labels, 2) < 2,
				XLabel = 'wave number(cm^-^1)';
			else
				XLabel = Labels{2};
			end
			if size(Labels, 2) < 3,
				YLabel = 'absorption';
			else
				YLabel = Labels{3};
			end
			xlabel(XLabel, FG.FontSize, FG.FontSizeVal);
			ylabel(YLabel, FG.FontSize, FG.FontSizeVal);
			if length(Legends)>0, legend(Series, Legends, FG.FontSize, FG.FontSizeVal); end
			set(gca, FG.FontSize, FG.FontSizeVal);
		end

		function fig = plotWithDiagnalLine(FG, X, Y, Labels, Legends, CustomPlotStyle, hasDiagnalLine, axises)
			% @formula: FigureGenerator.plotWithDiagnalLine(X, Y, Labels, Legends, CustomPlotStyle, hasDiagnalLine, axises)
			% @example:
			%       RMSE_R = struct('x1', RMSE_R1, 'x2', RMSE_R2, 'x3', RMSE_R3);
			%       RMSE_X = struct('x1', RMSE_X1, 'x2', RMSE_X2, 'x3', RMSE_X3);
			%       FG = FigureGenerator;
			%       FG.plotWithDiagnalLine(RMSE_R, RMSE_X, {'RMSEplot', 'RMSE-R', 'RMSE-Xmix'}, {'AllRegion', 'NewRegion', 'ReferenceRegion'});
			% @desc:  plot data to show original(correct answer) data and calculated data
			% @input: (Struct): structure containing regional data
			%           (1xd matrix).x*: incremental index for plotting figures
			%         (Struct)Y: spectral data to plot
			%           (nxd matrix)Y.x*: x* has the same length with xaxis.x*
			%                                   sometimes Spectra.x* are expected to have n components' data.
			%         (StrCell)Labels: {title, xlabel, ylabel}. Is added in this order.
			%         (StrCell)Legends: the permutation of legend names
			% @output:(Figure) return figure object to save outside of this function
			fig = figure;
			cnt = 1;
			Plot = struct();
            if nargin < 8, axises = [0 1]; end
			if nargin < 7,
				hasDiagnalLine = true;
			end
			if nargin < 6,
				PlotStyle = FG.PlotStyle;
				CustomPlotStyle = {};
			else
				PlotStyle = CustomPlotStyle;
			end
            if iscell(X)
                Xtmp = struct();
                for i = 1:length(X)
                    Xtmp.(['x' num2str(i)]) = X{i};
                end
                X = Xtmp;
            end
            if iscell(Y)
                Ytmp = struct();
                for i = 1:length(Y)
                    Ytmp.(['x' num2str(i)]) = Y{i};
                end
                Y = Ytmp;
            end
			if hasDiagnalLine,
                % draw the diagnal line
%                 for i = 1:length(fieldnames(X))
%                     minimum = mean(mean(X.(['x' num2str(i)])));
%                     maximum = mean(mean(X.(['x' num2str(i)])));
%                     maximum = max(maximum, max(max(X.(['x' num2str(i)]))) );
%                     maximum = max(maximum, max(max(Y.(['x' num2str(i)]))) );
%                     minimum = min(minimum, min(min(X.(['x' num2str(i)]))) );
%                     minimum = min(minimum, min(min(Y.(['x' num2str(i)]))) );
%                 end
% 				axises = [minimum maximum];
				plot(axises, axises, '-k');
                axis([axises axises])
				hold on;
            end

            for elem = fieldnames(Y)',
				Plot.(elem{1}) = plot(X.(elem{1}), Y.(elem{1}), PlotStyle{cnt});
				set(Plot.(elem{1}), FG.MarkerLineWidth, FG.MarkerLineWidthVal);
				set(Plot.(elem{1}), FG.MarkerSize, FG.MarkerSizeVal);
				hold on;
				cnt = cnt + 1;
			end
			Series = zeros(1, cnt-1);
			for I = 1:cnt-1,
				if size(Plot.([FG.IdxName num2str(I)]), 1) > 1,% if the plot has several curves,
					Series(I) = Plot.([FG.IdxName num2str(I)])(1, 1);
				else
					Series(I) = Plot.([FG.IdxName num2str(I)]);
				end
			end
			if size(Labels, 2) < 1,
				Title = 'RR plot';
			else
				Title = Labels{1};
			end
			if size(Labels, 2) < 2,
				XLabel = 'RMSE(R)';
			else
				XLabel = Labels{2};
			end
			if size(Labels, 2) < 3,
				YLabel = 'RMSE(Xmix)';
			else
				YLabel = Labels{3};
			end
			title(Title, FG.FontSize, FG.FontSizeVal);
			xlabel(XLabel, FG.FontSize, FG.FontSizeVal);
			ylabel(YLabel, FG.FontSize, FG.FontSizeVal);
			if length(Legends)~=0, legend(Series, Legends, FG.FontSize, FG.FontSizeVal); end
			set(gca, FG.FontSize, FG.FontSizeVal);
			axis square;
		end

		function rrplot(Struct, Rval, NumSpectra, LegendName)
			% @formula: FigureGenerator.yyplotSameStructure(Struct, NumSpectra, Legends)
			% @desc:
			% @input: (Struct)Struct: data set which have same size.
			%         (scholar)NumSpectra: indicator showing where you want to pick up
			%         (StrCell)Legends:
			% @output:(void) show yy plot.
			figure;
			Color = {'xm', '*r', '+b'};
			cnt = 1;% counter
			Graph = struct();
			Legends = zeros(1, numel(Struct));
			for elem = fieldnames(Struct)',
			Graph.(['h' cnt]) = plot(Struct.(elem{1})(NumSpectra, 1), Rval(NumSpectra, 1), Color{cnt});
			hold on;
			Legends(cnt) = Graph.(['h' cnt]);
			cnt = cnt + 1;
			end
			% draw a diagnal line.
			axises = 0:0.0001:.5;
			plot(axises, axises, '-k');
			% set graph labels.
			PropName = 'FontSize';PropVal = 12;
			xlabel('Rpred', PropName, PropVal);
			ylabel('R', PropName, PropVal);
			title('Mole Fraction Accuracy', PropName, PropVal);
			legend(Legends,LegendName);
			set(gca, FG.FontSize, FG.FontSizeVal);
			axis square;
		end
		function Object = generateInputStructure(FG, Data1, Data2, Data3, Data4, Data5, Data6)
			% @formula: Object = generateInputStructure(FG, Data1, Data2, Data3, Data4, Data5, Data6)
			% @example:
			% 		FG = FigureGenerator();% construct an object
			%       RMSE_R = FG.generateInputStructure(RMSE_R1, RMSE_R2, RMSE_R3);
			%       RMSE_X = FG.generateInputStructure(RMSE_X1, RMSE_X2, RMSE_X3);
			%       FG.plotWithDiagnalLine(RMSE_R, RMSE_X, {'RMSEplot', 'RMSE-R', 'RMSE-Xmix'}, {'AllRegion', 'NewRegion', 'ReferenceRegion'});
			%
			% @desc:  generate input structure for FigureGenerator's other methods.
			% @input: (Undefined)Data*: original data which is supposed to convert into structure.
			% @output:(Struct)Object: Input is held in Object.x1, Object.x2 and so on.
			MaxInputIdx = 7;
			if nargin < MaxInputIdx,
				for I=nargin:MaxInputIdx-1,
					eval(['Data' num2str(I) ' = [];']);
				end
			end
			Object = struct();
			for I = 1:nargin-1,
				eval(['Object.' FG.IdxName num2str(I) ' = Data' num2str(I) ';']);
			end
		end

	end% end methods

end% end classdef
