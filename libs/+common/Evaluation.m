classdef Evaluation
	%% DESCRIPTION:
	%    Module for evaluation methods
	%  EXAMPLE:
	%    eval = common.Evaluation;
	%    r = eval.calcRSquare(yobs, ypred);
	methods(Static)
		%%% EVAL FUNCTIONS %%%
        function rms = computeRMS(Y)
            %% computeRMS(Y)
            % 
            % @desc: calculate root mean square of a vector or matrix Y
            % Ex)
            % common.Evaluation.computeRMS( randn(2,100) )
            % ans =
            %     1.0782
            %     1.0086
            [samplenum, ~] = size(Y);
            rms = zeros(samplenum, 1);
            rootmeansquare = @(x) sqrt( mean( x.^2 ) );
            for s = 1:samplenum
                rms(s) = rootmeansquare( Y(s,:) );
            end
        end
            
        function density = computeProbability(X, Y, x, ythreshold)
            %% density = computeProbability(X, Y, x, ythreshold)
            %
            % @desc:
            % returns 0 <= P(x|y<=ythreshold) <= 1
            %
            % @input:
            % X, Y: x and y measurements of 2D data
            % x: value
            % ythreshold: a constraint of y
            if min(X) > x
                error('x is less than the minimum of X');
            end
            density = sum(Y <= ythreshold & X <= x) ./ sum(X <= x);
        end
        
        
		function rae = computeRelAbsErr( Yobst , Ycalct )
            % :::WARNING:::
            % both parameters are not treated as not equivalent.
			% @formula: s = computeRelAbsErr( Y , Ycalc , isNormal , showFigure )
			% @desc:  Method to compute RMSE(Root Mean Square Error).
			% @input: (nxd matrix)Yobs: observed value of Y.
			% 		  (nxd matrix)Ycalc: predicted value of Y.
			% @output:(scholar)rae: returns relative absolute errors
            %
            % @example:
            % common.Evaluation.computeRelAbsErr( ...
            %   .2*ones(2,5) ,...
            %   [.1 .1 .2 .2 .4; .2 .05 .05 .1 .1])
            % 
            % ans =
            % 
            %     2.0000
            %     2.5000
            Yobs = Yobst';
            Ycalc = Ycalct';
            samplenumber = size(Yobs, 1);
            rae = zeros(samplenumber, 1);
            for idx = 1:samplenumber
                vectmp = diag(abs(Yobs(idx,:) - Ycalc(idx,:))) ...
                         ./ diag(Yobs(idx,:));
                vec = vectmp(~isnan(vectmp));
                rae(idx,:) = sum(vec);
            end
		end
        
		function mae = computeMAE( Y , Ycalc , isNormal , showFigure )
			% @formula: s = computeRMSE( Y , Ycalc , isNormal , showFigure )
			% @desc:  Method to compute RMSE(Root Mean Square Error).
			% @input: (nxd matrix)Y: observed value of Y.
			% 		  (nxd matrix)Ycalc: predicted value of Y.
			% 		  (boolean)isNormal: if is normal RMSE, fraction is changed.
			% 		  (boolean)showFigure: if you want to show a figure of RMSE, input 1 in here.
			% @output:(scholar)s: calculated RMSE value.
			if nargin < 4, showFigure = false; end
			if nargin < 3, isNormal = true; end
			[ m, n ] = size( Ycalc );
            N = size( Y , 2 );
			mae = mean(abs(Y - Ycalc))';
			if showFigure
			    figure;
			    plot( mae , 'b' );
			    xlabel('num of components');ylabel( 'RMSE' );
			    set(gcf, 'Color' , 'w' ); set(gca, 'FontSize' ,20);
			end
		end
		function rmse = computeRMSE( Y , Ycalc , isNormal , showFigure )
			% @formula: s = computeRMSE( Y , Ycalc , isNormal , showFigure )
			% @desc:  Method to compute RMSE(Root Mean Square Error).
			% @input: (nxd matrix)Y: observed value of Y.
			% 		  (nxd matrix)Ycalc: predicted value of Y.
			% 		  (boolean)isNormal: if is normal RMSE, fraction is changed.
			% 		  (boolean)showFigure: if you want to show a figure of RMSE, input 1 in here.
			% @output:(scholar)s: calculated RMSE value.
			if nargin < 4, showFigure = false; end
			if nargin < 3, isNormal = true; end
			[ m, n ] = size( Ycalc );
            N = size( Y , 2 );
			%SSE
			if N == 1, Yobs = repmat( Y , 1 , n ); else Yobs = Y; end
			SSE = diag( ( Yobs - Ycalc )' * ( Yobs - Ycalc ) ); %'
			if isNormal, m = m + 1; end
			%MSE
			rmse = sqrt( SSE / ( m - 1 ) );
			if showFigure
			    figure;
			    plot( rmse , 'b' );
			    xlabel('num of components');ylabel( 'RMSE' );
			    set(gcf, 'Color' , 'w' ); set(gca, 'FontSize' ,20);
			end
		end

	end% end methods
end% end classdef
