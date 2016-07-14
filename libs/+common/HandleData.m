classdef HandleData
%  legacy source code which is not recommended to use.
%  EXAMPLE:
%    hdata = common.HandleData;
	methods(Static)
		%%% SMOOTHING %%%

		function [Q2q, lmnum, X, y] = smoothenWithSavitzkyGolay(winnum, eqorder,difforder, X, y, hasFigure)
			% @formula: [Q2q, lmnum, X, y] = smoothenWithSavitzkyGolay(winnum, eqorder,difforder, X, y, hasFigure)
			% @desc:  Method to compute RMSE(Root Mean Square Error).
			% @input: (nxd matrix)Y: observed value of Y.
			% 		  (nxd matrix)Ycalc: predicted value of Y.
			% 		  (boolean)isNormal: if is normal RMSE, fraction is changed.
			% 		  (boolean)hasFigure: if you want to show a figure of RMSE, input 1 in here.
			% @output:(scholar)s: calculated RMSE value.
			obsx = X(:,:);
			if nargin < 6, hasFigure = false; end
			obsx_sg = sgolay( obsx, eqorder, difforder, winnum );
			if hasFigure,
				% plot original data.
				figure;
				plot( lamda, obsx );
				figsyori_k( '波長 [nm]', '吸光度(処�?��)', 20 );
				% plot smoothened data.
				figure;
				plot( lamda, obsx_sg );
				figsyori_k( '波長 [nm]', '吸光度(処�?�?', 20 );
			end
			fold = 5;
			X = obsx_sg;
			[StdRegCoef, RegCoef, Const, Y_R2, Q2q, Ycalc, Ypred, T, P]...
				= pls_all2_k(X, y, min(30,rank(X)),fold, true, false);
			lmnum = find(Q2q==max(Q2q));% determine the optimal number of component(ONC).
			bb = RegCoef(:,lmnum);% pick regression coefficient at ONC
			cc = Const(:,lmnum);% pick constant at ONC
			ycalc = X*bb + cc;% calculate y at ONC
			yyplot_k( y, ycalc );% compare calculated and observed values.
		end

		function [Spectra] = calculateBoxCoxedSpectra(Xpure, SpectraRange, Ratio, Lambda, temp_times)
			% @formula: [Spectra] = calculateBoxCoxedSpectra(Xpure, SpectraRange, Ratio, Lambda, temp_times)
			% @desc:  Method to compute spectra based on Box-Cox Transformation.
			% @input: (2xd matrix)Xpure: spectra of pure components.
			% 							 2 components are available in binary mixture system.
			% 		  (1xd matrix)SpectraRange: selected region(cm-1)
			% 		  (nx2 matrix)Ratio: mole fractions of pure components.
			% 							 2 components are available in binary mixture system.
			% 		  (nxd matrix)Lambda: BoxCox parameter
			% 		  (scholar)temp_times: scaling factor to increase the sensitivity of optimization
			% @output:(nxd matrix)Spectra: calculated Spectra
			if nargin < 5, temp_times = 0; end
			Xpure_sc = Xpure + temp_times;
			RXpred_sc = Ratio * Xpure_sc(:,SpectraRange);
			[M,N]=size(RXpred_sc);
			RXpred = zeros(M, N);
			for I=1:N,
				for J=1:M,
					if size(Lambda, 1)~=1,
						ElemLambda = Lambda(J, I);
					else
						ElemLambda = Lambda(1, I);
					end
					RXpred(J, I) = ((RXpred_sc(J,I)).^ElemLambda-1.0)./ElemLambda;
				end
			end;
			Spectra = RXpred - temp_times;
		end
	end% end methods
end% end classdef