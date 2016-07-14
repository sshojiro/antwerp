classdef Utility
	%% DESCRIPTION:
	%   Utility common module
	% USAGE:
	% 	import common.*;
	% 	util = Utility;
	% 	util.foreach(params, @set, args);
	methods (Static)
        function index = input(var, customstr)
            %% util.input(var, customstr)
            %
            % @example:
            % util = common.Utility;
            % util.input({[1 2 4 5 7] [7 8 10 11 13] 4:3:16}, 'select the index of the calibration sets:');
            % `var` cell is expanded into:
            % 1)
            %      1     2     4     5     7
            % 
            % 2)
            %      7     8    10    11    13
            % 
            % 3)
            %      4     7    10    13    16
            % 
            % select the index of the calibration sets:
            %      1 % for example, press 1 and then, enter
            if nargin == 1,
                custromstr = '';
                finalprompt = 'select and type the number';
            else
                finalprompt = customstr;
            end
            if isstruct(var)
                % struct
                for name = fieldnames(var)'
                    disp([' ' name{1} ':']);
                    disp(var.(name{1}));
                end
            elseif iscell(var)
                % 1-row/column cell
                if size(var,1)==1 || size(var,2)==1
                    disp('`var` cell is expanded into:');
                    for idx = 1:length(var)
                        disp([num2str(idx) ')']);
                        disp(var{idx});
                    end
                end
            else
                error('`var` must be a struct or a cell');
            end
            index = input([finalprompt '\n     ']);
        end
        function ContinuousNumbers = getContinuousRegions(boolvec)
            %% ContinuousNumbers = getContinuousRegions(vec)
            %
            % @desc:
            % [1:10 100:110 200:220] => {1:10 100:110 200:220}
            % @input:
            % boolvec: an integer vector
            if isequal(unique(boolvec), [0 1])
                nonzero = find([0 boolvec]~=0);
            else
                error(' input must be a logical, or 0/1 vector');
            end
            findStart = @(v) v( find( diff([0 v])>1 ) ) - 1; %#ok<FNDSB>
            findEnd = @(v) [v( find( diff(nonzero)>1 ) )  v(end)] - 1; %#ok<FNDSB>
            if ~isempty(nonzero)
                vec_start = findStart(nonzero);
                vec_end = findEnd(nonzero);
                GroupNumber = length(vec_start);
                ContinuousNumbers = cell(1, GroupNumber);
                for groupnum = 1:GroupNumber
                    ContinuousNumbers{groupnum} = vec_start(groupnum):vec_end(groupnum);
                end
            else
                error('no region is selected');
            end
        end
        function [trimmed, trimmed_bool] = trimUnreliableNumber(boolvec, trimming_length)
            %% trimmed = trimUnreliableNumber(boolvec)
            %
            % @desc:
            % trim #{trimming_length}-length regions
            %
            % @input:
            % boolvec: 1/0 vectors
            % trimming_length: 1, 4, 9
            if ~exist('trimming_length', 'var'), trimming_length = 1; end
            nonzero = find([0 boolvec] ~= 0);% 0でないインデックス値を取得
            if ~isempty(nonzero)
                vs = nonzero(find(diff([0 nonzero]) > 1))-1;% ベクトルの連続データの開始点を取得
                ve = [nonzero(find(diff(nonzero) > 1)) nonzero(end)]-1;% ベクトルの連続データの終了点を取得
                vres = cell(1, length(vs));
                for n=1:length(vs)
                    vres{n} = boolvec(vs(n):ve(n));% ベクトルの作成
                end
                vlength = cellfun('length', vres);% 各ベクトルの長さ
                cnt = 1;
                ret = cell(1, sum(vlength > 1));
                for regionnum = 1:length(vs)
                    if vlength(regionnum) > trimming_length
                        ret{cnt} = vs(regionnum):ve(regionnum);
                        cnt = cnt + 1;
                    end
                end
                trimmed = cell2mat(ret);
                trimmed_bool = zeros(1, length(boolvec));
                trimmed_bool(trimmed) = 1;
                trimmed_bool = logical(trimmed_bool);
            else
                % raise exception
                trimmed = [];
                trimmed_bool = boolvec;
            end
        end
        function cel = splitVectorByBreaks(vec, breaks)
            %% cel = splitVectorByBreaks(vec, breaks)
            %
            % @desc:
            % a method to split a vector by a break vector
            % 
            % @input:
            % vec: d length
            % breaks: d-1 length
            vmax = length(vec);
            num_breaks = length(breaks);
            cel = cell(1, num_breaks+1);
            for i = 1:num_breaks+1
                if i == 1
                    cel{i} = 1:breaks(i);
                elseif i == num_breaks+1
                    cel{i} = breaks(i-1)+1:vmax;
                else
                    cel{i} = breaks(i-1)+1:breaks(i);
                end
            end
        end
        function cel = extractContinuousNumber(vec)
            cnt = 1;
            for i = min(vec):max(vec)
                if i == min(vec)
                    cel{cnt} = i;
                else
                    if ~isempty(find(vec==i)) && ~isempty(find(vec==i-1))
                        cel{cnt} = [cel{cnt} i];
                    elseif ~isempty(find(vec==i)) && isempty(find(vec==i-1))
                        cnt = cnt + 1;
                        cel{cnt} = i;
                    end
                end
            end
        end
        function ret = replaceSpecificLetters(str, from, to)
            %% ret = replaceSpecificLetters(str, from, to)
            %
            % @desc: a method to replace 'from' with 'to' in 'str' letters
            %
            % @usage:
            % str = 'test_text';
            % common.Utility.replaceSpecificLetters({str 'text_test'}, {'_'}, {'-'})
            % ans = 
            %       'test-text'    'text-test'
            %
            % @input:
            % str: (cell) string that is target
            % from: (cell) letters to replace
            % to: (cell) letters to replace sth with
            %
            % @output:
            % ret: (cell) replaced string
            
            %% set default
            if nargin < 2, from = {'_'}; end
            if nargin < 3, to = {'-'}; end
            %% procedure
            ret = cell(1, length(str));
            for j = 1:length(from)
                for i = 1:length(str)
                    tmp = strsplit(str{i}, from{j});
                    ret{i} = strjoin(tmp, to{j});
                end
            end
        end
        function fnames = grepFileNames(abspath_to_dir, words)
            %% fnames = grepFileNames(abspath_to_dir, words)
            %
            % @desc:
            % grep sepcific words in filenames
            %
            % @input:
            % abspath_to_dir: (string)absolute path to directory
            % words: (cell)the words which are included in filenames
            %
            % @output:
            % fnames: (cell)filenames
            FileInfo = dir(abspath_to_dir);
            cnt = 1;
            fnames = {};
            for w = 1:length(words)
                for i = 1:length(FileInfo)
                    if  strfind(FileInfo(i).name, words{w})
                        fnames{cnt} = FileInfo(i).name;
                        cnt = cnt + 1;
                    end
                end
            end
        end
        function idx = getSortIndex(vec)
            % @desc: return sorted vector's original index.
            %
            % Example)
            % r=[0.9980;0.9940];
            % idx = common.Utility.getSortIndex(r')
            % idx =
            %      2     1
            % r(idx, :)
            % ans =
            %     0.9940
            %     0.9980
            len = length(vec);
            idx = zeros(size(vec));
            for i=1:len
                idx(i) = find(sort(vec)==vec(i));
            end
        end
        function decimal = calcDecimal(C)
            %% decimal = calcDecimal(C)
            % @input:
            % C: cell candidates
            % @output:
            % decimal: decimal number that implies the whole combination numbers.
            % common.Utility.calcDecimal({0:3 0:3 0:3 0:3 0:3})
            % ans =
            %         1024
            len = length(C);
            decimal = 1;% initialize
            for I = 1:len
                clen = length(C{I});
                decimal = decimal * clen;
            end
        end
        
        function [combc, comb, ret] = calcCombination(candidates, cnt)
           %% comb = calcCombination(candidates, cnt)
            % @input:
            % C: cell candidates
            % @output:
            % decimal: decimal number that implies the whole combination numbers.
            % c = {0:3 0:3 0:3 0:3 0:3};
            % t = common.Utility.calcCombination(c, 33)
            % ans = 
            %     [0]    [0]    [2]    [0]    [1]
            % t{3}
            % ans =
            %      2
            m = common.Utility.calcDecimal(candidates);
            if cnt <= m
                len = length(candidates);
                c   = zeros(1, len);
                for i = 1:len
                    c(i) = length(candidates{i});
                end
                combc = cell(1, len);
                comb  = zeros(1, len);
                ret   = zeros(1, len);
                for i = 1:len
                    p = 1;
                    if i~=len
                        for j = i+1:len
                            p = p * c(j);
                        end
                    end
                    combc{i} = floor((cnt - mod(cnt, p))/p);
                    comb(i)  = floor((cnt - mod(cnt, p))/p);
                    ret(i)   = candidates{i}(comb(i) + 1);
                    cnt = cnt - combc{i} * p;
                end
            end
        end
        
		function exp_mat = expandMatrix(Mat)
			% @formula: exp_mat = expandMatrix(Mat)
			% @desc:  Expand a Matrix to a vector
			% @input: (nxd matrix)Mat: Matrix
			% @output:((nxd)x1 vector)exp_mat: expanded matlab
			hig = size(Mat, 1);
			wid = size(Mat, 2);
			item_num= hig * wid;
			exp_mat = zeros(item_num, 1);
			for I=1:wid
				for J=1:hig
					exp_mat((I-1)*hig+J) = Mat(J, I);
				end
			end
		end
		function [TargetDir] = setResultDir(FilePathName)
			% @formula: DirName = setResultDir(mfilename('fullpath'))
			% @desc:  set Result Directory, making that target directory
			% @input: (String)FilePathName: script file name
			% @output:(String)TargetDir: Target Directory
			% @example:
			% 	import common.*;
			% 	util = common.Utility;
			% 	DirName = util.setResultDir(mfilename('fullpath'));
			%	@detail:
			%	FilePathName: '{$ROOT_DIR}/experiments/sample_code/exp2015_03_30_validate_with_actual_binary.m'
			%	-> ./results/sample-code/2015-03-30-validate-with-actual-binary/
			%	FilePathName: 'exp2015_03_30_validate_with_actual_binary.m'
			%	-> ./results/2015-03-30-validate-with-actual-binary/
			%	FilePathName: '{$ROOT_DIR}/experiments/test_code.m'
			%	-> ./results/test-code/
			RESULT_DIR = './results/';
            if strmatch(pwd, FilePathName)
                FilePathName = strrep(FilePathName, [pwd '\'], '');
            end
			PathNames = strsplit(FilePathName, '\');
			len = length(PathNames);
			if len > 1% when using mfilename('fullpath') as FilePathName
				if strmatch('experiments', PathNames) > 0% auto detection
					Type = 'experiments';
                elseif strmatch('scripts', PathNames) > 0% auto detection
					Type = 'scripts';
				else
					Type = 'tests';
				end
				if strcmp(PathNames{len-1}, Type)% the file is located under the experiments/ directory
					FileName = PathNames{len};
					Target = RESULT_DIR;
					mkdir(Target);
				else% the file is located under the category directory
					cnt = strmatch(Type, PathNames);
					CategoryName = PathNames{cnt+1};
					strcell = strsplit(CategoryName, '_');
					CategoryName = strjoin(strcell, '-');
					FileName = PathNames{len};
					Target = [RESULT_DIR CategoryName '/'];
				end
			else% when using just mfilename or a simple filename as FilePathName
				FileName = FilePathName;
				Target = RESULT_DIR;
			end
			prefix = regexp(FileName, '^[a-z]+[0-9]+', 'match');
			if size(prefix, 2)>0
				% case: FileName == exp2015_03_31_noise_sensitivity.m
				% 			'exp' prefix is removed here.
				prefix = regexp(FileName, '^[a-z]+', 'match');
				trimmer = strsplit(FileName, prefix);
				FileName = trimmer{2};
			end
			strcell = strsplit(FileName, '_');
			FileName = strjoin(strcell, '-');
			TargetDir = [Target FileName];
			mkdir(TargetDir);
		end
		function saveVars(FolderName, FileName, VariablesCell)
			% @formula: DirName = saveVars(mfilename('fullpath'))
			% @desc:  set Result Directory, making that target directory
			% @input: (String)FolderName: Save Folder Name
			% 		  (String)FileName: mat file name you want to save variablews into
			% 		  (Cell)VariablesCell: Variable cell you want to save
			% @output:(String)TargetDir: Target Directory
			% import common.*;
			% util = common.Utility;
			% DirName = util.setResultDir(mfilename('fullpath'))
			% util.saveVars(DirName, 'file_name', {'Var1', 'Var2'});
			src = [FolderName '/' common.Utility.addPrefixTime(FileName) '.mat'];
			Variables = '';
			for idx = 1:numel(VariablesCell)
				if idx>=2,
					Variables = [Variables ''', ''' VariablesCell{idx}];
				else
					Variables = [Variables VariablesCell{idx}];
				end
			end
			cmd = ['save(''' src ''', ''' Variables ''');' ]
			eval(cmd);
		end
		function [Mat] = convertStructToMat(Struct)
			% @formula: [Mat] = convertStructToMat(struct('x1', M1, 'x2', M2, 'x3', M3))
			% @desc: To Convert structures which have exactly the same size matrix to Matrix
			% @input: (Structure)Struct: structures which have exactly the same size matrix as values
			% @output:(nxd matrix)Mat: result matrix
			cnt = 1;
			if size(Struct.x1, 1) > 1,
				RowNum = size(Struct.x1, 1);
			else
				RowNum = size(Struct.x2, 2);
			end
			Mat = zeros(size(fieldnames(Struct), 1), RowNum);
			for Elem = fieldnames(Struct)',
				if size(Struct.(Elem{1}), 1)>1,
					Mat(cnt, :) = Struct.(Elem{1})';
				else
					Mat(cnt, :) = Struct.(Elem{1});
				end
				cnt = cnt + 1;
			end
		end
		function Struct = convertMatToStruct(Mat)
			% @desc: compute the area of spectrum. Compute numerical integral.
			% @input: (nxd matrix)Mat: a matrix which is cut into n-children structure.
			% @output: (Struct)Struct: a structure
			% 		   (1xd matrix)Struct.x*: the element of the output structure
			Struct = struct();
			for I = 1:size(Mat, 1),
				Struct.(['x' I]) = Mat(I, :);
			end
		end

		function [Area] = computeSpectralArea(xaxis, Spectrum)
			% @desc: compute the area of spectrum. Compute numerical integral.
			% @input: (1xd matrix)xaxis: spectrum measure
			%		  (1xd matrix)Spectrum: object spectrum
			% @output: (scholar)Area: sum of approximate area.
			Area = 0.0;
			delta = abs(xaxis(2) - xaxis(1));
			if size(xaxis, 2) == size(Spectrum, 2),
				for I = 1:size(xaxis, 2),
					Area = Area + Spectrum(I) * delta;
				end
			else
				% raise exception
				sprintf('make sure the dimensions of input.')
			end
		end

		function [Mat] = concatenateVectors(m1, m2)
			% @desc: merge one matrix
			% @example:
			% 		[1 5], [2 3]
			% 		=> [1 2 3 5]
			% @input: (1xd matrix)v1: a vector supposed to be greater.
			%		  (1xd matrix)v2: a vector supposed to be smaller.
			% @output: (Struct)ComparisonResult: a structure.
			% 		     (1xd boolean)isGreater: boolean vector which contains the former is greater than another.
			% 		     (1xc integer vector)Index: a vector of the indexes where the former scholar is greater than the latter one.
			d1 = size(m1, 2);
			d2 = size(m2, 2);
			Mat = zeros(1, d1+d2);
			idx1 = 1;
			idx2 = 1;
			while idx1 ~= d1 | idx2 ~= d2,
				if m1(idx1) > m2(idx2),
					Mat(1, idx1+idx2-1) = m2(idx2);
					if idx2 ~= d2,
						idx2 = idx2 + 1
					end
				else
					Mat(1, idx1+idx2-1) = m1(idx1);
					if idx1 ~= d1,
						idx1 = idx1 + 1;
					end
				end
			end
		end

		function [ComparisonResult] = compareVectors(v1, v2)
			% @desc: return the elements' index where the scholar is greater than another.
			% 		 USE this to check a spectrum achieving constraints.
			% @input: (1xd matrix)v1: a vector supposed to be greater.
			%		  (1xd matrix)v2: a vector supposed to be smaller.
			% @output: (Struct)ComparisonResult: a structure.
			% 		     (1xd boolean)isGreater: boolean vector which contains the former is greater than another.
			% 		     (1xc integer vector)Index: a vector of the indexes where the former scholar is greater than the latter one.
			ComparisonResult = struct();
			ComparisonResult.isGreater = v1 > v2;
			ComparisonResult.Index = find(sum(ComparisonResult.isGreater)==size(ComparisonResult.isGreater,1));
		end

		function [ColorSpec] = generateColorSpec(Idx, Sort)
			% @desc: To generate color spec from index number
			% @input: (Integer)Idx: ColorIndex
			% 		  (String)Sort: 'Marker', 'LineMarker', 'Line' are available as 'Sort'
            %         'SolidLine' and 'DottedLine' are also avaiable now.
			% @output:(String)ColorSpec: '+b', 'om' and so on
			if nargin < 2,
				Sort = 'Marker';
				% LineMarker
				% Line
			end
			Color  = {'m', 'c', 'r', 'g', 'b', 'k'};
			Line   = {'-', '--', ':', '-.'};
			Marker = {'+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'hexagram'};
			if strcmp(Sort, 'Marker'),% 8 * 13 combinations
				ColorSpec = [Color{1+mod(Idx, size(Color, 2))} Marker{1+mod(Idx, size(Marker, 2))}];
			elseif strcmp(Sort, 'Line'),
				ColorSpec = [Color{1+mod(Idx, size(Color, 2))} Line{1+mod(Idx, size(Line, 2))}];
			elseif strcmp(Sort, 'SolidLine'),
				ColorSpec = [Color{1+mod(Idx, size(Color, 2))} '-'];
			elseif strcmp(Sort, 'DottedLine'),
				ColorSpec = [Color{1+mod(Idx, size(Color, 2))} '-.'];
			elseif strcmp(Sort, 'LineMarker'),
				ColorSpec = [Line{1+mod(Idx, size(Line, 2))} Color{1+mod(Idx, size(Color, 2))} Marker{1+mod(Idx, size(Marker, 2))}];
			end
		end

		function Str = addPrefixDate(String)
			% @desc: return a string which is added Date Prefix to the original one
			% 		 This avoid overriding when saving into .fig or .mat file.
			% @input: (String)String: the original string
			% @output:(String)Str: result string which has date prefix
			Str = [datestr(now, 'yyyy_mm_dd_') String];
		end

		function Str = addPrefixTime(String)
			% @desc: return a string which is added Time Prefix to the original one
			% 		 This avoid overriding when saving into .fig or .mat file.
			% @input: (String)String: the original string
			% @output:(String)Str: result string which has time prefix
			Str = [datestr(now, 'yyyy_mm_dd_hh_MM_') String];
		end

		function saveJpg(fig, Location, FileName, Number, inFig, permitSaving, omittimestamp)
			% @desc:  utility to make it easy to save a figure
			% @example:
			% 		common.Utility.saveJpg(fig, Location, FileName, Number, inFig, permitSaving);
			% @input: (Figure)fig: figure object to save
			% 		(String)Location: where the file to be saved
			% 		(String)FileName: the name of the file
			%		  (Number, optional)Number: the various number to get it into the file name
			%		  (boolean, optional)inFig: save in Fig file, else, save in Jpg file.
			%		  (boolean, optional)permitSaving: permit saving.
			if nargin < 6, permitSaving = true; end
            if nargin < 7, omittimestamp = false; end
            if ~omittimestamp% default
                FileName = common.Utility.addPrefixTime(FileName);
            end
			if nargin > 3 && Number~=0,
				fname = [Location '/' FileName num2str(Number)];
			else
				fname = [Location '/' FileName];
			end
			if nargin < 5,
				Ext = 'jpg';% file extension like '.jpg'
				inFig = false;
			elseif nargin >= 5 && inFig,
				Ext = 'fig';
            elseif ~inFig,
                Ext = 'jpg';
			end
			if permitSaving,
				saveas(fig, fname, Ext);
			end
        end

		function result = splitVector(v)
			% @desc:  Split a vector into the successive smaller ones.
			% @input: (1xd integer matrix)v: vector such as a spectrum region
			% @output:(Struct)result: a structure which has vectors as elements
			% 		  (1xk matrix)result.value.v*: smaller incremental value vectors are obtained.
			% 		  (1xk matrix)result.index.i*: smaller index vectors are obtained.
			% @example:
			% 		ret = Utility.splitVector([563:670 874:906])
			% 		ret.value.v1 => [563:670]
			% 		ret.value.v2 => [874:906]
			% 		ret.index.i1 => [1:108]
			% 		ret.index.i2 => [109:141]
			[n,d] = size(v);% v is expected to be a 1xd matrix.
			result=struct('index', struct(), 'value', struct());% result is a structure.
			cnt = 1;		% abbr of counter
			idx_start = 1;  % index at start point.
			for i = 1:d,
				if i ~= d & abs(v(i) - v(i+1)) > 1,
					result.value.(['v', num2str(cnt)]) = v(idx_start:i);
					result.index.(['v', num2str(cnt)]) = idx_start:i;
					idx_start = i+1;
					cnt = cnt + 1;% update counter
				elseif i == d,
					result.value.(['v', num2str(cnt)]) = v(idx_start:i);
					result.index.(['v', num2str(cnt)]) = idx_start:i;
				end
			end
		end

		function plotRealTime(params)
			% @desc:  Method via which dots are plotted continuously.
			% 		  This is used as a wrapper function.
			% @input: (Struct)params: structure
			% 		  (n length vector)params.xcomp: component on x axis
			% 		  (n length vector)params.ycomp: component on y axis
			% 		  (String)params.symbol: symbol string such as 'xr', '-b' and so on.
			% 		  (boolean)params.isBold: if this is "true", get line bold.
			% 		  (boolean)params.isSquare: if this is "true", rescale the figure in square.
			% @output:(void) new dots are put on the figure already shown.
			if isfield(params, 'isBold') & params.isBold,
				plot(params.xcomp, params.ycomp, params.symbol, 'LineWidth', 2);
			else
				plot(params.xcomp, params.ycomp, params.symbol);
			end
			drawnow;
			if isfield(params, 'isSquare') & params.isSquare, axis square; end
		end
	end% end methods
end% end classdef
