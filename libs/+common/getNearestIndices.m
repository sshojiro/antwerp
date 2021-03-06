function interval = getNearestIndices(c, target)
%% interval = getNearestIndices(c, target)
% 
% @desc:
% ある数字が含まれるvectorのmin, maxを引き取る
% 
% @input:
% c = 991;% 要素番号上での数字を入れる
% target = [1:15 985:999 1100:1111];% 要素番号へ変換してから使う
% 
% @output:
% interval: lower/upper boundsを得る % get 985, 999

%% implementation
for d = [1 2]
    cnt = 1;
    % d==1 => upward
    % d==2 => downward
    switch(d)
        case 1
            region = c:max(target);
        case 2
            region = c:-1:min(target);
    end
    
    for idx = region
        iscontinuous = length(unique(target==idx)) == 2;
        if idx == max(target) || idx == min(target)
            iscontinuous = false;
        end
        if ~iscontinuous
            i = find(region==idx);
            switch(d)
                case 1
                    interval.upper = region(i);
                case 2
                    interval.lower = region(i);
            end
            break;
        end
    end
end

