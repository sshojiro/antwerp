function interval = getNearestIndices(c, target)
%% interval = getNearestIndices(c, target)
% 
% @desc:
% ���鐔�����܂܂��vector��min, max���������
% 
% @input:
% c = 991;% �v�f�ԍ���ł̐���������
% target = [1:15 985:999 1100:1111];% �v�f�ԍ��֕ϊ����Ă���g��
% 
% @output:
% interval: lower/upper bounds�𓾂� % get 985, 999

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

