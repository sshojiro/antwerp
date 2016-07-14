function vec = cyclic(v)
    %% cyclic
    % vec = cyclic(v)
    % common.Utility.cyclic([1:3])
    % => [1 2 3 2 1 3 3 1 2]
    if size(v, 1) == length(v),
        vin = v';
    else
        vin = v;
    end
    vec = zeros(1, length(v).^2);
    for i = 1:length(v)
        cnt = i;
        for j = 1:length(v)
            vec((i-1)*length(v)+j) = v(cnt);
            cnt = cnt + 1;
            if cnt > length(v), cnt = 1; end
        end
    end
end