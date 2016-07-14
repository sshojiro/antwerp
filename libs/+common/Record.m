classdef Record
    properties(GetAccess=public, SetAccess=private)
        t;
    end
    methods
        function self = Record(X, y)
            self.t = table();
            self.t.X = X;
            self.t.y = y;
        end
    end
end