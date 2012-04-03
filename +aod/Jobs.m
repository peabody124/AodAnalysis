%{
aod.Jobs (computed) # A scan site
->aod.Cell
<<JobFields>>
%}

classdef Jobs < dj.Relvar
    properties(Constant)
        table = dj.Table('aod.Jobs');
    end
    
    methods 
        function self = Jobs(varargin)
            self.restrict(varargin{:})
        end
    end    
end
