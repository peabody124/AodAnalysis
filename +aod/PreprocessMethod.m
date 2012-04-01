%{
aod.PreprocessMethod (lookup) # A scan site

preprocess_method_num             : int     # readable format of session start
---
preprocess_method_name            : varchar(255) # path to the ephys data
%}

classdef PreprocessMethod < dj.Relvar
    properties(Constant)
        table = dj.Table('aod.PreprocessMethod');
    end
    
    methods 
        function self = PreprocessMethod(varargin)
            self.restrict(varargin{:})
        end
    end
end
