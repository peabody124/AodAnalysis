%{
aod.Mouse (manual) # A mouse

mouse_id        : int unsigned          # The mouse id
datetime        : datetime          # readable format of session start
---
%}

classdef Mouse < dj.Relvar
    properties(Constant)
        table = dj.Table('aod.Mouse');
    end
    
    methods 
        function self = Mouse(varargin)
            self.restrict(varargin{:})
        end 
    end
end
