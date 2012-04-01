%{
aod.Motion (imported) # A mouse

->aod.ScanData
---
x              : longblob              # The id of the cell
y              : longblob              # The id of the cell
z              : longblob              # The id of the cell
t              : longblob              # The time step of the traces (in ms)
%}

classdef Motion < dj.Relvar & dj.AutoPopulate
    properties(Constant)
        table = dj.Table('aod.Motion');
        popRel = aod.Scan;
    end
    
    methods 
        function self = Motion(varargin)
            self.restrict(varargin{:})
        end
        
        function makeTuples(self, key)
            tuple = key;
            [tuple.x tuple.y tuple.z tuple.t] = computeMotion(aod.Scan & key)
            insert(aod.Motion, tuple);
        end
    end
end
