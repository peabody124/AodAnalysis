%{
aod.MdsScale (computed) # Correlation matrix for a scan site

->aod.CorrelationMatrix
---
x                       : longblob         # The first dimension
y                       : longblob         # The second dimension
z                       : longblob         # The third dimension
%}

classdef MdsScale < dj.Relvar & dj.AutoPopulate
    properties(Constant)
        table = dj.Table('aod.MdsScale');
        popRel = aod.CorrelationMatrix;
    end
    
    methods 
        function self = MdsScale(varargin)
            self.restrict(varargin{:})
        end

        function makeTuples(self, key)
            tuple = key;
            
            c = fetch1(aod.CorrelationMatrix & key, 'corr');
            
            [p] = mdscale(abs(c),3);
            tuple.x = p(:,1);
            tuple.y = p(:,2);
            tuple.z = p(:,3);
            
            insert(aod.MdsScale,tuple);
        end
    end
end
