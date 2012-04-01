%{
aod.CorrelationMatrix (computed) # Correlation matrix for a scan site

->aod.PreprocessScan
---
corr                       : longblob                      # The preprocessed trace
cov                        : longblob                      # The preprocessed trace
%}

classdef CorrelationMatrix < dj.Relvar & dj.AutoPopulate
    properties(Constant)
        table = dj.Table('aod.CorrelationMatrix');
        popRel = aod.PreprocessScan;
    end
    
    methods 
        function self = CorrelationMatrix(varargin)
            self.restrict(varargin{:})
        end

        function makeTuples(self, key)
            tuple = key;
            
            [dat t] = getArray(aod.PreprocessScan & key);
            
            tuple.corr = corrcoef(dat);
            tuple.cov = cov(dat);
            
            insert(aod.CorrelationMatrix,tuple);
        end
    end
end
