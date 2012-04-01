%{
aod.PreprocessMean (computed) # A scan site

->aod.PreprocessScan
---
mean_trace                 : longblob                      # The mean preprocessed trace
mean_trace_time            : longblob                      # The mean preprocessed trace
%}

classdef PreprocessMean < dj.Relvar
    properties(Constant)
        table = dj.Table('aod.PreprocessMean');
    end
    
    methods 
        function self = PreprocessMean(varargin)
            self.restrict(varargin{:})
        end
        
        function makeTuples( this, key )
            % Not written yet
            tuple = key;
            
            [traces tuple.mean_trace_time] = getArray(aod.PreprocessScan(key));
            tuple.mean_trace = mean(traces,2);
            
            insert(this,tuple);
        end
        
    end
end
