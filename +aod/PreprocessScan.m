%{
aod.PreprocessScan (computed) # A scan site

->aod.Scan
->aod.PreprocessMethod
---
%}

classdef PreprocessScan < dj.Relvar & dj.AutoPopulate
    properties(Constant)
        table = dj.Table('aod.PreprocessScan');
        popRel = aod.Scan * aod.PreprocessMethod;
    end
    
    methods 
        function self = PreprocessScan(varargin)
            self.restrict(varargin{:})
        end
        
        function makeTuples( this, key )
            % Not written yet
            tuple = key;
            
            insert(this,tuple);

            cells = dj.struct.join(key, fetch(aod.Cell & key));
            for cell = cells'
                makeTuples(aod.PreprocessCell, cell);
            end
        end
        
        function [dat t] = getArray( this )
            % [dat t] = getArray( relvar ) -- get the preprocessed traces
            % and return as an array
            assert(count(this) == 1);
            key = fetch(this);
            [traces dt] = fetchn(aod.PreprocessCell(key), 'trace', 'dt');
            dat = cat(2,traces{:});
            t = (1:size(dat,1)) * dt(1);
        end
        
    end
end
