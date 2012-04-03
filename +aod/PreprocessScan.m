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
    
    methods (Access=protected)
        function makeTuples( this, key )
            % Not written yet
            tuple = key;
            
            insert(this,tuple);

            cells = dj.struct.join(key, fetch(aod.Cell & key));
            for cell = cells'
                makeTuples(aod.PreprocessCell, cell);
            end
        end
    end
    
    methods
        function self = PreprocessScan(varargin)
            self.restrict(varargin{:})
        end
        
        function [traces t] = getArray( this )
            % [dat t] = getArray( relvar ) -- get the preprocessed traces
            % and return as an array
            assert(count(this) == 1);
            key = fetch(this);
            dat = fetch(aod.PreprocessCell(key), 'trace', 'dt');
            dat = dj.struct.sort(dat,{'cell_num'});
            traces = cat(2,dat.trace);
            t = (1:size(traces,1)) * dat(1).dt;
        end
        
    end
end
