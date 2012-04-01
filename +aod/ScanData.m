%{
aod.ScanData (computed) # Imports the data from a scan

->aod.Scan
---
num_cells        : int unsigned          # The mouse number of cells
%}

classdef ScanData < dj.Relvar & dj.AutoPopulate
    properties(Constant)
        table = dj.Table('aod.ScanData');
        popRel = aod.Scan;
    end
    
    methods 
        function self = ScanData(varargin)
            self.restrict(varargin{:})
        end
                
        function makeTuples( this, key )
            % Not written yet
            tuple = key;
            
            [traces t coordinates] = getTraces(aod.Scan & key);
            tuple.num_cells = size(traces,2);
            insert(this,tuple);

            dt = mean(diff(t));
            for i = 1:size(traces,2)
                cell = key;
                cell.cell_num = i;
                cell.x = coordinates(i,1);
                cell.y = coordinates(i,2);
                cell.z = coordinates(i,3);
                cell.dt = dt;
                cell.trace = traces(:,i);
                insert(aod.Cell, cell);
            end
            	    
        end
    end
end
