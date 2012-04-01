%{
aod.Cell (imported) # A mouse

->aod.ScanData
cell_num       : int unsigned          # The id of the cell
---
x              : double                # The id of the cell
y              : double                # The id of the cell
z              : double                # The id of the cell
dt             : double                # The time step of the traces (in ms)
trace          : longblob              # The calcium trace
%}

classdef Cell < dj.Relvar
    properties(Constant)
        table = dj.Table('aod.Cell');
    end
    
    methods 
        function self = Cell(varargin)
            self.restrict(varargin{:})
        end 
    end
end
