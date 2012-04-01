%{
aod.Scan (manual) # A scan site

->aod.Mouse
scan_datetime             : datetime     # readable format of session start
---
aod_scan_path             : varchar(255) # path to the ephys data
%}

classdef Scan < dj.Relvar
    properties(Constant)
        table = dj.Table('aod.Scan');
    end
    
    methods 
        function self = Scan(varargin)
            self.restrict(varargin{:})
        end
        
        function fn = getFileName(self)
            % Determine the file name
            assert(count(self) == 1, 'Only for one file');
            fn = getLocalPath(fetch1(self,'aod_scan_path'));
        end
        
        function [traces t coordinates] = getTraces(self)
            % Try the various methods to load the data
            try
                [traces t coordinates] = loadTraces(getFileName(self));
            catch
                [traces t coordinates] = loadAODTraces(getFileName(self));
            end
        end
        
        function [x y z t] = computeMotion(self)
            % Compute the motion for this site
            [x y z t] = trackMotion(getFileName(self));
        end
    end
end
