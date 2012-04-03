%{
aod.ConnectedGLM (imported) # Information about which movie is played

->aod.PreprocessCell
->aod.MovieInformation
---
ones                : double       # The column of ones
mean                : double       # The mean population input
history             : longblob     # The trace indicating which presentation
interactions        : longblob     # The movie for each presentation
stimulus            : longblob     # Natural or phase scrambled
fitinfo             : longblob     # Information about the fit
dt                             : double                        # The time step of the traces (in ms)
%}

classdef ConnectedGLM < dj.Relvar & dj.AutoPopulate
    properties(Constant)
        table = dj.Table('aod.ConnectedGLM');
        popRel = aod.PreprocessCell .* (aod.Scan .* aod.MovieInformation) .* aod.PreprocessMethod('preprocess_method_name="fast_oopsi"');
    end
    
    methods 
        function self = ConnectedGLM(varargin)
            self.restrict(varargin{:})
        end
    end
    
    methods (Access=protected)
        
        function self = makeTuples(self, key)
            % Compute a GLM fit for a single cell including coupling terms
            % stimulus terms and history terms

            tuple = key;
            tuple.dt = fetch1(aod.PreprocessCell(key),'dt');

            % Get the data to fit
            [traces traces_t] = getArray(aod.PreprocessScan(key));
            stim = get_regressor(aod.MovieInformation(key), ...
                tuple.dt, 10);

            r = traces(:,key.cell_num);
            traces(:,key.cell_num) = [];

            mr = mean(traces,2);

            % Generate a history vector
            h = [circshift(r,1) circshift(r,2) circshift(r,3) circshift(r,4), ...
                circshift(r,5) circshift(r,6) circshift(r,7) circshift(r,8), ...
                mean([circshift(r,9) circshift(r,10) circshift(r,11) circshift(r,12)],2), ...
                mean([circshift(r,13) circshift(r,14) circshift(r,15) circshift(r,16)],2), ...
                mean([circshift(r,17) circshift(r,18) circshift(r,19) circshift(r,20)],2)];

            % Exclude the parts where history wraps arouund
            max_history = 20;
            X = [ones(length(mr),1) mr h traces stim];
            X(1:max_history,:) = [];
            r(1:max_history) = [];

            [b fitinfo] = lassoglm(X, r, 'poisson','CV',10);
            
            tuple.ones = b(1); b(1) = [];
            tuple.mean = b(2); b(2) = [];
            tuple.history = b(1:size(h,2)); b(1:size(h,2)) = [];
            tuple.interactions = b(1:size(traces,2)); b(1:size(traces,2)) = [];
            tuple.stimulus = stim;
            tuple.fitinfo = fitinfo;
            
            insert(aod.ConnectedGLM, tuple);
        end

    end
end
