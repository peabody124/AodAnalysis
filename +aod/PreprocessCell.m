%{
aod.PreprocessCell (computed) # A scan site

->aod.PreprocessScan
->aod.Cell
---
trace                       : longblob                      # The preprocessed trace
dt                          : double                        # The time step of the traces (in ms)
%}

classdef PreprocessCell < dj.Relvar
    properties(Constant)
        table = dj.Table('aod.PreprocessCell');
    end
    
    methods 
        function self = PreprocessCell(varargin)
            self.restrict(varargin{:})
        end

        function makeTuples(self, key)
            tuple = key;
            
            % Get the data
            method_name = fetch1(aod.PreprocessMethod & key,'preprocess_method_name');
            [raw_trace dt] = fetch1(aod.Cell & key, 'trace', 'dt');
            
            if strcmp(method_name,'raw') == 1
                tuple.trace = raw_trace;
                tuple.dt = dt;
            elseif strcmp(method_name,'ds20') == 1 % Downsample to 20 Hz
                ds = round((1/20)/dt);
                assert(ds > 1);
                
                tuple.trace = decimate(raw_trace,ds,'fir');
                tuple.trace = (tuple.trace - mean(tuple.trace)) / mean(tuple.trace);
                tuple.dt = dt * ds;
            elseif strcmp(method_name,'ds5') == 1 % Downsample to 5 Hz
                ds = round((1/5)/dt);
                assert(ds > 1);
                
                tuple.trace = decimate(raw_trace,ds,'fir');
                tuple.trace = (tuple.trace - mean(tuple.trace)) / mean(tuple.trace);
                tuple.dt = dt * ds;
            elseif strcmp(method_name,'ds20_mc') == 1 % Downsample to 20 Hz after removing mean
                ds = round((1/20)/dt);
                assert(ds > 1);
                
                raw_traces = fetchn(aod.Cell & aod.ScanData(key), 'trace');
                mean_trace = mean(cat(2,raw_traces{:}),2);

                tuple.trace = decimate(raw_trace - mean_trace,ds);
                tuple.trace = (tuple.trace - mean(tuple.trace)) / mean(tuple.trace);
                tuple.dt = dt * ds;
            elseif strcmp(method_name,'fast_oopsi') == 1 % Run the voegelstein fast oopsi method                
                ds = round((1/20)/dt);
                assert(ds > 1);
                
                ds_trace = decimate(raw_trace,ds);
                
                highPass = 0.1;

                k = hamming(round(1/(dt*ds)/highPass)*2+1);
                k = k/sum(k);
                trace = ds_trace - convmirr(ds_trace,k);  %  dF/F where F is low pass
  
                [tuple.trace P_best] = fast_oopsi(trace, struct('dt',dt * ds));
                tuple.dt = dt * ds;
            elseif strcmp(method_name,'bandpass') == 1 % Downsample to 20 Hz and apply a HPF
                ds = round((1/20)/dt);
                assert(ds > 1);

                highPass = 0.1;
                
                raw_trace = decimate(raw_trace,ds);
                
                k = hamming(round(1/(dt*ds)/highPass)*2+1);
                k = k/sum(k);
                tuple.trace = raw_trace - convmirr(raw_trace,k);
                tuple.dt = dt * ds;
            elseif strcmp(method_name,'extracted_events') == 1 % Take events where mean activity above threshold
                
                fastoopsi_key = setfield(key, 'preprocess_method_num', ...
                             fetch1(aod.PreprocessMethod('preprocess_method_name="fast_oopsi"'), ...
                             'preprocess_method_num'));
                
                assert(count(aod.PreprocessScan(fastoopsi_key))==1);
                
                % Preprocessed mean, use this to find active periods
                pm = fetch1(aod.PreprocessMean(fastoopsi_key),'mean_trace');
                thresh = quantile(pm, 0.9);
                thresholded = pm > thresh;
                
                trace = fetch1(aod.PreprocessCell(fastoopsi_key),'trace');
                
                ridx = find(diff(thresholded) == 1);
                fidx = find(diff(thresholded) == -1);
                
                if fidx(1) < ridx(1), fidx(1) = []; end
                if length(ridx) > length(fidx), ridx(end) = []; end
                
                for i = 1:length(ridx)
                    tuple.trace(i,1) = mean(trace(ridx(i):fidx(i)));
                end
                tuple.dt = 0;
            elseif strcmp(method_name,'above_threshold') == 1 % Take events where mean activity above threshold
                
                fastoopsi_key = setfield(key, 'preprocess_method_num', ...
                             fetch1(aod.PreprocessMethod('preprocess_method_name="fast_oopsi"'), ...
                             'preprocess_method_num'));
                
                assert(count(aod.PreprocessScan(fastoopsi_key))==1);
                
                % Preprocessed mean, use this to find active periods
                pm = fetch1(aod.PreprocessMean(fastoopsi_key),'mean_trace');
                thresh = quantile(pm, 0.8);
                
                [trace dt] = fetch1(aod.PreprocessCell(fastoopsi_key),'trace','dt');
                                
                tuple.trace(:,1) = trace(pm > thresh);
                
                tuple.dt = dt;
            else
                assert(false);
            end
            
            insert(aod.PreprocessCell,tuple);
        end
    end
end
