%{
aod.MovieInformation (imported) # Information about which movie is played

->aod.Scan
---
stim_trace                       : longblob     # The trace indicating which presentation
movie_num                        : longblob     # The movie for each presentation
natural_movie                    : longblob     # Natural or phase scrambled
dt                               : double       # The time step of the traces (in ms)
%}

classdef MovieInformation < dj.Relvar
    properties(Constant)
        table = dj.Table('aod.MovieInformation');
    end
    
    methods 
        function self = MovieInformation(varargin)
            self.restrict(varargin{:})
        end
    
        function X = get_regressor(this,dt,bins_per_movie)
            
            assert(count(this) == 1, 'Only for one object');
            
            dat = fetch(this,'*');
            
            % Resample stimulus types
            t_orig = (0:length(dat.stim_trace)-1) * dat.dt;
            t_new = 0:dt:t_orig(end);
            stimTrace = interp1(t_orig,dat.stim_trace,t_new,'nearest');
            
            % Determine block of regressors for each movie
            [~,~,idx] = unique(dat.movie_num);
            movieNum = idx + double(dat.natural_movie) * max(idx);

            % Find when movies start and stop
            a = stimTrace > 0;
            fidx = find(diff(a)==-1);
            ridx = find(diff(a)==1);
            
            % Work out bin size to give appropriate number of regressors 
            % per movie
            bin_size = ceil(max(fidx-ridx) / bins_per_movie);
            
            X = zeros(length(stimTrace),0);
            for i = 1:length(stimTrace)
                if stimTrace(i) > 0
                    stimType = movieNum(stimTrace(i));
                    stimTime = min(floor((i - ridx(find(ridx < i, 1, 'last'))) / bin_size), bins_per_movie-1);
                    X(i,1+(stimType-1) * bins_per_movie + stimTime) = 1;
                    a(i) = stimTime;
                end
            end

        end

    end
    
    methods(Static)
        function to_insert = formatData(dat)
            % which movies are natural
            to_insert.natural_movie = cellfun(@(x) strcmp(x, 'natural')==1, dat.movieTypes);
            to_insert.movie_num = dat.movieNum;
            to_insert.stim_trace = dat.stimTrace;
            to_insert.dt = 1/dat.fps;
        end
    end
end
