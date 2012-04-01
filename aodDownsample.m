function [traces t] = aodDownsample(traces,t,ds,dfof)

if nargin < 4; dfof = true; end

f = factor(ds);

for i = 1:size(traces,2)
    temp = traces(:,i);
    for j = 1:length(f)
        temp = decimate(temp,f(j));
    end
    traces2(:,i) = temp;
end

t = t(1:ds:end);
% temp = t;
% for j = 1:length(f)
%     temp = decimate(temp,f(j));
% end
% 
% t = temp;
traces = traces2;

if dfof
    m = mean(traces,1);
    traces = bsxfun(@rdivide,bsxfun(@minus,traces,m),m);
end