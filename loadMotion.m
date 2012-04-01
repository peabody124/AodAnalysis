function [mot t coordinates] = loadMotion(fn);

mode = loadHWS(fn,'config','mode');
assert(~isempty(mode), 'Cannot find mode information, invalid AOD scan file');
assert(mode(1) == 0, 'Motion not found in volume scans');

coordinates = loadPoints(fn);
coordinates = coordinates(end-mode(4)+1:end,:);

[mot1 dt t0] = loadHWS(fn,'MotionData','ImCh1');
[mot2 dt t0] = loadHWS(fn,'MotionData','ImCh2');

np = double(mode(4));
columns = double(mode(5));
rows = np/columns;
frames = floor(numel(mot1)/columns/rows);

mot(1,:,:,:) = reshape(mot1(1:columns*rows*frames),[1 columns rows frames]);
mot(2,:,:,:) = reshape(mot2(1:columns*rows*frames),[1 columns rows frames]);

t = (1:size(mot,4)) * dt * np;
