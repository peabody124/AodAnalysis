function coords = loadPoints(fn)

points = loadHWS(fn,'config','points');
points = double(reshape(points,3,[])');
coords = bsxfun(@rdivide,points,[1460000 1460000 600]);

