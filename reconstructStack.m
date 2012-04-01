function im = reconstructStack(stack,shifts,range)

%[foo idx] = max(reshape(L1(:,:,:,1)+L2(:,:,:,1),size(L1,1) * size(L2,2), []));
%[shifts(:,1) shifts(:,2)] = ind2sub([size(L1,1) size(L1,2)],idx);

if nargin < 3; range = max(abs(shifts(:))); end

% for now not  trying to interpolate image
shifts = round(shifts);

im = zeros(size(stack,1)+2*range,size(stack,2)+2*range);
count = im;

for i = 1:size(stack,3)
    idxI = range + shifts(i,1) + (1:size(stack,1));
    idxJ = range + shifts(i,2) + (1:size(stack,2));
    im(idxI,idxJ) = im(idxI,idxJ) + stack(:,:,i);
    count(idxI,idxJ) = count(idxI,idxJ) + 1;
end

im = im ./ count;
im(count < 100) = mean(im(count(:) > 100));
im(im(:) == 0) = 1;