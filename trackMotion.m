function [xpos ypos zpos t details] = trackMotion(fn)
% Track motion in AOD traces
%
% [xpos ypos zpos t details] = trackMotion(fn)
%
% JC 2010-07-12


%% load and normalize stack
[mot t coordinates] = loadMotion(fn);

mot = double(mot);

m = mean(reshape(mot(1,:,:,:),[],size(mot,4))');
n = var(reshape(mot(1,:,:,:),[],size(mot,4))');
[b1 stats] = robustfit(m(:),n(:));

m = mean(reshape(mot(2,:,:,:),[],size(mot,4))');
n = var(reshape(mot(2,:,:,:),[],size(mot,4))');
[b2 stats] = robustfit(m(:),n(:));

% convert to photon rates
mot(1,:,:,:) = ( mot(1,:,:,:) + b1(1)/b1(2) ) / b1(2);
mot(2,:,:,:) = ( mot(2,:,:,:) + b2(1)/b2(2) ) / b1(2);

mot(mot(:) < 0) = 0;

% shifts to compute for
DeltaX = linspace(-3,3,200);
DeltaY = linspace(-3,3,200);
DeltaZ = linspace(-6,6,600);

%% compute the likelihood of the stack at each position
gridSize = size(mot,2);
offset = 3;
for i = 1:size(mot,3)/gridSize
    stack1 = squeeze(mot(1,1:gridSize,1+(i-1)*gridSize:gridSize+(i-1)*gridSize,:));
    stack2 = squeeze(mot(2,1:gridSize,1+(i-1)*gridSize:gridSize+(i-1)*gridSize,:));
    
    % first pass compute offsets from patch
    P1(:,:,:,i) = trackStack(stack1,[],-offset:offset,-offset:offset);
    P2(:,:,:,i) = trackStack(stack2,[],-offset:offset,-offset:offset);
    %P1 = 0*P2;

    % compute reference frame
    [foo idx] = max(reshape(P1(:,:,:,i)+P2(:,:,:,i),size(P1,1) * size(P2,2), []));
    [shifts(:,1) shifts(:,2)] = ind2sub([size(P1,1) size(P1,2)],idx);
    ref1 = reconstructStack(stack1,shifts-offset-1,offset);
    ref2 = reconstructStack(stack2,shifts-offset-1,offset);

    %plot(t,shifts); drawnow
    
    % interpolate ref for computing finer shifts
    up = 1;
    ref1 = interp2(ref1,up,'*cubic');
    ref2 = interp2(ref2,up,'*cubic');
    
    % compute motion against bigger reference
    up = 2^up;
    L1(:,:,:,i) = trackStack(stack1,ref1,-offset*up:offset*up,-offset*up:offset*up,up);
    L2(:,:,:,i) = trackStack(stack2,ref2,-offset*up:offset*up,-offset*up:offset*up,up);
    %L1 = 0*L2;
    
    % estimate how well model did
    m1(i,:) = max(reshape(L1(:,:,:,i),size(L1,1)*size(L1,2),[]));
    m2(i,:) = max(reshape(L2(:,:,:,i),size(L1,1)*size(L1,2),[]));
    
    % compute marginal likelihoods across shifts, interpolated up
    [likelihood_i(:,:,i) likelihood_j(:,:,i)] = likelihoodToPos(L1(:,:,:,i)+L2(:,:,:,i));

    DeltaI = linspace(-offset,offset,size(likelihood_i,2));
    DeltaJ = linspace(-offset,offset,size(likelihood_j,2));

    % project likelihoods into real coordinates
    dz_dj = mean(diff(coordinates(gridSize^2 * (i-1) + (1:gridSize:gridSize^2),3)));
    if dz_dj == 0 % horizontal plane
        dx_di = mean(diff(coordinates(gridSize^2 * (i-1) + (1:gridSize),1)));
        dy_dj = mean(diff(coordinates(gridSize^2 * (i-1) + (1:gridSize:gridSize^2),2)));
        
        likelihood_x(:,:,i) = interp1(DeltaI * dx_di,likelihood_i(:,:,i)',DeltaX,'*cubic',-100000)';
        likelihood_y(:,:,i) = interp1(DeltaJ * dy_dj,likelihood_j(:,:,i)',DeltaY,'*cubic',-100000)';        
        likelihood_z(:,:,i) = zeros([size(likelihood_i,1) length(DeltaZ)]);
    else          % vertical plane
        dx_di = mean(diff(coordinates(gridSize^2 * (i-1) + (1:gridSize),1)));
        dy_di = mean(diff(coordinates(gridSize^2 * (i-1) + (1:gridSize),2)));
        
        if max(DeltaI * dx_di) < 4
            likelihood_x(1:size(likelihood_i,1),1:length(DeltaX),i) = 0; 
        else
            likelihood_x(:,:,i) = interp1(DeltaI * dx_di,likelihood_i(:,:,i)',DeltaX,'*cubic',-100000)';
        end
        
        if max(DeltaI * dy_di) < 4
            likelihood_y(1:size(likelihood_i,1),1:length(DeltaY),i) = 0;
        else
            likelihood_y(:,:,i) = interp1(DeltaI * dy_di,likelihood_i(:,:,i)',DeltaY,'*cubic',-100000)';        
        end
        likelihood_z(:,:,i) = interp1(DeltaJ * dz_dj,likelihood_j(:,:,i)',DeltaZ,'*cubic',-100000)';        
        
    end
        
end

% for some reason the interpolation sometimes creates NaNs
likelihood_x(isnan(likelihood_x(:))) = -inf;
likelihood_y(isnan(likelihood_y(:))) = -inf;
likelihood_z(isnan(likelihood_z(:))) = -inf;

% extract maximum likelihood estimates for each plane
[foo xpos] = max(likelihood_x,[],2); xpos = squeeze(xpos);
[foo ypos] = max(likelihood_y,[],2); ypos = squeeze(ypos);
[foo zpos] = max(likelihood_z,[],2); zpos = squeeze(zpos);

% combine all planes for a single position inference
lX = sum(likelihood_x,3); [foo xpos(:,end+1)] = max(lX,[],2);
lY = sum(likelihood_y,3); [foo ypos(:,end+1)] = max(lY,[],2);
lZ = sum(likelihood_z,3); [foo zpos(:,end+1)] = max(lZ,[],2);

% convert to microns
xpos = DeltaX(xpos);
ypos = DeltaY(ypos);
zpos = DeltaZ(zpos);

% store details in case requested
details.fit = m1 + m2;
details.xpos = xpos;
details.ypos = ypos;
details.zpos = zpos;
details.likelihood_x = likelihood_x;
details.likelihood_y = likelihood_y;
details.likelihood_z = likelihood_z;
details.likelihood_i = likelihood_i;
details.likelihood_j = likelihood_j;

details.DeltaX = DeltaX;
details.DeltaY = DeltaY;
details.DeltaZ = DeltaZ;

% update: compute expected value
p = exp(bsxfun(@minus,lX,max(lX,[],2)));
p = bsxfun(@rdivide,p,sum(p,2));
xpos = p * DeltaX';

p = exp(bsxfun(@minus,lY,max(lY,[],2)));
p = bsxfun(@rdivide,p,sum(p,2));
ypos = p * DeltaY';

p = exp(bsxfun(@minus,lZ,max(lZ,[],2)));
p = bsxfun(@rdivide,p,sum(p,2));
zpos = p * DeltaZ';

return;

subplot(211);
plot(t,xpos,t,ypos+3,t,zpos+6)

subplot(212);
mot2 = reshape(permute(mot(:,1:20,1:40,:),[2 1 3 4]),[40 40 size(mot,4)]);
for i = 1:size(mot2,3)-4; 
    cla; 
    imagesc(mean(mot2(:,:,i:i+3),3)); 
    hold on; 
    plot(gridSize/2-ypos(i,3),gridSize/2-xpos(i,3),'.w',3*gridSize/2-zpos(i,3),gridSize/2+xpos(i,3)/dx_di+ypos(i,3)/dy_di,'.k'); 
    drawnow; 
end

function [likelihood_i likelihood_j] = likelihoodToPos(L)
up = 3;
likelihood_i = zeros(size(L,3),0);
likelihood_j = zeros(size(L,3),0);
for i = 1:size(L,3)
    a = interp2(L(:,:,i),up,'cubic');
    likelihood_i(i,1:size(a,1)) = log(sum(exp(a - max(a(:))),2))' + max(a(:));
    likelihood_j(i,1:size(a,2)) = log(sum(exp(a - max(a(:))),1)) + max(a(:)); 
    if mod(i,1000) == 0 || i == size(L,3)
        disp(sprintf('Computing position marginals: [%u/%u]',i,size(L,3))); 
    end
    
end

function L = trackStack(stack,ref,x,y,up)

if nargin < 2 || isempty(ref)
    ref = mean(stack(:,:,round((end/2):(end/2)+100)),3);
end

if nargin < 3 || isempty(x), x = -3:3; end
if nargin < 4 || isempty(y), y = -3:3; end
if nargin < 5 || isempty(up), up = 1; end

% if need trim edge of stack to shift enough against ref
% assumes ref and stack size the same
if range(x) + size(stack,2) > size(ref,2)
    stack(:,1:-x(1),:) = [];
    stack(:,end-x(end)+1:end,:) = [];
end
if range(y) + size(stack,1) > size(ref,1)
    stack(1:-y(1),:,:) = [];
    stack(end-y(end)+1:end,:,:) = [];
end

% indicies of reference
idxX = 1-x(1):up:size(ref,2)-x(end);
idxY = 1-y(1):up:size(ref,1)-y(end);

% preallocate memory
L = zeros([length(y) length(x) size(stack,3)]);
for i = 1:length(x)
    for j = 1:length(y)
        r = ref(idxY + y(i), idxX + x(j));
        const(i,j) = sum(r(:));
        logGen(i,j,:) = log(r(:));
    end
end
for k = 1:size(stack,3)
    sample = stack(:,:,k);
    gl = sum(gammaln(1+sample(:)));
    for i = 1:length(x); 
        for j = 1:length(y); 
            L(i,j,k) = -const(i,j) - gl + reshape(logGen(i,j,:),1,[]) * sample(:);
        end
    end
    if mod(k,1000) == 0 || k == size(stack,3)
        disp(sprintf('Computing shift likelihood: [%u/%u]',k,size(stack,3))); 
    end
end

function L = trackShifts(stack)

x = -2:2;
y = -2:2;

idxX = 1-x(1):size(stack,2)-x(end);
idxY = 1-y(1):size(stack,1)-y(end);

for k = 1:size(stack,3)
    ref = mean(stack(:,:,max(1,k-10):max(1,k-1)),3);
    sample = stack(idxY,idxX,k);
    gl = gammaln(1+sample(:));
    for i = 1:length(x); 
        for j = 1:length(y);
            L(i,j,k) = imageLikelihood(ref(idxY + y(i), idxX + x(j)),sample,gl); 
        end
    end
    if mod(k,500) == 0 || k == size(stack,3)
        disp(sprintf('[%u/%u]',k,size(stack,3))); 
    end
end

function L = shiftedMotion(generative,sample,x,y)
% Compute the image likelihood at various shifts
sample = sample(1-y(1):end-y(end),1-x(1):end-x(end));

if x > 0
    generative = generative(:,1:end-x);
    sample = sample(:,x+1:end);
elseif x < 0
    generative = generative(:,1-x:end);
    sample = sample(:,1:end+x);
end

if y > 0
    generative = generative(1:end-y,:);
    sample = sample(1+y:end,:);
elseif y < 0
    generative = generative(1-y:end,:);
    sample = sample(1:end+y,:);
end

L = imageLikelihood(generative,sample);

function L = imageLikelihood(generative,sample,gl)
if nargin < 3, gl = gammaln(1+sample(:)); end
% Compute the likelihood of seeing a given image given the rates
L = sum(-generative(:) -gl  + sample(:) .* log(generative(:)));
