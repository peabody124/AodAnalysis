function n = aod_fastoopsi(traces,time)

V.dt = mean(diff(time));
V.Ncells = 1;
V.T = length(time);
V.est_b = 1;
V.est_gam = 1;
V.est_a = 1;
V.est_lam = 1;
V.est_sig = 1;

P.gam = 1-mean(diff(time))/1.5;

for i = 1:size(traces,2)
    % initially want to estimate baseline, can't use poissson
    V.fast_iter_max = 3;
    V.fast_poiss = 0;

    F = traces(:,i);
    F = F / mean(F);
    P.sig = std(diff(F))/2;
    P.b = min(F);
    [n(:,i) P_best V_best C]  = fast_oopsi(F,V,P);
    
    % parse reconstruction (hacky)
    n(n(:,i) < 0.04,i) = 0;
    C = filter(1,[1, -P_best.gam],n(:,i));
    
    bl = F - C;        
    h = hamming(round(30 / mean(diff(time))));
    h = h / sum(h);    
    bl_est = imfilter(bl,h,'replicate');
    bl_est = bl_est - max(bl_est);    

    V.fast_iter_max = 1;
    V.fast_poiss = 1;
    %P = P_best;
    P.b = .9*median(F-bl_est);    
    [n(:,i) P_best V_best C]  = fast_oopsi(F,V,P);
    
%     subplot(311);
%     plot(time,F,time,C);
%     subplot(312);
%     plot(time,n);
%     subplot(313);
%     imagesc(n);
%     drawnow

end