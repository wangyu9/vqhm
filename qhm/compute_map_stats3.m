function [stats] = compute_map_stats3(V,F,UV)

u = UV(:,1);
v = UV(:,2);

stats = struct();

G = intrinsic_grad(V,F); 
Area = abs(doublearea(V,F))/2;

% n = size(UV,1);
f = size(G,1)/2;

Gx = G(1:f,:);
Gy = G(f+1:2*f,:);

Js = zeros(f,2,2);

Js(:,1,1) = Gx * u;
Js(:,1,2) = Gy * u;

Js(:,2,1) = Gx * v;
Js(:,2,2) = Gy * v;

stats.sigmas = zeros(f,2);
stats.sigma_ratio = zeros(f,1);

for j=1:f
    Jj = squeeze(Js(j,:,:));
    [uu,ss,vv] = svd(Jj);
    ss = diag(ss);
    stats.sigmas(j,:) = ss;
    stats.sigma_ratio(j) = ss(1)./ss(2);
end

stats.MIPS_old = (stats.sigma_ratio + 1./stats.sigma_ratio);

stats.wavg_MIPS = (stats.sigma_ratio + 1./stats.sigma_ratio) .* Area / sum(Area);


flipped = doublearea([u,v], F) < 0;

stats.num_flips = numel(find(flipped));

stats.Area = Area;

% stats.get_wavg_ARAP = @()  
% mean(sum((stats.sigmas-1).^2,2))

stats.arap = sum((stats.sigmas-1).^2,2);

% mean(   max([stats.sigmas(:,1),1./stats.sigmas(:,2)],[],2)  )

stats.iso = max([stats.sigmas(:,1),1./stats.sigmas(:,2)],[],2); 
