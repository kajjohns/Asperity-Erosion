function [i_lock, ring_tau] = get_locked_and_ring_tau(...
        ptsDZ, radii, ring_tau_asp, D, centroids, minrad, ...
        cos_t, sin_t, tan_dips, KDTree)

nAsp   = numel(radii);
nTheta = numel(cos_t);
nCent  = size(centroids,1);

i_lock   = false(nCent,1);
ring_tau = zeros(nCent,1);

% Nearest dip index for all asperities
inds = knnsearch(KDTree, ptsDZ);

% Filter asperities above min radius
valid = radii > minrad;
ptsDZ = ptsDZ(valid,:);
radii = radii(valid);
inds  = inds(valid);
ring_tau_asp = ring_tau_asp(valid);
nAsp = numel(radii);

% Build polygons in bulk (nAsp × nTheta)
aspX  = ptsDZ(:,1) + radii.*cos_t';                     
aspY  = ptsDZ(:,2) + radii.*tan_dips(inds).*sin_t';     
ringX = ptsDZ(:,1) + (radii+D).*cos_t';                 
ringY = ptsDZ(:,2) + (radii+D).*tan_dips(inds).*sin_t'; 

% Loop only over polygons for inpolygon
for k = 1:nAsp
    % Locked asperity polygon
    insideAsp = inpolygon(centroids(:,2), -centroids(:,3), aspX(k,:), aspY(k,:));
    i_lock = i_lock | insideAsp;
    
    % Ring polygon
    insideRing = inpolygon(centroids(:,2), -centroids(:,3), ringX(k,:), ringY(k,:));
    ring_tau(insideRing & ~i_lock) = ring_tau_asp(k);
end
