function i_lock = get_locked_radii_indices(ptsDZ,radii,centroids,dips,minrad)


% Precompute
thetas = linspace(0,2*pi,50);
cos_t = cos(thetas(:));
sin_t = sin(thetas(:));
tan_dips = tan(dips*pi/180);

% Build KDTree once
persistent Mdl
if isempty(Mdl)
    Mdl = KDTreeSearcher([centroids(:,2), -centroids(:,3)]);
end

% Initialize locked vector
i_lock = false(size(centroids,1),1);

for k = 1:size(ptsDZ,1)
    if radii(k) > minrad
        % Nearest centroid
        ind = knnsearch(Mdl, ptsDZ(k,:));

        % Polygon
        X = ptsDZ(k,1) + radii(k)*cos_t;
        Y = ptsDZ(k,2) + tan_dips(ind)*radii(k)*sin_t;

        % Vectorized inpolygon
        locked = inpolygon(centroids(:,2), -centroids(:,3), X, Y);
        i_lock = i_lock | locked;
    end
end



thetas = linspace(0,2*pi,50);
i_lock = false(size(centroids,1),1);
%i_ring = i_lock;


for k=1:size(ptsDZ,1)

    if radii(k)>minrad

    ind = dsearchn([centroids(:,2) -centroids(:,3)],ptsDZ(k,:));

    asp = [ptsDZ(k,1)+radii(k)*cos(thetas)' ptsDZ(k,2)+radii(k)*tan(dips(ind)*pi/180)*sin(thetas)'];
    locked = inpoly([centroids(:,2) -centroids(:,3)],asp);

    i_lock = i_lock | locked; 

    %ring= [ptsDZ(k,1)+(radii(k)+D)*cos(thetas) ptsDZ(k,2)+(radii(k)+D)*sin(thetas)];
    %i_ring = inpoly([centroids(:,2) -centroids(:,3)],ring) & ~i_lock;

    end

end





