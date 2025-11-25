function i_lock = get_locked_radii_indices2(node_centroid,ptsDZ,radii,centroids,dips,minrad)


thetas = linspace(0,2*pi,50);
i_lock = false(size(centroids,1),1);
%i_ring = i_lock;

for k=1:size(ptsDZ,1)

    if radii(k)>minrad

    ind = dsearchn([centroids(:,2) -centroids(:,3)],ptsDZ(k,:));

    asp = [ptsDZ(k,1)+radii(k)*cos(thetas)' ptsDZ(k,2)+radii(k)*tan(dips(ind)*pi/180)*sin(thetas)'];
    locked = inpoly([centroids(:,2) -centroids(:,3)],asp);

    %check that horizontal distance is also within radius (for flat areas
    %of interface)
    j = node_centroid(k);
    asp_xy = [centroids(j,1)+radii(k)*cos(thetas)' centroids(j,2)+radii(k)*sin(thetas)'];
    locked_xy = inpoly([centroids(:,1) centroids(:,2)],asp_xy);

    locked = locked & locked_xy;

    i_lock = i_lock | locked; 

    %ring= [ptsDZ(k,1)+(radii(k)+D)*cos(thetas) ptsDZ(k,2)+(radii(k)+D)*sin(thetas)];
    %i_ring = inpoly([centroids(:,2) -centroids(:,3)],ring) & ~i_lock;

    end

end





