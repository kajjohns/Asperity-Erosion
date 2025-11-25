function ring_tau = get_ring_tau(ring_tau_asp, ptsDZ, D, radii, centroids, i_lock, dips, minrad)



thetas = linspace(0,2*pi,50);
ring_tau = zeros(size(centroids,1),1);

for k=1:size(ptsDZ,1)

    if radii(k)>minrad

    ind = dsearchn([centroids(:,2) -centroids(:,3)],ptsDZ(k,:));


    ring= [ptsDZ(k,1)+(radii(k)+D)*cos(thetas)' ptsDZ(k,2)+(radii(k)+D)*tan(dips(ind)*pi/180)*sin(thetas)'];
    i_ring = inpoly([centroids(:,2) -centroids(:,3)],ring) & ~i_lock;

    ring_tau(i_ring) = ring_tau_asp(k);

    end

end





