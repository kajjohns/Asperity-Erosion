%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script builes a trianglur mesh of a surface defined
% by depth contours based on code of Yo Fukushima, 9 Jul 2010


% load contour file with columns: `lon, lat, height (km, negative numbers)
contr_file = 'hik_contours.txt';

% columns: lon, lat, rake (angle, degrees), rate (mm/yr)
% Note: rake is diretion of overriding plate relative to subducting plate,
% assuming right-hand rule (thumb is strike vector, index finger is dip
% vecteor, palm down)
conv_rate_file = 'hik_rake_rate.txt'

%origin to convert from llh to local Caresian coordinates
origin = [176, -40];  

%nominal side length of patches (km)
intv = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../tools
addpath ../ne_10m_coastline

contr = load(contr_file);

%% Durga -- comment this line out the first time
rakes_rates = load(conv_rate_file);


%% convert to local km
contrkm = llh2local([contr'],origin)';
%contrkm2 = llh2local([contr(:,2)'; contr(:,1)'],origin2)';
contrkm = [contrkm, contr(:,3)];


figure;scatter(contrkm(:,1),contrkm(:,2),10,contrkm(:,3),'filled'); colorbar; 
xlabel('Easting'); ylabel('Northing'); daspect([1,1,1]);
hold on
bbox = [170 -60; 180 -30]; plot_coast_xy(bbox,origin,'k')


[el,nd] = make_tri_mesh(contrkm,intv);


figure;
h = trisurf(el,nd(:,1),nd(:,2),nd(:,3)); colorbar; view(2); daspect([1,1,1]);
xlabel('Easting'); ylabel('Northing'); daspect([1,1,1]);
hold on
bbox = [170 -60; 180 -30]; plot_coast_xy(bbox,origin,'k')

%generate triangular mesh geometry parameters
patch_stuff=make_triangular_patch_stuff(el,nd);


%% Durga -- run to this line to compute triangular patch geometry. Then run other codes to make rakes and rates.

%interpolate rakes and rates to patch centers
xyrake = llh2local(rakes_rates(:,1:2)',origin)';
centroids = patch_stuff.centroids_faces;
rakes = griddata(xyrake(:,1),xyrake(:,2),rakes_rates(:,3),centroids(:,1),centroids(:,2));
rakes(isnan(rakes)) = griddata(xyrake(:,1),xyrake(:,2),rakes_rates(:,3),centroids(isnan(rakes),1),centroids(isnan(rakes),2),'nearest');

rates = griddata(xyrake(:,1),xyrake(:,2),rakes_rates(:,4),centroids(:,1),centroids(:,2));
rates(isnan(rates)) = griddata(xyrake(:,1),xyrake(:,2),rakes_rates(:,4),centroids(isnan(rates),1),centroids(isnan(rates),2),'nearest');


figure;
h = trisurf(el,nd(:,1),nd(:,2),nd(:,3),rakes); colorbar; view(2); daspect([1,1,1]);
xlabel('Easting'); ylabel('Northing'); daspect([1,1,1]);
hold on
bbox = [170 -60; 180 -30]; plot_coast_xy(bbox,origin,'k')
title('convergence rake')

figure;
h = trisurf(el,nd(:,1),nd(:,2),nd(:,3),rates); colorbar; view(2); daspect([1,1,1]);
xlabel('Easting'); ylabel('Northing'); daspect([1,1,1]);
hold on
bbox = [170 -60; 180 -30]; plot_coast_xy(bbox,origin,'k')
title('plate convergence rate')


figure;
h = trimesh(el,nd(:,1),nd(:,2),nd(:,3)); view(2); daspect([1,1,1]);
hold on
quiver3(centroids(:,1),centroids(:,2),centroids(:,3),patch_stuff.strikevec_faces(:,1),patch_stuff.strikevec_faces(:,2),patch_stuff.strikevec_faces(:,3))
quiver3(centroids(:,1),centroids(:,2),centroids(:,3),patch_stuff.dipvec_faces(:,1),patch_stuff.dipvec_faces(:,2),patch_stuff.dipvec_faces(:,3))
title('patch strike and dip vectors')