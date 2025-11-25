function [A,lats,lons] = plot_global_etopo1(lat_range,lon_range)

[A,R] = geotiffread('ETOPO1_Ice_g_geotiff.tif');
lons = linspace(R.XWorldLimits(1),R.XWorldLimits(2),size(A,2));
lats = fliplr(linspace(R.YWorldLimits(1),R.YWorldLimits(2),size(A,1)));


%specify lat and lon range for map
%lat_range = [35 37];
%lon_range = [-119 -117];


i = lats<lat_range(2) & lats>lat_range(1);
j = lons>lon_range(1) & lons<lon_range(2);

lats = lats(i);
lons = lons(j);

A = A(i,:);
A = A(:,j);



%plot dem
%figure
%[h,I,z]=dem(lons,lats,A,'LatLon','Legend','LandColor',.8*ones(359,3),'SeaColor',.6*ones(359,3),'Zlim',[-2000 2000],'Contrast',1);
