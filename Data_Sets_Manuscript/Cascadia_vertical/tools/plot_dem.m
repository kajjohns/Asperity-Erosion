addpath ./dem


%create dem 
lat_range = [30 50];
lon_range = [-129 -101];
[A,lats,lons] = plot_global_etopo1(lat_range,lon_range);


%plot shaded hillslope dem
figure
[h,I,z]=dem(lons,lats,A,'LatLon','Legend','LandColor',.8*ones(359,3),'SeaColor',.6*ones(359,3),'Zlim',[-2000 2000],'Contrast',1);
hold on

%add state borders
bbox_ll = [-125 -108  32 50];
plot_state_borders_llh_box('k',bbox_ll)
