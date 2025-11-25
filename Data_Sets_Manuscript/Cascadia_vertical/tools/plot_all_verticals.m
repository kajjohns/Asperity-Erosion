
%% download midas vels
addpath MIDAS
%minimum length of time series to use
minT = 5;
[sta,lat,lon,h,lab,t1,tm,dT,m,ngood,numsol,ve,vn,vu,se,sn,su,xe50,xn50,xu50,rho]=GetMIDASVelocities('IGS14',[]);

minlon = -127;
maxlon = -115;
maxlat = 50;
minlat = 40;

ind = dT>=minT & lon>minlon & lon<maxlon & lat<maxlat & lat>minlat & dT>=minT;
lat = lat(ind); lon=lon(ind); ve=ve(ind); vn=vn(ind); vu=vu(ind); se=se(ind); sn=sn(ind); su=su(ind); 

midas_vels = [lon lat vu su];

%% load other data sets
krogstad = load('./Krogstad_ScienceDirect_files/Krogstad_leveling.txt');
newton = load('./Newton_Supplemental/Newton_Data_S2.txt');
burgett = load("Burgett_2009.txt");

%GIA models of Yousefi et al. 2021
yousefi = load('./Yousefi_2021/all_models.txt');

%
%load Lau netcdf variables into a structure array
fn = './Lau/CONUS_gps_verticals.nc';
ni = ncinfo(fn);
for i=1:length(ni.Variables)
    vn = ni.Variables(i).Name;
    lau.(vn) = ncread(fn, vn);  % The result is a structure 
end


addpath ./dem

%create dem 
lat_range = [38 50];
lon_range = [-128 -115];
[A,lats,lons] = plot_global_etopo1(lat_range,lon_range);


%plot shaded hillslope dem
figure
[h,I,z]=dem(lons,lats,A,'LatLon','Legend','LandColor',.8*ones(359,3),'SeaColor',.6*ones(359,3),'Zlim',[-2000 2000],'Contrast',1);
hold on

%add state borders
bbox_ll = [-128 -115  38 50];
plot_state_borders_llh_box('k',bbox_ll)
xlim(lon_range)
ylim(lat_range)


hold on
scatter(midas_vels(:,1),midas_vels(:,2),30,midas_vels(:,3),'fill')
scatter(krogstad(:,2),krogstad(:,1),30,krogstad(:,3),'fill','d')
scatter(newton(:,1),newton(:,2),60,newton(:,3))
scatter(burgett(:,2),burgett(:,1),30,burgett(:,6),'fill','s')

caxis([-4 4])
colorbar
make_cmap
colormap(cmap)
title('all data, no alignment, no corrections')


%% Yousefi predicted

for k=1:11

%plot shaded hillslope dem
figure
[h,I,z]=dem(lons,lats,A,'LatLon','Legend','LandColor',.8*ones(359,3),'SeaColor',.6*ones(359,3),'Zlim',[-2000 2000],'Contrast',1);
hold on

%add state borders
bbox_ll = [-128 -115  38 50];
plot_state_borders_llh_box('k',bbox_ll)
xlim(lon_range)
ylim(lat_range)


hold on
scatter(yousefi(:,1)-360,yousefi(:,2),30,yousefi(:,2+k),'fill')
caxis([-4 4])
colorbar
colormap(cmap)

end

%yousefi average

%plot shaded hillslope dem
figure
[h,I,z]=dem(lons,lats,A,'LatLon','Legend','LandColor',.8*ones(359,3),'SeaColor',.6*ones(359,3),'Zlim',[-2000 2000],'Contrast',1);
hold on

%add state borders
bbox_ll = [-128 -115  38 50];
plot_state_borders_llh_box('k',bbox_ll)
xlim(lon_range)
ylim(lat_range)

mean_yousefi = mean(yousefi(:,3:end),2);

hold on
scatter(yousefi(:,1)-360,yousefi(:,2),30,mean_yousefi,'fill')
caxis([-4 4])
colorbar
colormap(cmap)
title('Average of Yousefi models')

%plot shaded hillslope dem
figure
[h,I,z]=dem(lons,lats,A,'LatLon','Legend','LandColor',.8*ones(359,3),'SeaColor',.6*ones(359,3),'Zlim',[-2000 2000],'Contrast',1);
hold on

%add state borders
bbox_ll = [-128 -115  38 50];
plot_state_borders_llh_box('k',bbox_ll)
xlim(lon_range)
ylim(lat_range)


hold on
scatter(yousefi(:,1)-360,yousefi(:,2),30,std(yousefi(:,3:end),[],2),'fill')
caxis([-4 4])
colorbar
colormap(cmap)
title('STD of Yousefi models')

%% Lau contributions
%plot shaded hillslope dem
figure
[h,I,z]=dem(lons,lats,A,'LatLon','Legend','LandColor',.8*ones(359,3),'SeaColor',.6*ones(359,3),'Zlim',[-2000 2000],'Contrast',1);
hold on

%add state borders
bbox_ll = [-128 -115  38 50];
plot_state_borders_llh_box('k',bbox_ll)
xlim(lon_range)
ylim(lat_range)


hold on
scatter(lau.lon(:),lau.lat(:),30,lau.net_vu(:),'fill')
caxis([-4 4])
colorbar
colormap(cmap)
title('Lau total (gps + all corrections)')


figure
[h,I,z]=dem(lons,lats,A,'LatLon','Legend','LandColor',.8*ones(359,3),'SeaColor',.6*ones(359,3),'Zlim',[-2000 2000],'Contrast',1);
hold on

%add state borders
bbox_ll = [-128 -115  38 50];
plot_state_borders_llh_box('k',bbox_ll)
xlim(lon_range)
ylim(lat_range)

lau_correction = lau.gps_vu_smooth(:)-lau.net_vu_smooth(:);

hold on
scatter(lau.lon(:),lau.lat(:),30,lau_correction,'fill')
caxis([-4 4])
colorbar
colormap(cmap)
title('Lau, all corrections')




%% combine and filter
%shift Newton and Burgett velocities into alignment with MIDAS 
all_vels = midas_vels(:,1:3);

%NOTE: Krogstat velocities are including in Newton
%all_vels = [all_vels; krogstad(:,2) krogstad(:,1) krogstad(:,3)];

%Bring Newton vels into alignment with MIDAS
for k=1:size(newton,1)

    dist =  sqrt((all_vels(:,1)-newton(k,1)).^2 + (all_vels(:,2)-newton(k,2)).^2 );
    [y,i] = min(dist);
    newton_map_to_GPS(k) = all_vels(i,3);

end

adjust_n = mean(newton(:,3)-newton_map_to_GPS');


%need to bring Burgett into alignment with newton and midas -- burgetter
%velocities are higher
for k=1:size(burgett,1)

    dist =  sqrt((all_vels(:,1)-burgett(k,2)).^2 + (all_vels(:,2)-burgett(k,1)).^2 );
    [y,i] = min(dist);
    burgett_map_to_GPS(k) = all_vels(i,3);

end

adjust_b = mean(burgett(:,6)-burgett_map_to_GPS');



all_vels = [all_vels; newton(:,1:3)-adjust_n];
all_vels = [all_vels; burgett(:,2) burgett(:,1) burgett(:,6)-adjust_b];

all_sigs = [midas_vels(:,4); newton(:,4);burgett(:,7)];



% %median filter -- at each observation point, take median of all values
% %within specified radius
% filter_radius = 0.25;  %degress
% 
% for k=1:size(all_vels,1)
% 
%     dist = sqrt((all_vels(:,1)-all_vels(k,1)).^2 + (all_vels(:,2)-all_vels(k,2)).^2 );
%     filt_vel(k) = median(all_vels(dist<filter_radius,3));
% 
% end

%outlier detection -- at each observation point compute  take median of all values
%within specified radius, take difference of vels and median, replace any velocities excreeding a threshold
%with the median

filter_radius = 0.25;  %degress
threshold = 2;

for k=1:size(all_vels,1)

    dist = sqrt((all_vels(:,1)-all_vels(k,1)).^2 + (all_vels(:,2)-all_vels(k,2)).^2 );
    med = median(all_vels(dist<filter_radius,3));

    if abs(all_vels(k,3)-med) > threshold
        filt_vel(k) = med;
    else
        filt_vel(k) = all_vels(k,3);
    end
end



figure
[h,I,z]=dem(lons,lats,A,'LatLon','Legend','LandColor',.8*ones(359,3),'SeaColor',.6*ones(359,3),'Zlim',[-2000 2000],'Contrast',1);
hold on

%add state borders
bbox_ll = [-128 -115  38 50];
plot_state_borders_llh_box('k',bbox_ll)
xlim(lon_range)
ylim(lat_range)


hold on
scatter(all_vels(:,1),all_vels(:,2),30,filt_vel,'fill')
caxis([-4 4])
colorbar
colormap(cmap)
title(['Data, median filter, radius = ' num2str(filter_radius) ' degrees'])


%subtract GIA corrections

mean_yousefi_interp = griddata(yousefi(:,1)-360,yousefi(:,2),mean_yousefi,all_vels(:,1),all_vels(:,2),'nearest');
lau_correction_interp = griddata(lau.lon(:),lau.lat(:),lau_correction,all_vels(:,1),all_vels(:,2),'nearest');




figure
[h,I,z]=dem(lons,lats,A,'LatLon','Legend','LandColor',.8*ones(359,3),'SeaColor',.6*ones(359,3),'Zlim',[-2000 2000],'Contrast',1);
hold on

%add state borders
bbox_ll = [-128 -115  38 50];
plot_state_borders_llh_box('k',bbox_ll)
xlim(lon_range)
ylim(lat_range)

hold on
scatter(all_vels(:,1),all_vels(:,2),30,filt_vel'-mean_yousefi_interp,'fill')
caxis([-4 4])
colorbar
colormap(cmap)
title(['Data minus mean of Yousefi models'])





figure
[h,I,z]=dem(lons,lats,A,'LatLon','Legend','LandColor',.8*ones(359,3),'SeaColor',.6*ones(359,3),'Zlim',[-2000 2000],'Contrast',1);
hold on

%add state borders
bbox_ll = [-128 -115  38 50];
plot_state_borders_llh_box('k',bbox_ll)
xlim(lon_range)
ylim(lat_range)

hold on
scatter(all_vels(:,1),all_vels(:,2),30,filt_vel'-lau_correction_interp,'fill')
caxis([-4 4])
colorbar
colormap(cmap)
title(['Data minus Lau correction'])


vert = filt_vel'-lau_correction_interp;

%% Plot profiles
lats = [41 43 45 47];
figure
for k=1:length(lats)

    subplot(2,2,k)

    ind = abs((all_vels(:,2)-lats(k)))<0.15;

    %errorbar(all_vels(ind,1),vert(ind),all_sigs(ind),'o')
     plot(all_vels(ind,1),vert(ind),'o')
    xlabel('longitude')
    ylabel('vertical velocity, mm/yr')
    title(['Latitude = ' num2str(lats(k))])

    ylim([-5 5])
    grid on

end



   


