
A = load('etopo1_Cascadia.xyz');


Nlon=800;
Nlat=800;


%lons = linspace(min(A(:,1)),max(A(:,1)),Nlon);
%lats = linspace(min(A(:,2)),max(A(:,2)),Nlat);

lons = linspace(-128,-117.5,Nlon);
lats = linspace(40,50,Nlat);




[LON,LAT] = meshgrid(lons,lats);
Z = griddata(A(:,1),A(:,2),A(:,3),LON,LAT);

figure
%[h,I,z]=dem(lons,lats,Z,'LatLon','Legend','LandColor',.5*ones(359,3),'SeaColor',.25*ones(359,3),'Zlim',[000 2000],'Contrast',1);
%[h,I,z]=dem(lons,lats,Z,'LatLon','Legend','LandColor',.5*ones(359,3),'SeaColor',1*ones(359,3),'Zlim',[000 2000],'Contrast',1,'Azimuth',45);
[h,I,z]=dem(lons,lats,Z,'LatLon','Legend','LandColor',.5*ones(359,3),'SeaColor',.7*ones(359,3),'Zlim',[-2000 2000],'Contrast',1,'Azimuth',45);

