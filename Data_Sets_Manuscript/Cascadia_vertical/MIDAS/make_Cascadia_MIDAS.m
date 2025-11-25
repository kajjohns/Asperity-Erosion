make_cmap

%minimum length of time series to use
minT = 5;

[sta,lat,lon,h,lab,t1,tm,dT,m,ngood,numsol,ve,vn,vu,se,sn,su,xe50,xn50,xu50,rho]=GetMIDASVelocities('IGS14',[]);

minlon = -127;
maxlon = -115;
maxlat = 50;
minlat = 40;


ind = dT>=minT & lon>minlon & lon<maxlon & lat<maxlat & lat>minlat & dT>=minT;

lat = lat(ind); lon=lon(ind); ve=ve(ind); vn=vn(ind); vu=vu(ind); se=se(ind); sn=sn(ind); su=su(ind); 
dT = dT(ind);

%median filter
 x=linspace(-127,-115,40);
 y=linspace(40,50,40);
 [X,Y]=meshgrid(x,y);
 V = griddata(lon,lat,vu,X,Y,'nearest');


Vf = medfilt2(V);
Vf_pts = interp2(X,Y,Vf,lon,lat,'nearest');




figure
worldmap([40 50],[-127+360 -115+360 ])
hold on
scatterm(lat,lon,50,vu,'fill')
load coastlines
geoshow(coastlat,coastlon)
colorbar
colormap(cmap)
caxis([-3 3])


figure
worldmap([40 50],[-127+360 -115+360 ])
hold on
scatterm(lat,lon,50,Vf_pts,'fill')
load coastlines
geoshow(coastlat,coastlon)
colorbar
colormap(cmap)
caxis([-3 3])




plot_dem_Cascadia
hold on
scatter(lon,lat,50,vu,'fill')
colormap(cmap)
caxis([-6 6])
colorbar


plot_dem_Cascadia
hold on
scatter(lon,lat,50,Vf_pts,'fill')
colormap(cmap)
caxis([-6 6])
colorbar


