addpath ..


make_cmap

%minimum length of time series to use
minT = 5;

[sta,lat,lon,h,lab,t1,tm,dT,m,ngood,numsol,ve,vn,vu,se,sn,su,xe50,xn50,xu50,rho]=GetMIDASVelocities('IGS08',[]);

minlon = -127;
maxlon = -117;
maxlat = 50;
minlat = 40;


ind = dT>=minT & lon>minlon & lon<maxlon & lat<maxlat & lat>minlat & dT>=minT;

lat = lat(ind); lon=lon(ind); ve=ve(ind); vn=vn(ind); vu=vu(ind); se=se(ind); sn=sn(ind); su=su(ind); 
dT = dT(ind);

%median filter
 x=linspace(-127,-117,40);
 y=linspace(40,50,40);
 [X,Y]=meshgrid(x,y);
 V = griddata(lon,lat,vu,X,Y,'nearest');


Vf = medfilt2(V);
Vf_pts = interp2(X,Y,Vf,lon,lat,'nearest');




figure
worldmap([40 50],[-127+360 -117+360 ])
hold on
scatterm(lat,lon,50,vu,'fill')
load coastlines
geoshow(coastlat,coastlon)
colorbar
colormap(cmap)
caxis([-3 3])


figure
worldmap([40 50],[-127+360 -117+360 ])
hold on
scatterm(lat,lon,50,Vf_pts,'fill')
load coastlines
geoshow(coastlat,coastlon)
colorbar
colormap(cmap)
caxis([-3 3])




%northern profile
minlon = -127;
maxlon = -117;
maxlat = 48.2;
minlat = 47.3;


% %central profile
% minlon = -127;
% maxlon = -117;
% maxlat = 45;
% minlat = 44;

% %southern profile
% minlon = -127;
% maxlon = -117;
% maxlat = 43.5;
% minlat = 41.5;


ind = lon>minlon & lon<maxlon & lat<maxlat & lat>minlat & dT>=minT;
lat = lat(ind); lon=lon(ind); ve=ve(ind); vn=vn(ind); vu=vu(ind); se=se(ind); sn=sn(ind); su=su(ind); 
Vf_pts = Vf_pts(ind);

origin = [-122.5 46];
xy = llh2local([lon';lat'],origin)';

figure
worldmap([40 50],[-127+360 -117+360 ])
hold on
scatterm(lat,lon,50,vu,'fill')
geoshow(coastlat,coastlon)


%vertical color scale
vumin = -3;
vumax = 3;
NC=64;
z1 = linspace(0,1,NC/2)';
z1f = flipud(z1);
cmap = [ones(NC/2,1) z1 z1; 
        z1f          z1f  ones(NC/2,1)];
cmap = flipud(cmap);

colormap(cmap)
colorbar




%load leveling data
lev = load('leveling_Bruhat.txt');
xy_lev = llh2local([lev(:,1)';lev(:,2)'],origin)';



%reference
[y,ref] = max(xy(:,1));

vhoriz = sqrt( (ve-ve(ref)).^2 + (vn-vn(ref)).^2 );



%load slab
S = load('CascadiaSlabxyz.txt');


figure
worldmap([min(S(:,2)) max(S(:,2))],[-127+360 -117+360 ])
scatterm(S(:,2),S(:,1),10,S(:,3),'fill')
geoshow(coastlat,coastlon)
colorbar



ind = S(:,1)>minlon & S(:,1)<maxlon & S(:,2)<maxlat & S(:,2)>minlat;
S = S(ind,:);

xy_slab = llh2local([S(:,1)';S(:,2)'],origin)';


trench = min(xy_slab(~isnan(S(:,3)),1));
xmax = max(xy(:,1)-trench);

figure

subplot(311)
errorbar(xy(:,1)-trench,vhoriz,se,'o')
xlim([0 xmax])
grid on

subplot(312)
errorbar(xy(:,1)-trench,vu,su,'o')
hold on
errorbar(xy(:,1)-trench,Vf_pts,su,'r.')
xlim([0 xmax])
grid on

subplot(313)
plot(xy_slab(:,1)-trench,S(:,3),'.')
axis equal
xlim([0 xmax])
grid on


Xslab = xy_slab(:,1)-trench;
Dslab = S(:,3);

Xdata = xy(:,1)-trench;
Vh = vhoriz;
sigh = se;

Vu = Vf_pts;
sigu = su;








