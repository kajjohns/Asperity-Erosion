function [sta,lat,lon,h,lab,t1,tm,dT,m,ngood,numsol,ve,vn,vu,se,sn,su,xe50,xn50,xu50,rho]=GetMIDASVelocities(frame,velfile)
% [sta,lat,lon,h,lab,t1,tm,dT,m,ngood,numsol,ve,vn,vu,se,sn,su,xe50,xn50,xu50,rho]=GetMIDASVelocities(frame,velfile)
%
% frame = 'IGS14'
%
% if velfile is empty gets from http://geodesy.unr.edu
% otherwise uses file. 
% 



if isempty(velfile)
    velfile = ['https://geodesy.unr.edu/velocities/midas.' char(frame) '.txt'];
    [S,stat]=urlread(velfile);

 

    if stat==1
        C=textscan(S,'%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    else
        error('Could not get midas velocity file from http://geodesy.unr.edu');
    end
else
    fid=fopen(velfile,'r');
    C=textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    fclose(fid);
end

sta=C{1};
lab=C{2};
t1=C{3};
tm=C{4};
dT=C{5};
m=C{6};
ngood=C{7};
numsol=C{8};
ve=1000*C{9};
vn=1000*C{10};
vu=1000*C{11};
se=1000*C{12};
sn=1000*C{13};
su=1000*C{14};
xe50=C{15};
xn50=C{16};
xu50=C{17};
% what are columns 18-24?
rho = ngood./(dT*365);




llhfile='https://geodesy.unr.edu/NGLStationPages/llh.out';

[L,stat]=urlread(llhfile);
if stat==1
    C=textscan(L,'%s %f %f %f');
    stnllh=C{1};
    latllh=C{2};
    lonllh=C{3};
    hllh=C{4};
else
    error('Could not get station coordinates.');
end

% hdir=GetHomePath;
% [stnllh,latllh,lonllh,hllh]=...
%     Readllh([chomp(hdir) '/Science/XYZfiles/tsdata/llh/llh.all']);
[~,jl,il]=intersect(sta,stnllh);

sta=sta(jl);
lab=lab(jl);
lon = lonllh(il);
lat = latllh(il);
h=hllh(il);
t1=t1(jl);
tm=tm(jl);
dT=dT(jl);
m=m(jl);
ngood=ngood(jl);
numsol=numsol(jl);
ve=ve(jl);
vn=vn(jl);
vu=vu(jl);
se=se(jl);
sn=sn(jl);
su=su(jl);
xe50=xe50(jl);
xn50=xn50(jl);
xu50=xu50(jl);
rho=rho(jl);

lon(lon<-180)=lon(lon<-180)+360;



