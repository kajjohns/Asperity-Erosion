%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT 

%invert velocities or baselines?
%true for velocities, false for baselines
invert_vel = false; 

%compute stress H-matrix?
% false uses existing h-matrix
compute_hmat = false; 


%compute displacement Greens functions?
%false uses existing displacement GFs save in 
%matlab file with name disp_filename
compute_disp = true;
disp_filename = 'test_disp';


%path to build mesh script
mesh_path = './mesh/make_mesh.m';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath ../tools
addpath ../ne_10m_coastline/

%build mesh
run(mesh_path)


%% load and plot data
C = readtable('NZ_GNSS_comb_v9_sillremoved_17-Apr-2024.dat');
xysites_rollins = llh2local([C.lon';C.lat'],origin)';
%Ve = C.ve; Sige = C.se;
%Vn = C.vn; Sign = C.sn;
Vu = C.vu; Sigu = C.su;

load gnss_data.txt
xysites = llh2local(gnss_data(:,1:2)',origin)';
Ve = gnss_data(:,5);
Vn = gnss_data(:,6);
Sige = gnss_data(:,7);
Sign = gnss_data(:,8);



ind = xysites(:,2)<-300 | xysites(:,2)>700 | xysites(:,1)<-600 | xysites(:,1)>400;
xysites(ind,:) = [];
Ve(ind) = [];
Vn(ind) = [];
Sige(ind) = [];
Sign(ind) = [];

vert_xy = xysites_rollins;
ind = vert_xy(:,2)<-300 | vert_xy(:,2)>700 | vert_xy(:,1)<-600 | vert_xy(:,1)>400;
Vu(ind) = [];
Sigu(ind) = [];
vert_xy(ind,:) = [];

%remove nans

ind = isnan(Vu);
Vu(ind) = [];
Sigu(ind) = [];
vert_xy(ind,:) = [];

ind = (isnan(Vn) | isnan(Ve));
Ve(ind) = [];
Vn(ind) = [];
Sige(ind) = [];
Sign(ind) = [];
xysites(ind,:) = [];



%compute baselines if requested
if ~invert_vel

    %find baselines that cross elements breaking the ground surface (need
    %to integrate strain rates across these baselines)
    nodes_z = [nd(el(:,1),3) nd(el(:,2),3) nd(el(:,3),3)];  %depth of nodes in all triangles
    nodes_x = [nd(el(:,1),1) nd(el(:,2),1) nd(el(:,3),1)];  %x position of nodes in all triangles
    nodes_y = [nd(el(:,1),2) nd(el(:,2),2) nd(el(:,3),2)];  %y position of nodes in all triangles
    zero_nodes = nodes_z==0;
    surf_break = sum(zero_nodes,2)==2; %find triangles with two nodes at surface
    %mak patch endpoints for elements that break surface
    PatchEnds = nan(size(el,1),4);
    for j=1:size(PatchEnds,1)
        if surf_break(j)
            PatchEnds(j,[1 3]) = nodes_x(j,zero_nodes(j,:));
            PatchEnds(j,[2 4]) = nodes_y(j,zero_nodes(j,:));
        end
    end
        
    [Vbase,Sigbase,BaseEnds,L,baselines,Vec_unit,crossind] = make_baseline_rate_changes(xysites,Ve,Vn,Sige,Sign,PatchEnds);

end

%clean up baselines -- remove long ones
ind = L>50;
Vbase(ind) = [];
Sigbase(ind) = [];
BaseEnds(ind,:) = [];
L(ind) = [];
baselines(ind,:) = [];
Vec_unit(ind,:) = [];

%load NZNSHM slip deficit strain rates (crustal faults) and convert to
%baselines.  Remove this contribution from the observed baselines.

load backslip_strainrates_NZ
obs_coords_xy = llh2local(obs_coords_llh(:,1:2)',origin)';
Vbase_crustal = compute_baselines_from_strainrates(xysites,obs_coords_xy,mean(Exx_bs,2),mean(Exy_bs,2),mean(Eyy_bs,2),Vec_unit,baselines,L);

Vbase = Vbase - Vbase_crustal';


%plot velocities
figure; hold on;  
quiver(xysites(:,1),xysites(:,2),Ve,Vn)
axis equal
title('Observed velocities')
bbox = [170 -60; 180 -30]; plot_coast_xy(bbox,origin,'k')


%plot baselines (if computed)
if ~invert_vel
    figure; hold on;  
    axis equal
    for k=1:size(BaseEnds,1)
        cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[Vbase(k)/L(k) Vbase(k)/L(k)])
    end
    title('Baseline elongation rates')
    bbox = [170 -60; 180 -30]; plot_coast_xy(bbox,origin,'k')
    load cmap
    colormap(cmap)
    colorbar
    caxis([-max(abs(Vbase./L)) max(abs(Vbase./L))])
end





%% compute displacement GFs or load file

if compute_disp  

    if invert_vel
        %velocity GFs
        [Ge,Gn,Gu] = make_dispG_triangular(el,nd,[xysites zeros(size(xysites,1),1)]',rakes,[]);
    else
        %baseline GFs
        Gbase = Get_Gs_baselines(el,nd,xysites,rakes,baselines,Vec_unit,crossind,L);
    end

else

    load(disp_filename)

end

%vertical
[Gedummy,Gndummy,Gu] = make_dispG_triangular(el,nd,[vert_xy zeros(size(vert_xy,1),1)]',rakes,[]);





%% Hmatrix setup
setup_Hmat   %this is a script

