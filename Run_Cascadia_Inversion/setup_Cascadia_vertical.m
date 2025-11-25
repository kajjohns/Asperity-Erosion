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

%name of file containing data
%format of columns: lon, lat, Ve(mm/yr), Vn(mm/yr), Sige, Sign
data_filename = 'cascadia_horizontal_gps_240412.txt';
vertical_data_filename = 'Cascadia_vertical_homogenized.txt';


%path to build mesh script
mesh_path = './make_mesh.m';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath ../tools
addpath ../ne_10m_coastline

%build mesh
run(mesh_path)


%% load and plot data


%load GPS data -- convert to mm/yr
data = load(data_filename);
xysites = llh2local(data(:,1:2)',origin)';
Ve = data(:,3)*1000;
Vn = data(:,4)*1000;
if size(data,2)==6
    Sige = data(:,5)*1000;
    Sign = data(:,6)*1000;

else  %if no uncertainties, assign value of 2
    Sige = 2*ones(size(Ve));
    Sign = 2*ones(size(Vn));
end



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


%plot velocities
figure; hold on;  
quiver(xysites(:,1),xysites(:,2),Ve,Vn)
axis equal
title('Observed velocities')
plot_coast_xy(bbox,origin,'k')


%plot baselines (if computed)
if ~invert_vel
    figure; hold on;  
    axis equal
    for k=1:size(BaseEnds,1)
        cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[Vbase(k)/L(k) Vbase(k)/L(k)])
    end
    title('Baseline elongation rates')
    plot_coast_xy(bbox,origin,'k')
    load cmap
    colormap(cmap)
    colorbar
    caxis([-max(abs(Vbase./L)) max(abs(Vbase./L))])
end


drawnow


%% compute displacement GFs or load file

if compute_disp  

  
        %velocity GFs
        [Ge,Gn,Gu] = make_dispG_triangular(el,nd,[xysites zeros(size(xysites,1),1)]',rakes,[]);
    
        %rakes are in wrong direction for GFs - switch
        Ge = -Ge;
        Gn = -Gn;
        Gu = -Gu;
    
    
        %baseline GFs
        Gbase = Get_Gs_baselines(el,nd,xysites,rakes,baselines,Vec_unit,crossind);
    
        %rakes are in wrong direction for GFs - switch
        Gbase = -Gbase;

    

else

    load(disp_filename)

end






%% Hmatrix setup
setup_Hmat   %this is a script
patch_stuff = make_triangular_patch_stuff(el,nd);
centroids = patch_stuff.centroids_faces;

%vertical
vert = load(vertical_data_filename);
xysites_vert = llh2local(vert(:,1:2)',origin)';
ind = xysites_vert(:,1)<400;
xysites_vert = xysites_vert(ind,:);
Vu = vert(ind,3); %convert from m to mm
Sigu = vert(ind,4);

[Ge_dummy,Gn_dummy,Gu] = make_dispG_triangular(el,nd,[xysites_vert zeros(size(xysites_vert,1),1)]',rakes,[]);
Gu = -Gu;   

