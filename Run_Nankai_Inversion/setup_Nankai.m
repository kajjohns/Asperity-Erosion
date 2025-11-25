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
%data_filename = 'hikurangi_elastic_nosmooth_hetgf.txt';


%path to build mesh script
%mesh_path = './build_mesh/make_mesh.m';
%load mesh file
load nankai_mesh_JAMSTEC_intv7.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath ../tools
addpath ../ne_10m_coastline


%build mesh
%need to interpolate rakes to mesh
patch_stuff = make_triangular_patch_stuff(el,nd);
centroids = patch_stuff.centroids_faces;
load Elizabeth_mesh
rakes = griddata(centroids_Eliz(:,1),centroids_Eliz(:,2),rake,centroids(:,1),centroids(:,2));
rakes(isnan(rakes)) = griddata(centroids_Eliz(:,1),centroids_Eliz(:,2),rake,centroids(isnan(rakes),1),centroids(isnan(rakes),2),'nearest');
rates = griddata(centroids_Eliz(:,1),centroids_Eliz(:,2),converge_rate,centroids(:,1),centroids(:,2));
rates(isnan(rates)) = griddata(centroids_Eliz(:,1),centroids_Eliz(:,2),converge_rate,centroids(isnan(rates),1),centroids(isnan(rates),2),'nearest');

%convert to mm/yr
rates = rates*1000;

%% load and plot data
load baseline_data
%convert from m/yr to mm/yr
Vbase = Vbase*1000;
Sigbase = Sigbase*1000;

max_lat = 38;
min_lat = 30;
max_lon = 138.7;
min_lon = 130.0;
bbox = [min_lon min_lat; max_lon max_lat];

%load vertical 
load interseismic_GPS_vertical
xy_vert = llh2local(llh',fliplr(origin))';
%toss out data outside of model region
ind = xy_vert(:,1)<300 & xy_vert(:,1)>-250 & xy_vert(:,2)<450;
xy_vert = xy_vert(ind,:);
Vu = Vu_filter(ind,:);


%plot baselines (if computed)
if ~invert_vel
    figure; hold on;  
    axis equal
    for k=1:size(BaseEnds,1)
        cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[Vbase(k)/L(k) Vbase(k)/L(k)])
    end
    title('Baseline elongation rates')
    load cmap
    colormap(cmap)
    colorbar
    caxis([-max(abs(Vbase./L)) max(abs(Vbase./L))])
    plot_coast_xy(bbox,fliplr(origin),'k')
end





%% compute displacement GFs or load file

if compute_disp  

    if invert_vel
        %velocity GFs
        [Ge,Gn,Gu] = make_dispG_triangular(el,nd,[xysites zeros(size(xysites,1),1)]',rakes,[]);
    else
        %baseline GFs
        Gbase = Get_Gs_baselines(el,nd,xy_hz,rakes,baselines,Vec_unit,crossind,L);
    end

else

    load(disp_filename)

end

Gbase = -Gbase; %change to backslip

%vertical
[Gedummy,Gndummy,Gu] = make_dispG_triangular(el,nd,[xy_vert zeros(size(xy_vert,1),1)]',rakes,[]);
Gu = -Gu;




%% Hmatrix setup
setup_Hmat   %this is a script

