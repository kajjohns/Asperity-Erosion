%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build_mesh.m
%
% Construct a triangular mesh of a subduction interface from depth contours
% and plate convergence information, and save the geometry for use in
% forward modeling and inversions of eroding asperities.
%
% INPUT FILES (expected in ./data or on MATLAB path):
%   - contr_file: subduction interface contours
%         'Cascadia_contours.txt'
%         columns: lon (deg), lat (deg), depth (km, negative)
%
%   - conv_rate_file: plate convergence rake and rate
%         'Cascadia_rake_rate.txt'
%         columns: lon (deg), lat (deg),
%                  rake (deg; direction of overriding plate wrt subducting
%                        plate, right-hand rule),
%                  rate (mm/yr)
%
%   - origin.txt:
%         1×3 row with reference lon, lat, elevation used by llh2local
%         to define the local Cartesian coordinate system (km).
%
% USER PARAMETERS:
%   max_depth : maximum depth (km; positive scalar) of the mesh (e.g., 46)
%   intv      : nominal side length (km) of triangular patches (e.g., 8)
%   bbox      : [lon_min lat_min; lon_max lat_max] bounding box for plots
%
% DEPENDENCIES (on MATLAB path):
%   - tools/llh2local.m
%   - tools/make_tri_mesh.m
%   - tools/make_triangular_patch_stuff.m
%   - tools/plot_coast_xy.m
%   - ne_10m_coastline/* (coastline data)
%
% MAIN STEPS:
%   1) Load depth contours and convert (lon,lat,depth) -> local Cartesian.
%   2) Remove contours deeper than max_depth.
%   3) Build a triangular mesh of the interface.
%   4) Compute patch geometry (centroids, strike/dip vectors, etc.).
%   5) Interpolate convergence rake and rate onto patch centroids.
%   6) Plot contours, mesh, rake, rate, strike/dip vectors.
%   7) Save all mesh products to mesh_file.mat.
%
% OUTPUT:
%   mesh_file.mat containing:
%       el, nd, patch_stuff, rakes, rates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------
%  USER-SPECIFIED INPUT FILES AND PARAMETERS
% -------------------------------------------------------------------------

% File: lon, lat, depth (km, negative)
contr_file = 'Cascadia_contours.txt';

% Maximum depth to include in the mesh (km)
max_depth = 46;

% File: lon, lat, rake (deg), convergence rate (mm/yr)
conv_rate_file = 'Cascadia_rake_rate.txt';

% Nominal side length (km) for triangular patches
intv = 8;

% Bounding box for coastline plot [lon_min lat_min; lon_max lat_max]
bbox = [-128 38; -120 54];

%% ------------------------------------------------------------------------
%  ADD PATHS TO SUPPORTING FUNCTIONS AND DATA
% -------------------------------------------------------------------------
addpath tools
addpath ne_10m_coastline
addpath data

%% ------------------------------------------------------------------------
%  LOAD ORIGIN FOR LOCAL COORDINATE SYSTEM
% -------------------------------------------------------------------------
load origin.txt   % [lon0 lat0 elev0], used by llh2local

%% ------------------------------------------------------------------------
%  LOAD DEPTH CONTOURS
% -------------------------------------------------------------------------
contr = load(contr_file);   % [lon, lat, depth]

% Convert to local Cartesian coordinates (km)
% llh2local returns [x; y; z] in km relative to origin
contrkm = llh2local([contr'], origin)';
contrkm = [contrkm, contr(:,3)];  % append depth column

% Keep only contours shallower than max_depth
contrkm = contrkm(contrkm(:,3) > -max_depth, :);

%% ------------------------------------------------------------------------
%  PLOT DEPTH CONTOURS
% -------------------------------------------------------------------------
figure
scatter(contrkm(:,1), contrkm(:,2), 10, contrkm(:,3), 'filled');
colorbar
xlabel('Easting (km)'); ylabel('Northing (km)');
daspect([1 1 1]);
hold on
plot_coast_xy(bbox, origin, 'k')
title('Subduction interface depth contours')
xlim([min(contrkm(:,1)) max(contrkm(:,1))])
ylim([min(contrkm(:,2)) max(contrkm(:,2))])

%% ------------------------------------------------------------------------
%  CONSTRUCT TRIANGULAR MESH FROM DEPTH CONTOURS
% -------------------------------------------------------------------------
[el, nd] = make_tri_mesh(contrkm, intv);
%  nd: node coordinates (x,y,depth)
%  el: element connectivity (triangles)

%% ------------------------------------------------------------------------
%  PLOT TRIANGULAR MESH
% -------------------------------------------------------------------------
figure;
trisurf(el, nd(:,1), nd(:,2), nd(:,3), 'edgecolor', [.5 .5 .5]);
view(2); colorbar; daspect([1 1 1]);
xlabel('Easting'); ylabel('Northing');
hold on
plot_coast_xy(bbox,origin,'k')
title('Triangular mesh of subduction interface')
xlim([min(contrkm(:,1)) max(contrkm(:,1))])
ylim([min(contrkm(:,2)) max(contrkm(:,2))])

%% ------------------------------------------------------------------------
%  COMPUTE GEOMETRIC INFORMATION FOR EACH TRIANGLE
% -------------------------------------------------------------------------
patch_stuff = make_triangular_patch_stuff(el, nd);
% provides centroids, normals, strike/dip vectors, etc.

centroids = patch_stuff.centroids_faces;

%% ------------------------------------------------------------------------
%  LOAD CONVERGENCE RAKES & RATES AND INTERPOLATE TO PATCH CENTROIDS
% -------------------------------------------------------------------------
% Columns: lon, lat, rake, rate
rakes_rates = load(conv_rate_file);

% Convert station locations to local Cartesian
xyrake = llh2local(rakes_rates(:,1:2)', origin)';

% Interpolate rake
rakes = griddata(xyrake(:,1), xyrake(:,2), ...
                 rakes_rates(:,3), centroids(:,1), centroids(:,2));

% Fill NaNs via nearest neighbor
rakes(isnan(rakes)) = griddata(xyrake(:,1), xyrake(:,2), ...
                rakes_rates(:,3), centroids(isnan(rakes),1), ...
                centroids(isnan(rakes),2), 'nearest');

% Interpolate convergence rate (mm/yr)
rates = griddata(xyrake(:,1), xyrake(:,2), ...
                 rakes_rates(:,4), centroids(:,1), centroids(:,2));

% Fill NaNs via nearest neighbor
rates(isnan(rates)) = griddata(xyrake(:,1), xyrake(:,2), ...
                rakes_rates(:,4), centroids(isnan(rates),1), ...
                centroids(isnan(rates),2), 'nearest');

%% ------------------------------------------------------------------------
%  PLOT INTERPOLATED FIELDS (RAKE AND RATE)
% -------------------------------------------------------------------------
figure;
trisurf(el, nd(:,1), nd(:,2), nd(:,3), rakes, 'edgecolor', 'none');
view(2); colorbar; daspect([1 1 1]);
xlabel('Easting'); ylabel('Northing');
hold on; plot_coast_xy(bbox,origin,'k');
title('Interpolated convergence rake')

figure;
trisurf(el, nd(:,1), nd(:,2), nd(:,3), rates, 'edgecolor', 'none');
view(2); colorbar; daspect([1 1 1]);
xlabel('Easting'); ylabel('Northing');
hold on; plot_coast_xy(bbox,origin,'k');
title('Interpolated plate convergence rate (mm/yr)')

%% ------------------------------------------------------------------------
%  PLOT STRIKE AND DIP VECTORS FOR QA/QC
% -------------------------------------------------------------------------
figure;
trimesh(el, nd(:,1), nd(:,2), nd(:,3), 'edgecolor', 'none');
view(2); daspect([1 1 1]); hold on;

% Plot strike and dip vectors at centroids
quiver3(centroids(:,1), centroids(:,2), centroids(:,3), ...
        patch_stuff.strikevec_faces(:,1), ...
        patch_stuff.strikevec_faces(:,2), ...
        patch_stuff.strikevec_faces(:,3))

quiver3(centroids(:,1), centroids(:,2), centroids(:,3), ...
        patch_stuff.dipvec_faces(:,1), ...
        patch_stuff.dipvec_faces(:,2), ...
        patch_stuff.dipvec_faces(:,3))

title('Patch strike and dip vectors')
xlim([min(contrkm(:,1)) max(contrkm(:,1))])
ylim([min(contrkm(:,2)) max(contrkm(:,2))])

%% ------------------------------------------------------------------------
%  SAVE MESH AND GEOMETRIC INFORMATION
% -------------------------------------------------------------------------
save mesh_file el nd patch_stuff rakes rates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
