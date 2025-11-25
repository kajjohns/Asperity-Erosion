
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build_GFs.m
%
% Compute (or load) elastic half-space Green’s functions (GFs) for a
% triangular mesh of a subduction interface, for use in asperity inversions.
% Optionally:
%   - Use station *velocities* or *baseline* length changes as data.
%   - Compute displacement Green’s functions (surface displacements per unit
%     slip on each triangular patch).
%   - Set up (or reuse) an H-matrix representation of stress kernels.
%
% Typical workflow:
%   1) Run build_mesh.m to generate:
%        mesh_file.mat -> el, nd, patch_stuff, rakes, rates
%   2) Run this script to:
%        - load geodetic data (horizontal and vertical)
%        - compute or load displacement GFs (Ge, Gn, Gu, Gbase)
%        - set up H-matrix (via setup_Hmat script)
%        - save everything into save_filename.mat
%
% USER OPTIONS:
%   invert_vel    : if true, invert GPS velocities; if false, invert baselines.
%   compute_hmat  : if true, compute H-matrix (slow, one-time per mesh).
%                   if false, existing H-matrix in ./hmat/ is reused.
%   compute_disp  : if true, compute displacement GFs and save to disk.
%                   if false, load an existing disp_filename.mat.
%   disp_filename : name of .mat file with existing displacement GFs.
%   horizontal_filename : text file with horizontal GPS data
%         columns: lon, lat, Ve, Vn, Sige, Sign
%   is_hz_meter   : true if Ve/Vn/Sige/Sign are in m/yr (convert to mm/yr).
%   vert_filename : text file with vertical GPS data
%         columns: lon, lat, Vu, Sigu
%                   (empty '' if no vertical data)
%   is_vert_meter : true if Vu/Sigu are in m/yr (convert to mm/yr).
%   bbox          : [lon_min lat_min; lon_max lat_max] for plotting coastline.
%   save_filename : base name of output .mat file for model setup.
%
% INPUT FILES (expected on path):
%   - mesh_file.mat: triangular mesh and rake/rate (from build_mesh.m)
%   - ./data/origin.txt: reference point for llh2local (lon, lat, elev).
%   - GPS data files in ./data (or on path).
%
% DEPENDENCIES (on MATLAB path):
%   - tools/llh2local.m
%   - tools/Nikkhoo/* (triangular dislocation / displacement kernels)
%   - tools/make_baseline_rate_changes.m
%   - tools/Get_Gs_baselines.m
%   - tools/make_dispG_triangular.m
%   - tools/plot_coast_xy.m
%   - setup_Hmat.m (script for H-matrix construction)
%   - ne_10m_coastline/* (coastline data)
%   - cline.m, cmap.mat (for colored baseline plots)
%
% OUTPUT:
%   save_filename.mat containing (at minimum):
%       el, nd, patch_stuff, rakes, rates
%       Ge, Gn, Gu, Gbase (if computed)
%       data, xysites, (and vertical/baseline quantities if used)
%   plus H-matrix-related variables if setup_Hmat is run.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------
%  USER CONFIGURATION: WHAT TO COMPUTE / LOAD
% -------------------------------------------------------------------------

% Invert velocities (true) or baselines (false)?
invert_vel = false;

% Compute stress H-matrix?
%   true  = compute and save (slow, but only once per mesh)
%   false = use existing H-matrix saved in 'hmat' folder
compute_hmat = false;

% Compute displacement Green’s functions?
%   true  = compute GFs and (optionally) save them
%   false = load GFs from disp_filename
compute_disp = true;
disp_filename = 'test_disp';

% Data files (in ./data or on MATLAB path)

% Horizontal GPS data:
%   columns: lon, lat, Ve(mm/yr), Vn(mm/yr), Sige(mm/yr), Sign(mm/yr)
horizontal_filename = 'cascadia_horizontal_gps_240412.txt';
is_hz_meter = true;   % true if in m/yr and need conversion to mm/yr

% Vertical GPS data:
%   columns: lon, lat, Vu(mm/yr), Sigu(mm/yr)
%   If none, provide empty string: vert_filename = '';
vert_filename = 'Cascadia_vertical_homogenized.txt';
is_vert_meter = false;  % true if in m/yr and need conversion to mm/yr

% Bounding box for plotting coastline
%   first row: lon,lat of bottom left corner
%   second row: lon,lat of top right corner
bbox = [-128 38; -120 54];

% Base name for output .mat file containing setup results
save_filename = 'setup_cascadia';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  END OF CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------
%  ADD PATHS TO TOOLS AND DATA
% -------------------------------------------------------------------------
addpath tools
addpath tools/Nikkhoo/
addpath data
addpath ne_10m_coastline

%% ------------------------------------------------------------------------
%  LOAD MESH AND ORIGIN
% -------------------------------------------------------------------------
load mesh_file.mat         % loads el, nd, patch_stuff, rakes, rates
load ./data/origin.txt     % reference point for llh2local

%% ------------------------------------------------------------------------
%  LOAD AND PREPARE HORIZONTAL GPS DATA
% -------------------------------------------------------------------------
% Load horizontal GPS site locations and velocities
data = load(horizontal_filename);
% Convert lon/lat to local Cartesian (km)
xysites = llh2local(data(:,1:2)', origin)';
% Extract east, north velocities and uncertainties
Ve   = data(:,3);
Vn   = data(:,4);
Sige = data(:,5);
Sign = data(:,6);

% Convert units from m/yr to mm/yr if requested
if is_hz_meter
    Ve   = Ve   * 1000;
    Vn   = Vn   * 1000;
    Sige = Sige * 1000;
    Sign = Sign * 1000;
end

%% ------------------------------------------------------------------------
%  LOAD AND PREPARE VERTICAL GPS DATA (OPTIONAL)
% -------------------------------------------------------------------------
if ~isempty(vert_filename)
    vert = load(vert_filename);
    xysites_vert = llh2local(vert(:,1:2)', origin)';  % local (x,y)
    Vu   = vert(:,3);  % vertical rate
    Sigu = vert(:,4);  % uncertainty

    % Convert units from m/yr to mm/yr if requested
    if is_vert_meter
        Vu   = Vu   * 1000;
        Sigu = Sigu * 1000;
    end
end

%% ------------------------------------------------------------------------
%  CONSTRUCT BASELINES IF INVERTING BASELINE DATA
% -------------------------------------------------------------------------
if ~invert_vel

    % For baselines that cross elements breaking the surface, we must
    % integrate along the baseline, accounting for discontinuities.
    %
    % Identify triangles that break the ground surface (two nodes at z=0).
    nodes_z = [nd(el(:,1),3) nd(el(:,2),3) nd(el(:,3),3)];  % depths at nodes
    nodes_x = [nd(el(:,1),1) nd(el(:,2),1) nd(el(:,3),1)];  % x positions
    nodes_y = [nd(el(:,1),2) nd(el(:,2),2) nd(el(:,3),2)];  % y positions

    zero_nodes = (nodes_z == 0);              % nodes on the surface
    surf_break = sum(zero_nodes,2) == 2;      % triangles with 2 surface nodes

    % Build segment endpoints for elements that break the surface
    PatchEnds = nan(size(el,1), 4);  % [x1 y1 x2 y2] for each triangle
    for j = 1:size(PatchEnds,1)
        if surf_break(j)
            PatchEnds(j,[1 3]) = nodes_x(j,zero_nodes(j,:));
            PatchEnds(j,[2 4]) = nodes_y(j,zero_nodes(j,:));
        end
    end

    % Construct baseline observables (length changes) and their geometry
    [Vbase, Sigbase, BaseEnds, L, baselines, Vec_unit, crossind] = ...
        make_baseline_rate_changes(xysites, Ve, Vn, Sige, Sign, PatchEnds);

end

%% ------------------------------------------------------------------------
%  PLOT OBSERVED HORIZONTAL VELOCITIES
% -------------------------------------------------------------------------
figure; hold on;
quiver(xysites(:,1), xysites(:,2), Ve, Vn)
axis equal
title('Observed horizontal velocities')
plot_coast_xy(bbox, origin, 'k')
xlim([min(xysites(:,1)) max(xysites(:,1))])
ylim([min(xysites(:,2)) max(xysites(:,2))])

%% ------------------------------------------------------------------------
%  PLOT BASELINE ELONGATION RATES (IF USING BASELINES)
% -------------------------------------------------------------------------
if ~invert_vel
    figure; hold on;
    axis equal
    % Plot each baseline colored by elongation rate
    for k = 1:size(BaseEnds,1)
        cline([BaseEnds(k,1) BaseEnds(k,3)], ...
              [BaseEnds(k,2) BaseEnds(k,4)], ...
              [Vbase(k)/L(k) Vbase(k)/L(k)]);
    end
    title('Observed baseline elongation rates')
    plot_coast_xy(bbox, origin, 'k')
    xlim([min(xysites(:,1)) max(xysites(:,1))])
    ylim([min(xysites(:,2)) max(xysites(:,2))])

    load cmap
    colormap(cmap)
    colorbar
    caxis([-max(abs(Vbase./L)) max(abs(Vbase./L))])
end

drawnow

%% ------------------------------------------------------------------------
%  COMPUTE OR LOAD DISPLACEMENT GREEN’S FUNCTIONS
% -------------------------------------------------------------------------
if compute_disp

    %----------------------------------------------------------------------
    % Velocity GFs at horizontal GPS sites
    %   Ge, Gn, Gu: [n_sites x n_patches]
    %   Each column = displacement due to 1 mm/yr slip on a patch
    %----------------------------------------------------------------------

    [Ge, Gn, Gu] = make_dispG_triangular( ...
        el, nd, [xysites zeros(size(xysites,1),1)]', rakes, []);

    % Rakes are defined opposite to the convention in GFs -> flip sign
    Ge = -Ge;
    Gn = -Gn;
    % Gu is computed here but will be recomputed for vertical sites below

    %----------------------------------------------------------------------
    % Baseline GFs (if using baselines)
    %   Gbase: [n_baselines x n_patches] elongation rate per unit slip
    %----------------------------------------------------------------------
    Gbase = Get_Gs_baselines(el, nd, xysites, rakes, ...
                             baselines, Vec_unit, crossind, L);

    % Flip sign to match rake convention
    Gbase = -Gbase;

    %----------------------------------------------------------------------
    % Vertical velocity GFs at vertical GPS sites (if provided)
    %   Gu: [n_vert_sites x n_patches]
    %----------------------------------------------------------------------
    if ~isempty(vert_filename)
        [~, ~, Gu] = make_dispG_triangular( ...
            el, nd, [xysites_vert zeros(size(xysites_vert,1),1)]', rakes, []);

        % Flip sign to match rake convention
        Gu = -Gu;
    end

    % (You may want to save Ge, Gn, Gu, Gbase to disp_filename here, if desired.)
    % save(disp_filename,'Ge','Gn','Gu','Gbase','-v7.3');  % example

else
    % Load precomputed displacement Green’s functions
    load(disp_filename)
end

%% ------------------------------------------------------------------------
%  QUICK QA/QC: CHECK SIGNS OF DISPLACEMENTS (UNIFORM BACKSLIP)
% -------------------------------------------------------------------------

% Approximate predicted horizontal velocities for uniform slip
figure; hold on;
quiver(xysites(:,1), xysites(:,2), sum(Ge,2), sum(Gn,2))
axis equal
title('Computed velocities for uniform backslip (1 mm/yr)')
plot_coast_xy(bbox, origin, 'k')
xlim([min(xysites(:,1)) max(xysites(:,1))])
ylim([min(xysites(:,2)) max(xysites(:,2))])

% Vertical velocities at vertical sites (requires vertical data)
figure; hold on;
scatter(xysites_vert(:,1), xysites_vert(:,2), 50, sum(Gu,2), 'fill')
load cmap
colormap(cmap)
colorbar
caxis([-0.2 0.2])
axis equal
title('Computed vertical velocities for 1 mm/yr uniform backslip')
plot_coast_xy(bbox, origin, 'k')
xlim([min(xysites(:,1)) max(xysites(:,1))])
ylim([min(xysites(:,2)) max(xysites(:,2))])

% Baseline elongation rates under uniform backslip (if using baselines)
figure; hold on;
Vb = sum(Gbase,2);  % predicted baseline elongation rates
axis equal
for k = 1:size(BaseEnds,1)
    cline([BaseEnds(k,1) BaseEnds(k,3)], ...
          [BaseEnds(k,2) BaseEnds(k,4)], ...
          [Vb(k)/L(k) Vb(k)/L(k)]);
end
title('Baseline elongation rates for uniform 1 mm/yr backslip')
plot_coast_xy(bbox, origin, 'k')
xlim([min(xysites(:,1)) max(xysites(:,1))])
ylim([min(xysites(:,2)) max(xysites(:,2))])

load cmap
colormap(cmap)
colorbar
caxis([-max(abs(Vb./L)) max(abs(Vb./L))])

%% ------------------------------------------------------------------------
%  HMATRIX SETUP (STRESS GREEN’S FUNCTIONS, ETC.)
% -------------------------------------------------------------------------
% setup_Hmat is a script that uses the current workspace (mesh, etc.) to
% construct and/or load an H-matrix representation of stress kernels.
setup_Hmat   % this is a script and may take significant time if compute_hmat=true

%% ------------------------------------------------------------------------
%  SAVE COMPLETE SETUP FOR INVERSION
% -------------------------------------------------------------------------
% Save all relevant variables for later inversion scripts
save(save_filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%