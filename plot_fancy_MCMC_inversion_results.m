%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_fancy_mcmc_inversion_results.m
%
% “Fancy” post-processing of MCMC inversion results:
%   - Overlays inferred coupling, locking probability, and process-zone
%     stressing rates on shaded-relief topography (ETOPO1).
%   - Plots observed and modeled baseline elongation rates and vertical
%     velocities on a DEM base.
%   - Examines misfit vs. distance from the trench and computes variance
%     reduction metrics.
%   - Computes total and “locked-only” moment rates and corresponding
%     equivalent 500-year moment magnitudes (Mw).
%   - Plots histograms of moment rates with a second x-axis in Mw.
%
% NOTE: Before you can make the DEM plots you must obtain ETOPO1_Ice_g_geotiff.tif from: 
% https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/georeferenced_tiff/?utm_source=chatgpt.com
%
% REQUIRED VARIABLES / INPUTS (from previous scripts/workspace):
%   From build_mesh.m / build_GFs.m / mcmc_inversion.m / plot_mcmc_inversion_results.m:
%     - el, nd, patch_stuff, strikes, srate, bbox, origin
%     - ptsDZ, centroids_rot, minrad, D
%     - baselines, Vbase, L, xysites, xysites_vert, Vu
%     - vert_filename (to know if vertical data exist)
%
%   From MCMC output folder (folder_name):
%     - logrho.txt        : log-likelihood trace
%     - M_radii.txt       : asperity radii samples
%     - M_ring_tau.txt    : ring stressing rate samples
%     - dhat.txt          : predicted data for each sample
%     - locked_index.txt  : locked patch indices/flags per sample
%     - creep_rates.txt   : creep rates per sample
%
% DEPENDENCIES:
%   - ./tools/dem/, ../dem: DEM and shaded-relief plotting (dem, plot_global_etopo1)
%   - tools/slanCM/slanCM.m: perceptual colormaps
%   - cmap.mat           : custom red-blue colormap
%   - get_locked_radii_indices.m, get_ring_tau.m
%   - local2llh.m, plot_coast.m, plot_coast_xy.m, cline.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add DEM tools for shaded-relief base maps
addpath ./tools/dem/

% Folder containing MCMC chain output
folder_name = 'Cascadia_outputs';

%----------------------------------------------------------------------
% Load and inspect log-likelihood to visualize convergence
%----------------------------------------------------------------------
eval(['load ./' folder_name '/logrho.txt'])
figure
plot(logrho)
xlabel('Sample index')
ylabel('Log-likelihood')
title('Log-likelihood trace')

%----------------------------------------------------------------------
% Discard burn-in samples
%----------------------------------------------------------------------
discard = 50;   % number of initial samples to remove as burn-in

% Load chain outputs
eval(['load ./' folder_name '/M_radii.txt'])
eval(['load ./' folder_name '/M_ring_tau.txt'])
eval(['load ./' folder_name '/dhat.txt'])
eval(['load ./' folder_name '/locked_index.txt'])
eval(['load ./' folder_name '/creep_rates.txt'])

% Remove burn-in from all chains
M_radii(1:discard,:)      = [];
M_ring_tau(1:discard,:)   = [];
dhat(1:discard,:)         = [];
locked_index(1:discard,:) = [];
creep_rates(1:discard,:)  = [];

%----------------------------------------------------------------------
% Reconstruct process-zone stressing rates Ring_Taus for each sample
%----------------------------------------------------------------------
clear Ring_Taus
for k = 1:size(M_radii,1)

    radii     = M_radii(k,:)';
    ring_taus = M_ring_tau(k,:)';

    % Determine locked patches for these radii
    i_locked = get_locked_radii_indices( ...
        ptsDZ, radii, centroids_rot, patch_stuff.dip_faces, minrad);

    % Accumulated process-zone stress on patches for this sample
    Ring_Taus(:,k) = get_ring_tau( ...
        ring_taus, ptsDZ, D, radii, centroids_rot, ...
        i_locked, patch_stuff.dip_faces, minrad);

end

%----------------------------------------------------------------------
% Colormap for coupling (custom blend of Blues + Inferno)
%----------------------------------------------------------------------
c1 = colormap(flipud(slanCM('inferno')));
c2 = colormap(flipud(slanCM('Blues')));
mycmap = [c2(1:2:end,:); c1(1:2:end,:)];

%----------------------------------------------------------------------
% Rotate interface back to original orientation (if needed downstream)
%----------------------------------------------------------------------
str = mean(strikes)*pi/180;  % mean strike in radians
R   = [cos(str) -sin(str); ...
       sin(str)  cos(str)];
centroids_rot      = (R * centroids(:,1:2)')';
centroids_rot(:,3) = centroids(:,3);
centroids_rot(:,2) = centroids_rot(:,2) - min(centroids_rot(:,2));

%----------------------------------------------------------------------
% Convert mesh nodes to geographic (lon,lat)
%----------------------------------------------------------------------
nd_ll = local2llh(nd(:,1:2)', origin)';

% Load red-blue colormap used for several plots
load cmap

% DEM utilities (if also located one directory up)
addpath ../dem

% Light azimuth for shaded-relief
clear alpha
LA = 45;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 1: Coupling ratio over shaded-relief topography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lat_range = [40 50];
lon_range = [-128 -121];

% Load ETOPO1 subset
[A,lats,lons] = plot_global_etopo1(lat_range, lon_range);

% DEM base map with shaded relief
hfig = figure;
[h,I,z] = dem(lons, lats, A, 'LatLon', ...
    'LandColor', .8*ones(359,3), ...
    'SeaColor',  .95*ones(359,3), ...
    'Zlim',      [-2000 2000], ...
    'Contrast',  1, ...
    'Azimuth',   LA);

hold on

% Coupling ratio: 1 - (mean creep rate / plate rate)
h = trisurf(el, nd_ll(:,1), nd_ll(:,2), nd(:,3)+50, ...
            1 - mean(creep_rates,1)' ./ (srate*1000), ...
            'edgecolor','none');
colorbar
title('Coupling Ratio')
caxis([-1 1])
colormap(mycmap)

% Coastline and map limits
plot_coast(bbox)
xlim(lon_range)
ylim(lat_range)

% Transparent overlay so DEM shows through
alpha(h,0.5)

ylim(lat_range)
xlim(lon_range)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 2: Probability of locking on shaded-relief topography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lat_range = [40 50];
lon_range = [-128 -121];
[A,lats,lons] = plot_global_etopo1(lat_range, lon_range);

hfig = figure;
[h,I,z] = dem(lons, lats, A, 'LatLon', ...
    'LandColor', .8*ones(359,3), ...
    'SeaColor',  .9*ones(359,3), ...
    'Zlim',      [-2000 2000], ...
    'Contrast',  1, ...
    'Azimuth',   LA);
hold on

h = trisurf(el, nd_ll(:,1), nd_ll(:,2), nd(:,3)+50, ...
            mean(locked_index), 'edgecolor','none');
colormap(flipud(slanCM('heat')));
colorbar
caxis([0 1])
title('Probability of Locking')

plot_coast(bbox)
xlim(lon_range)
ylim(lat_range)
alpha(h,0.7)

ylim(lat_range)
xlim(lon_range)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 3: Mean ring stressing rate on shaded-relief topography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lat_range = [40 50];
lon_range = [-128 -121];
[A,lats,lons] = plot_global_etopo1(lat_range, lon_range);

hfig = figure;
[h,I,z] = dem(lons, lats, A, 'LatLon', ...
    'LandColor', .8*ones(359,3), ...
    'SeaColor',  .9*ones(359,3), ...
    'Zlim',      [-2000 2000], ...
    'Contrast',  1, ...
    'Azimuth',   LA);
hold on

h = trisurf(el, nd_ll(:,1), nd_ll(:,2), nd(:,3)+50, ...
            mean(Ring_Taus,2), 'edgecolor','none');
colormap(flipud(slanCM('heat')));
colorbar
title('Mean Ring Tau')

plot_coast(bbox)
xlim(lon_range)
ylim(lat_range)
alpha(h,0.7)

ylim(lat_range)
xlim(lon_range)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 4: Baseline elongation rates (observed vs. model) on DEM base
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Split dhat into baseline and vertical components
Nbase    = size(baselines,1);
dhat_base = dhat(:,1:Nbase);
dhat_u    = dhat(:,1+Nbase:end);

% Convert station and baseline endpoints to geographic coordinates
sites_ll = local2llh(xysites',      origin)';
vert_ll  = local2llh(xysites_vert', origin)';

BaseEnds_ll = [ ...
    local2llh(BaseEnds(:,1:2)', origin)' ...
    local2llh(BaseEnds(:,3:4)', origin)'];

lat_range = [40 50];
lon_range = [-127 -119];
[A,lats,lons] = plot_global_etopo1(lat_range, lon_range);

% Observed baseline elongation rates
hfig = figure;
[h,I,z] = dem(lons, lats, A, 'LatLon', ...
    'LandColor', .8*ones(359,3), ...
    'SeaColor',  .9*ones(359,3), ...
    'Zlim',      [-2000 2000], ...
    'Contrast',  1, ...
    'Azimuth',   LA);
hold on

for k = 1:size(BaseEnds,1)
    cline([BaseEnds_ll(k,1) BaseEnds_ll(k,3)], ...
          [BaseEnds_ll(k,2) BaseEnds_ll(k,4)], ...
          [50 50], ...
          [Vbase(k)/L(k) Vbase(k)/L(k)]);
end

title('Observed baseline elongation rates')
colormap(cmap)
colorbar
caxis([-0.1 0.1])
plot_coast(bbox)
ylim(lat_range)
xlim(lon_range)

% Model baseline elongation rates (posterior mean)
hfig = figure;
[h,I,z] = dem(lons, lats, A, 'LatLon', ...
    'LandColor', .8*ones(359,3), ...
    'SeaColor',  .9*ones(359,3), ...
    'Zlim',      [-2000 2000], ...
    'Contrast',  1, ...
    'Azimuth',   LA);
hold on

v = mean(dhat_base,1);
for k = 1:size(BaseEnds,1)
    cline([BaseEnds_ll(k,1) BaseEnds_ll(k,3)], ...
          [BaseEnds_ll(k,2) BaseEnds_ll(k,4)], ...
          [50 50], ...
          [v(k)/L(k) v(k)/L(k)]);
end

title('Model baseline elongation rates')
colormap(cmap)
colorbar
caxis([-0.1 0.1])
plot_coast(bbox)
ylim(lat_range)
xlim(lon_range)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 5: Vertical rates (observed vs. model) on DEM base
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Observed vertical rates
hfig = figure;
[h,I,z] = dem(lons, lats, A, 'LatLon', ...
    'LandColor', .8*ones(359,3), ...
    'SeaColor',  .9*ones(359,3), ...
    'Zlim',      [-2000 2000], ...
    'Contrast',  1, ...
    'Azimuth',   LA);
hold on
scatter(vert_ll(:,1), vert_ll(:,2), 50, Vu, 'fill')

title('Observed vertical rates')
colormap(cmap)
colorbar
caxis([-5 5])
plot_coast(bbox)
ylim(lat_range)
xlim(lon_range)

% Model vertical rates (posterior mean)
hfig = figure;
[h,I,z] = dem(lons, lats, A, 'LatLon', ...
    'LandColor', .8*ones(359,3), ...
    'SeaColor',  .9*ones(359,3), ...
    'Zlim',      [-2000 2000], ...
    'Contrast',  1, ...
    'Azimuth',   LA);
hold on

vu = mean(dhat_u,1);
scatter(vert_ll(:,1), vert_ll(:,2), 50, vu, 'fill')

title('Model vertical rates')
colormap(cmap)
colorbar
caxis([-5 5])
plot_coast(bbox)
ylim(lat_range)
xlim(lon_range)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Misfit vs distance from trench (baselines and verticals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Approximate trench trace: nodes at maximum depth
i      = (nd(:,3) == max(nd(:,3)));
trench = nd(i,1:2);   % trench polyline in local coordinates

% Baseline centers (midpoints in local coordinates)
BaseCenters = [ ...
    mean(BaseEnds(:,[1,3]),2) ...
    mean(BaseEnds(:,[2,4]),2)];

% Distance from each baseline center to trench
[k,dist] = dsearchn(trench, BaseCenters);

% Plot observed and modeled baseline elongation rates vs distance
figure
plot(dist, abs(Vbase)./L, 'bo')
hold on
v = mean(dhat_base,1);
plot(dist, abs(v')./L, 'ro')
title('Fit to baseline elongation rates vs distance from trench')
xlabel('Distance from trench (km)')
ylabel('Elongation rate (1/yr)')
legend('Observed','Model','Location','best')

% Distance from vertical GPS sites to trench
[k,dist_u] = dsearchn(trench, xysites_vert);

figure
plot(dist_u, Vu, 'bo')
hold on
plot(dist_u, vu, 'ro')
title('Fit to vertical rates vs distance from trench')
xlabel('Distance from trench (km)')
ylabel('Vertical rate (mm/yr)')
legend('Observed','Model','Location','best')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variance reduction for baselines and verticals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

variance_reduction_baselines = 1 - norm(Vbase - v').^2 / norm(Vbase).^2
variance_reduction_vertical  = 1 - norm(Vu - vu').^2   / norm(Vu).^2

% Variance reduction within 300 km of the trench
i = dist < 300;
variance_reduction_baselines_300 = ...
    1 - norm(Vbase(i) - v(i)').^2 / norm(Vbase(i)).^2

i = dist_u < 300;
variance_reduction_vertical_300 = ...
    1 - norm(Vu(i) - vu(i)').^2 / norm(Vu(i)).^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute moment rates and equivalent 500-year Mw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch areas (km^2 -> m^2)
A = patch_stuff.area_faces * 1e6;

% Mean slip deficit rate (m/yr) from posterior mean creep_rates
sr = (srate) - mean(creep_rates)'/1000;   % srate is in m/yr, creep_rates in mm/yr

% Total moment rate (N·m/yr) assuming shear modulus mu = 30 GPa
Mo_total = 30e9 * sum(sr .* A);

% Equivalent Mw for 500-year accumulation (convert to dyne-cm first)
Mw_500_total = (2/3) * log10(500 * Mo_total * 1e7) - 10.7

% Locking-probability-weighted moment rate
locking_prob = mean(locked_index)';
Mo_locking   = 30e9 * sum(sr .* A .* locking_prob);
Mw_500_locking = (2/3) * log10(500 * Mo_locking * 1e7) - 10.7

%----------------------------------------------------------------------
% Histograms of total and “locked-only” moment rates from all samples
%----------------------------------------------------------------------

% sr: [n_patches x n_samples] slip deficit rates (m/yr)
sr = repmat(srate, 1, size(creep_rates,1)) - creep_rates'/1000;

% Moment rate for each sample (N·m/yr)
Mo_total   = 30e9 * sum(sr .* repmat(A,1,size(creep_rates,1)));

% Locking-weighted moment rate for each sample
Mo_locking = 30e9 * sum(sr .* repmat(A,1,size(creep_rates,1)) .* locked_index');

figure
histogram(Mo_total,   'Normalization','pdf')
hold on
histogram(Mo_locking, 'Normalization','pdf')

xlabel('Moment rate Mo [N·m/yr]')
ylabel('Probability density')
legend('total','locked areas')
set(gca,'FontSize',15)

%----------------------------------------------------------------------
% Add second x-axis on top for Mw (equivalent 500-year Mw)
%----------------------------------------------------------------------

T = 500;  % time interval in years
% Convert moment rate (N·m/yr) to Mw for T-year accumulation:
%   Mw = (2/3)(log10(M0[dyn·cm]) - 10.7), with M0 = T * Mo * 1e7 (dyn·cm)
Mo2Mw = @(Mo) (2/3) * log10(T * Mo * 1e7) - 10.7;

ax1 = gca;
ax1.XAxisLocation = 'bottom';
ax1.Box  = 'off';
ax1.Layer = 'top';
drawnow

% Transparent axes on top for Mw
ax2 = axes('Position', ax1.Position, ...
           'Color', 'none', ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right');

ax2.YLim  = ax1.YLim;
ax2.YTick = [];  % no y-ticks on the top axis

% Match Mw range to Mo range via the transform
ax2.XLim = Mo2Mw(ax1.XLim);

xlabel(ax2, 'Equivalent 500-year Mw')
set(ax2,'FontSize',15,'XColor','k')


