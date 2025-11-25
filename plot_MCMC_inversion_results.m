
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_mcmc_inversion_results.m
%
% Post-processing and visualization of MCMC inversion results for
% asperity geometry and process-zone stressing rates.
%
% This script:
%   1) Loads MCMC chain outputs from folder_name:
%        - logrho.txt        : log-likelihood per stored sample
%        - M_radii.txt       : asperity radii samples
%        - M_ring_tau.txt    : ring stressing rate samples
%        - dhat.txt          : predicted data for each sample
%        - locked_index.txt  : locked-patch indicator per sample
%        - creep_rates.txt   : creep rates per sample
%
%   2) Plots log-likelihood vs sample index to help identify burn-in.
%
%   3) Discards a specified number of burn-in samples.
%
%   4) Reconstructs process-zone stressing rates on patches (Ring_Taus)
%      for each retained sample using:
%        - get_locked_radii_indices
%        - get_ring_tau
%
%   5) Computes posterior means and plots:
%        - Probability of locking
%        - Mean process-zone stressing rate
%        - Mean creep rate
%        - Coupling ratio
%        - Slip deficit rate
%
%   6) Compares observed and predicted data:
%        - Baseline elongation rates (if invert_vel == false)
%        - Vertical velocities: observed, predicted, residuals (if vertical
%          data are present via vert_filename)
%
% REQUIRED VARIABLES / INPUTS (from prior scripts/workspace):
%   - folder_name     : name of output folder used by mcmc_inversion.m
%   - ptsDZ           : asperity node coordinates (strike-depth grid)
%   - centroids_rot   : rotated centroids of triangular patches
%   - patch_stuff     : struct with dip_faces and mesh geometry
%   - minrad, D       : minimum radius, ring width (km)
%   - el, nd          : mesh connectivity and nodes
%   - origin          : reference for local <-> geographic conversion
%   - srate           : background plate rate (m/yr)
%   - bbox            : [lon_min lat_min; lon_max lat_max] for plots
%   - invert_vel      : true if velocities inverted, false if baselines
%   - baselines, Vbase, L, xysites, xysites_vert, Vu, vert_filename
%
% DEPENDENCIES:
%   - tools/slanCM/slanCM.m : colormap generation
%   - cmap.mat              : custom red-blue colormap
%   - get_locked_radii_indices.m
%   - get_ring_tau.m
%   - local2llh.m, plot_coast.m, plot_coast_xy.m, cline.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Folder containing MCMC output text files
folder_name = 'Cascadia_outputs';

% Load log-likelihood history and plot to identify burn-in period
eval(['load ./' folder_name '/logrho.txt'])
figure
plot(logrho)
xlabel('Sample index')
ylabel('Log-likelihood')
title('Log-likelihood trace (identify burn-in)')

% Discard burn-in samples:
%   The first "discard" samples will be removed from each chain file.
discard = 65;

%% End input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add path for custom colormap tools
addpath ./tools/slanCM/

% Load MCMC chain outputs from folder
eval(['load ./' folder_name '/M_radii.txt'])
eval(['load ./' folder_name '/M_ring_tau.txt'])
eval(['load ./' folder_name '/dhat.txt'])
eval(['load ./' folder_name '/locked_index.txt'])
eval(['load ./' folder_name '/creep_rates.txt'])

% Toss out burn-in samples from all chains
M_radii(1:discard,:)       = [];
M_ring_tau(1:discard,:)    = [];
dhat(1:discard,:)          = [];
locked_index(1:discard,:)  = [];
creep_rates(1:discard,:)   = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct process-zone stressing rate (Ring_Taus) for each sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear Ring_Taus
for k = 1:size(M_radii,1)

    % Extract radii and ring_taus for sample k
    radii     = M_radii(k,:)';
    ring_taus = M_ring_tau(k,:)';

    % Determine locked patches for these radii
    i_locked = get_locked_radii_indices( ...
        ptsDZ, radii, centroids_rot, patch_stuff.dip_faces, minrad);

    % Compute accumulated process-zone stress (per patch) for this sample
    Ring_Taus(:,k) = get_ring_tau( ...
        ring_taus, ptsDZ, D, radii, centroids_rot, ...
        i_locked, patch_stuff.dip_faces, minrad);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup colormaps and coordinate transforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct custom colormap for coupling (blend of Blues and Inferno)
c1 = colormap(flipud(slanCM('inferno')));
c2 = colormap(flipud(slanCM('Blues')));
mycmap = [c2(1:2:end,:); c1(1:2:end,:)];

% Convert mesh nodes from local (x,y) to (lon,lat)
nd_ll = local2llh(nd(:,1:2)', origin)';   % [lon, lat, (optional z)]

% Load custom red-blue colormap for other plots
load cmap



addpath ./ne_10m_coastline/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot probability of locking (posterior mean of locked_index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
h = trisurf(el, nd_ll(:,1), nd_ll(:,2), nd(:,3), ...
            mean(locked_index), 'edgecolor','none');
colormap(flipud(pink))
colorbar
set(gca,'fontsize',15)
title('Probability of Locking')
view(2)     % plan view
hold on
plot_coast(bbox)
xlim([-129 -120])
ylim([40 50])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot process zone stressing rate (posterior mean of Ring_Taus)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
h = trisurf(el, nd_ll(:,1), nd_ll(:,2), nd(:,3), ...
            mean(Ring_Taus,2), 'edgecolor','none');
colormap(flipud(pink))
colorbar
set(gca,'fontsize',15)
title('Mean Negative Stressing Rate')
view(2)
hold on
plot_coast(bbox)
xlim([-129 -120])
ylim([39 51])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot creep rate (posterior mean)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
h = trisurf(el, nd_ll(:,1), nd_ll(:,2), nd(:,3), ...
            mean(creep_rates), 'edgecolor','none');
colormap(turbo)
colorbar
set(gca,'fontsize',15)
title('Mean Creep Rate')
view(2)
hold on
plot_coast(bbox)
xlim([-129 -120])
ylim([39 51])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot coupling ratio: 1 - (mean creep rate / plate rate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
h = trisurf(el, nd_ll(:,1), nd_ll(:,2), nd(:,3), ...
            1 - mean(creep_rates,1)' ./ (srate*1000), 'edgecolor','none');
colormap(cmap)
colorbar
set(gca,'fontsize',15)
title('Coupling Ratio')
view(2)
hold on
plot_coast(bbox)
xlim([-129 -120])
ylim([40 50])
caxis([-1 1])       % symmetric color scale
colormap(mycmap)    % use custom coupling colormap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot slip deficit rate: (plate rate - mean creep rate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
h = trisurf(el, nd_ll(:,1), nd_ll(:,2), nd(:,3), ...
            (srate*1000) - mean(creep_rates)', 'edgecolor','none');
colormap(cmap)
colorbar
set(gca,'fontsize',15)
title('Slip Deficit Rate')
view(2)
hold on
plot_coast(bbox)
xlim([-129 -120])
ylim([39 51])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot fit to data (baselines and verticals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Baseline elongation rates (if inverting baselines)
if ~invert_vel  % baselines case

    % Number of baselines
    Nbase = size(baselines,1);

    % Split dhat into baseline and vertical contributions
    dhat_base = dhat(:,1:Nbase);
    dhat_u    = dhat(:,1+Nbase:end);

    % Convert site and endpoint coordinates to lat/lon
    sites_ll  = local2llh(xysites', origin)';
    vert_ll   = local2llh(xysites_vert', origin)';

    BaseEnds_ll = [ ...
        local2llh(BaseEnds(:,1:2)', origin)' ...
        local2llh(BaseEnds(:,3:4)', origin)'];

    %---------------------------------------------------
    % Observed baseline elongation rates
    %---------------------------------------------------
    figure;
    subplot(1,2,1)
    hold on;
    axis equal
    for k = 1:size(BaseEnds,1)
        cline([BaseEnds_ll(k,1) BaseEnds_ll(k,3)], ...
              [BaseEnds_ll(k,2) BaseEnds_ll(k,4)], ...
              [Vbase(k)/L(k) Vbase(k)/L(k)]);
    end
    title('Observed baseline elongation rates')
    colormap(cmap)
    colorbar
    caxis([-0.1 0.1])
    plot_coast(bbox)
    xlim([-129 -120])
    ylim([40 50])

    %---------------------------------------------------
    % Predicted baseline elongation rates (posterior mean)
    %---------------------------------------------------
    subplot(1,2,2)
    hold on;
    axis equal
    v = mean(dhat,1);     % mean predicted data over samples
    for k = 1:size(BaseEnds,1)
        cline([BaseEnds_ll(k,1) BaseEnds_ll(k,3)], ...
              [BaseEnds_ll(k,2) BaseEnds_ll(k,4)], ...
              [v(k)/L(k) v(k)/L(k)]);
    end
    title('Predicted baseline elongation rates')
    colormap(cmap)
    colorbar
    caxis([-0.1 0.1])
    plot_coast(bbox)
    xlim([-127 -117])
    ylim([40 50])

else
    % Horizontal-velocity inversion case:
    % (Add plotting of observed vs predicted horizontal velocities here if desired.)
end

%% Vertical velocities (if vertical data exist)
if ~isempty(vert_filename)

    %------------------------------------------------------------------
    % Observed vertical velocities
    %------------------------------------------------------------------
    figure;
    subplot(1,3,1)
    hold on
    scatter(xysites_vert(:,1), xysites_vert(:,2), 50, Vu, 'fill')
    axis equal
    title('Observed Vertical')
    plot_coast_xy(bbox, origin, 'k')
    colormap(cmap)
    colorbar
    grid on
    set(gca,'fontsize',15)
    xlim([-200 600])
    ylim([-600 600])
    caxis([-5 5])

    %------------------------------------------------------------------
    % Predicted vertical velocities (posterior mean)
    %------------------------------------------------------------------
    vu = mean(dhat(:, length(Vbase)+1:end));   % mean predicted verticals

    subplot(1,3,2)
    hold on
    scatter(xysites_vert(:,1), xysites_vert(:,2), 50, vu, 'fill')
    axis equal
    title('Predicted Vertical')
    plot_coast_xy(bbox, origin, 'k')
    colormap(cmap)
    colorbar
    grid on
    set(gca,'fontsize',15)
    xlim([-200 600])
    ylim([-600 600])
    caxis([-5 5])

    %------------------------------------------------------------------
    % Vertical residuals: observed - predicted
    %------------------------------------------------------------------
    subplot(1,3,3)
    hold on
    scatter(xysites_vert(:,1), xysites_vert(:,2), 50, Vu - vu', 'fill')
    axis equal
    title('Vertical Residual')
    plot_coast_xy(bbox, origin, 'k')
    colormap(cmap)
    colorbar
    grid on
    set(gca,'fontsize',15)
    xlim([-200 600])
    ylim([-600 600])
    caxis([-5 5])

end



