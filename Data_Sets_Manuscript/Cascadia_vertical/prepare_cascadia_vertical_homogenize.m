
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cascadia vertical velocity compilation and GIA-corrected field
%
% This script:
%   (1) Downloads and filters GNSS vertical velocities from the MIDAS
%       solution (IGS14) provided by the Nevada Geodetic Laboratory.
%   (2) Loads vertical velocities from leveling-based studies in Cascadia:
%       - Newton et al. (2021) vertical land-motion rates (Washington coast)
%       - Burgette et al. (2009) interseismic uplift rates (Oregon)
%   (3) Loads and applies GIA / vertical-motion corrections from Lau et al.
%       (2020) to remove long-wavelength isostatic and loading signals.
%   (4) Aligns the reference frames of the various datasets to MIDAS,
%       filters outliers using a local median filter, homogenizes the
%       spatial sampling on a regular lon–lat grid, and applies the Lau
%       correction.
%   (5) Produces maps and profiles of the combined vertical velocities,
%       and writes out a combined data file:
%           Cascadia_vertical_combined.txt
%       with columns: [lon  lat  v_vertical(mm/yr)  sigma(mm/yr)].
%
% DATA SOURCES AND CITATIONS
% --------------------------
% GNSS (MIDAS velocities):
%   Blewitt, G., Kreemer, C., Hammond, W. C., & Gazeaux, J. (2016).
%   MIDAS robust trend estimator for accurate GPS station velocities
%   without step detection. J. Geophys. Res. Solid Earth, 121(3),
%   2054–2068. https://doi.org/10.1002/2015JB012552
%   (Data from Nevada Geodetic Laboratory, http://geodesy.unr.edu)
%
% Leveling-based vertical velocities:
%   Burgette, R. J., Weldon II, R. J., & Schmidt, D. A. (2009).
%   Interseismic uplift rates for western Oregon and along-strike
%   variation in locking on the Cascadia subduction zone. J. Geophys.
%   Res. Solid Earth, 114, B01408.
%
%   Newton, T. J., Weldon, R., Miller, I. M., Schmidt, D., Mauger, G.,
%   Morgan, H., & Grossman, E. (2021).
%   An assessment of vertical land movement to support coastal hazards
%   planning in Washington State. Water, 13(3), 281.
%
% GIA / large-scale vertical corrections:
%   Lau, N., Blewitt, G., & Becker, T. W. (2020).
%   Present-Day Crustal Vertical Velocity Field for the Contiguous
%   United States. J. Geophys. Res. Solid Earth, 125, e2020JB020066.
%
% NOTE:
%   - This script assumes that:
%       * MIDAS toolbox / functions are in the "MIDAS" directory.
%       * Newton and Burgette level data are stored in:
%             ./Newton_Supplemental/Newton_Data_S2.txt
%             Burgett_2009.txt   (file name uses "Burgett" but data are
%                                from Burgette et al., 2009.)
%       * Lau vertical field is stored in:
%             ./Lau/CONUS_gps_verticals.nc
%       * DEM and plotting helper functions (plot_global_etopo1,
%         dem, plot_state_borders_llh_box, make_cmap) are in ./dem
%
%   - Uncertainties are floored at 0.5 mm/yr following a similar approach
%     to that adopted by Zeng for horizontal GNSS velocities.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ./MIDAS
addpath ./dem
addpath ./tools

%% Download and filter MIDAS GNSS vertical velocities

% Minimum length of GNSS time series (years) to use
minT = 5;

% Get MIDAS velocities in the IGS14 frame
[sta,lat,lon,h,lab,t1,tm,dT,m,ngood,numsol,ve,vn,vu,se,sn,su,xe50,xn50,xu50,rho] = ...
    GetMIDASVelocities('IGS14',[]);

% Geographic bounds for Cascadia subset (approx.)
minlon = -127;
maxlon = -115;
maxlat = 50;
minlat = 40;

% Select stations within time and geographic constraints
ind = dT >= minT & lon > minlon & lon < maxlon & ...
      lat < maxlat & lat > minlat;

lat = lat(ind);
lon = lon(ind);
ve  = ve(ind);
vn  = vn(ind);
vu  = vu(ind);
se  = se(ind);
sn  = sn(ind);
su  = su(ind);

% Store MIDAS vertical velocities and uncertainties
% Columns: [lon lat vu su]
midas_vels = [lon lat vu su];

%% Load leveling data sets (Newton and Burgette)

% Newton et al. (2021) vertical land-motion data
%   Columns assumed: [lon lat v_vertical sigma ...]
newton  = load('./Newton_Supplemental/Newton_Data_S2.txt');

% Burgette et al. (2009) leveling verticals
%   File name uses "Burgett_2009.txt" but data from Burgette et al. (2009).
%   Columns assumed: [lat lon ... v_vertical sigma ...] — see your own file spec.
burgett = load("./Burgett_data/Burgett_2009.txt");

%% Load Lau et al. (2020) GIA / large-scale vertical correction field

% Use Lau corrections to remove long-wavelength GIA / loading from the data
fn = './Lau/CONUS_gps_verticals.nc';
ni = ncinfo(fn);

% Load all variables in the netCDF into a structure
for i = 1:length(ni.Variables)
    vn = ni.Variables(i).Name;
    lau.(vn) = ncread(fn, vn);
end

% Compute correction to apply at each grid node:
%   gps_vu_smooth = observed GNSS vertical field
%   net_vu_smooth = modeled GIA + loading contribution
%   lau_correction = GIA+loading component to subtract
lau_correction = lau.gps_vu_smooth(:) - lau.net_vu_smooth(:);

%% Combine and frame-align vertical velocity data sets

% Start with MIDAS velocities: columns [lon lat vu]
all_vels = midas_vels(:,1:3);

% NOTE: Krogstad velocities are included in Newton, so we do not add them
% separately:
%   all_vels = [all_vels; krogstad(:,2) krogstad(:,1) krogstad(:,3)];

%----------------------------------------------------------------------
% Bring Newton velocities into alignment with MIDAS
%   For each Newton site, find the nearest MIDAS site, compare velocities,
%   and compute a mean offset to tie Newton to MIDAS vertical reference.
%----------------------------------------------------------------------

for k = 1:size(newton,1)

    % Distance between this Newton point and all existing (MIDAS) sites
    dist = sqrt((all_vels(:,1) - newton(k,1)).^2 + ...
                (all_vels(:,2) - newton(k,2)).^2 );

    % Find nearest MIDAS point and map its vertical velocity
    [y,i] = min(dist);
    newton_map_to_GPS(k) = all_vels(i,3);

    dist_Newton(k) = min(dist);
end

% Average vertical offset between Newton and MIDAS
adjust_n = mean(newton(:,3) - newton_map_to_GPS');

% Remove Newton velocities that are essentially redundant with MIDAS
% (very small separation in lon–lat space)
ind = dist_Newton < 10^-5;
newton(ind,:) = [];

%----------------------------------------------------------------------
% Bring Burgette velocities into alignment with MIDAS
%   Same idea: nearest-neighbor matching to estimate mean offset.
%----------------------------------------------------------------------

for k = 1:size(burgett,1)

    % Note: burgett(:,2) = lon, burgett(:,1) = lat (based on use here)
    dist = sqrt((all_vels(:,1) - burgett(k,2)).^2 + ...
                (all_vels(:,2) - burgett(k,1)).^2 );
    [y,i] = min(dist);
    burgett_map_to_GPS(k) = all_vels(i,3);

end

% Average vertical offset between Burgette and MIDAS
adjust_b = mean(burgett(:,6) - burgett_map_to_GPS');

% Append Newton and Burgette velocities after applying offsets
% Newton: [lon lat v_vertical - offset]
all_vels = [all_vels; [newton(:,1:2) newton(:,3) - adjust_n]];

% Burgette: note the (lon,lat) ordering in the file
all_vels = [all_vels; burgett(:,2) burgett(:,1) burgett(:,6) - adjust_b];

% Collect uncertainties (MIDAS + Newton + Burgette)
all_sigs = [midas_vels(:,4); newton(:,4); burgett(:,7)];

% Place a floor on uncertainties of 0.5 mm/yr
% (arbitrary, but consistent with the approach of Zeng for horizontal vels)
all_sigs(all_sigs < 0.5) = 0.5;

%% Outlier detection / filtering

% At each observation point, compute the median vertical velocity of all
% data within a specified radius. Replace values that deviate by more
% than a threshold with the local median.

filter_radius = 0.25;   % degrees (approx. neighborhood radius)
threshold     = 2;      % mm/yr: outlier threshold

for k = 1:size(all_vels,1)

    dist = sqrt((all_vels(:,1) - all_vels(k,1)).^2 + ...
                (all_vels(:,2) - all_vels(k,2)).^2);
    med  = median(all_vels(dist < filter_radius,3));

    if abs(all_vels(k,3) - med) > threshold
        filt_vel(k) = med;   % replace outlier with local median
    else
        filt_vel(k) = all_vels(k,3);
    end
end
filt_vel = filt_vel';

%% Homogenize spatial sampling on a regular lon–lat grid

% Create a grid of cells. Within each cell, retain the median velocity
% (and its location) and discard the rest of the points in that cell.

cell_dimension = 0.20;  % degrees (grid cell size)
lonmin = -128; lonmax = -114;
latmin = 39;  latmax = 50;

num_lon = round((lonmax - lonmin) / cell_dimension);
num_lat = round((latmax - latmin) / cell_dimension);

lon_edges = linspace(lonmin, lonmax, num_lon);
lat_edges = linspace(latmin, latmax, num_lat);

vels_homog   = [];
all_vels_filt = [all_vels(:,1:2) filt_vel];

for k = 1:num_lat-1
    for j = 1:num_lon-1

        ind = all_vels(:,1) > lon_edges(j)   & ...
              all_vels(:,1) < lon_edges(j+1) & ...
              all_vels(:,2) > lat_edges(k)   & ...
              all_vels(:,2) < lat_edges(k+1);

        if sum(ind) > 0
            % Keep the median velocity (and location) in this cell
            vels_homog(end+1,:) = median(all_vels_filt(ind,:),1);
        end

    end
end

%% Interpolate and subtract GIA corrections from Lau et al. (2020)

% Only use grid nodes where the correction is defined (non-NaN)
ind = ~isnan(lau_correction);

lon = lau.lon(:); lon = lon(ind);
lat = lau.lat(:); lat = lat(ind);

% Interpolate the Lau correction to the full set of observation locations
lau_correction_interp = griddata(lon, lat, lau_correction(ind), ...
                                 all_vels(:,1), all_vels(:,2), 'nearest');

% Interpolate the Lau correction to the homogenized grid
lau_correction_interp_homog = griddata(lon, lat, lau_correction(ind), ...
                                       vels_homog(:,1), vels_homog(:,2), ...
                                       'nearest');

% Apply the corrections: data minus GIA / loading correction
vert_vels       = filt_vel - lau_correction_interp;
vert_vels_homog = vels_homog(:,3) - lau_correction_interp_homog;

%% Plot maps: corrected vertical velocities (all data and homogenized)

% Common plotting region
lat_range = [38 50];
lon_range = [-128 -115];

%-----------------------------
% All data minus Lau correction
%-----------------------------
[A,lats,lons] = plot_global_etopo1(lat_range, lon_range);

figure
[h,I,z] = dem(lons, lats, A, 'LatLon', 'Legend', ...
              'LandColor', .8*ones(359,3), ...
              'SeaColor',  .6*ones(359,3), ...
              'Zlim',      [-2000 2000], ...
              'Contrast',  1);
hold on

% Add state borders
bbox_ll = [-128 -115 38 50];
plot_state_borders_llh_box('k', bbox_ll)
xlim(lon_range)
ylim(lat_range)

% Plot corrected vertical velocities
scatter(all_vels(:,1), all_vels(:,2), 30, vert_vels, 'fill')
caxis([-4 4])
colorbar
make_cmap
colormap(cmap)
title('Data minus Lau correction')

%----------------------------------
% Homogenized data minus Lau correction
%----------------------------------
[A,lats,lons] = plot_global_etopo1(lat_range, lon_range);

figure
[h,I,z] = dem(lons, lats, A, 'LatLon', 'Legend', ...
              'LandColor', .8*ones(359,3), ...
              'SeaColor',  .6*ones(359,3), ...
              'Zlim',      [-2000 2000], ...
              'Contrast',  1);
hold on

plot_state_borders_llh_box('k', bbox_ll)
xlim(lon_range)
ylim(lat_range)

scatter(vels_homog(:,1), vels_homog(:,2), 30, vert_vels_homog, 'fill')
caxis([-4 4])
colorbar
make_cmap
colormap(cmap)
title('Homogenized data minus Lau correction')

%% Plot map of vertical velocities before Lau correction (for comparison)

[A,lats,lons] = plot_global_etopo1(lat_range, lon_range);

figure
[h,I,z] = dem(lons, lats, A, 'LatLon', 'Legend', ...
              'LandColor', .8*ones(359,3), ...
              'SeaColor',  .6*ones(359,3), ...
              'Zlim',      [-2000 2000], ...
              'Contrast',  1);
hold on

plot_state_borders_llh_box('k', bbox_ll)
xlim(lon_range)
ylim(lat_range)

scatter(all_vels(:,1), all_vels(:,2), 30, filt_vel, 'fill')
caxis([-4 4])
colorbar
make_cmap
colormap(cmap)
title('Vertical data before applying Lau correction')

%% Plot along-strike profiles at selected latitudes

lats_profile = [41 43 45 47];

figure
for k = 1:length(lats_profile)

    subplot(2,2,k)

    % Select points within ±0.15° of the target latitude
    ind = abs(all_vels(:,2) - lats_profile(k)) < 0.15;

    % Plot Lau-corrected vertical velocities
    plot(all_vels(ind,1), vert_vels(ind), 'o')
    xlabel('Longitude')
    ylabel('Vertical velocity (mm/yr)')
    title(['Latitude = ' num2str(lats_profile(k))])

    ylim([-5 5])
    grid on

end

%% Write combined, corrected data file

% Columns: [lon, lat, vertical_velocity(Lau-corrected), sigma]
data = [all_vels(:,1:2) vert_vels all_sigs];

dlmwrite('Cascadia_vertical_combined.txt', data)
