
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mcmc_inversion.m
%
% Metropolis MCMC inversion for:
%   (1) Asperity geometry on a subduction interface, parameterized by
%       circular asperities with radii at fixed strike–depth node locations.
%   (2) Process-zone stressing rates (ring_taus) around each asperity.
%
% The model uses:
%   - A triangular mesh of the plate interface (el, nd, patch_stuff)
%   - Plate interface convergence rates "rates" (mm/yr)
%   - Elastic Green’s functions (G) relating backslip to geodetic data
%   - An H-matrix-based stress kernel (hmmvp) to compute creep rates
%
% At each MCMC step:
%   1) Propose a small perturbation to either:
%        - an asperity radius, or
%        - a ring stressing rate (ring_tau).
%   2) Compute which patches are locked / ring-loaded and the associated
%      backslip rate using H-matrix stress computations.
%   3) Forward model geodetic predictions (dhat) via Green’s functions.
%   4) Compute the log-likelihood relative to observed data and accept/
%      reject the proposal via Metropolis.
%   5) Save the chain of model parameters, likelihood, and predictions.
%
% REQUIRED PRE-COMPUTED INPUTS (from other scripts/workflow):
%   - build_mesh.m:
%       el, nd, patch_stuff, rakes, rates
%   - build_GFs.m and setup_Hmat.m:
%       Ge, Gn, Gu, Gbase, G, GG, dd, Vbase, Sigbase, Ve, Vn, Vu,
%       Sige, Sign, Sigu, invert_vel, vert_filename
%       gamb.hmat.savefn (file with H-matrix stored)
%
% KEY VARIABLES / NOTATION:
%   num_Dpts         : number of asperity nodes along strike
%   Zpts             : depths of asperity nodes (km)
%   radii            : radius of each asperity (km), one per node in ptsDZ
%   ring_taus        : process-zone stressing rate (Pa/yr) around asperities
%   D                : ring width (km)
%   minrad           : minimum allowable radius; radii below this set to zero
%   srate            : imposed background stressing / creep rate
%   bs_rate          : backslip rate (mm/yr) on each patch from asperities
%   creep_rate       : residual creep rate (data-derived) = rates - bs_rate
%   G, GG, dd        : Green’s matrices and normalized data
%   i_locked         : indices of locked (or asperity) patches
%
% OUTPUT (written incrementally to text files in folder_name/):
%   - M_radii.txt       : sampled asperity radii
%   - M_ring_tau.txt    : sampled ring stressing rates
%   - logrho.txt        : log-likelihood for each stored sample
%   - dhat.txt          : predicted data for each stored sample
%   - locked_index.txt  : locked patch flags / indices for each sample
%   - creep_rates.txt   : creep rates for each sample
%
% NOTES:
%   - This script is long-running (1e6 iterations by default).
%   - Use "continuing" and "startval" flags to resume a previous run or
%     start from a different folder’s final model.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  BEGIN SETUP (MCMC hyperparameters and asperity parameterization)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of asperity nodes along strike
num_Dpts = 30;

% Depths (km) of asperity nodes (Z direction)
Zpts = [4 7 10 13 15 20 25 30];

% Starting values for asperity radii (km)
radii = 20*ones(num_Dpts*length(Zpts),1);

% Proposal step size for radii (km)
stepsize_radii = 10*ones(size(radii));

% Minimum allowable radius (km) - smaller radii are effectively turned off
minrad = 10;

% Starting values for stressing rate on rings around asperities (Pa/yr)
ring_taus = 0.15e5*ones(size(radii));

% Proposal step size for ring stressing rates (Pa/yr)
stepsize_ring_tau = 0.15e5*ones(size(ring_taus));

% Width of process zone ring around asperity (km)
D = 8;

% Control chain start / continuation:
%   continuing = true  -> continue a previous chain in folder_name
%   startval   = true  -> initialize from the end of starting_folder
continuing       = true;            % continue sampling or start new
startval         = false;           % if true, use starting_folder values
starting_folder  = 'Cascadia_outputs';
folder_name      = 'Cascadia_outputs';

%%%%  END SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------
%  SOFTWARE CONFIGURATION AND OUTPUT FOLDER
% -------------------------------------------------------------------------
addpath hmmvp0.16
addpath tools
addpath ./Cascadia_obs

% Make output folder if it does not exist
if ~exist(folder_name, 'dir')
    mkdir(folder_name)
end

% Initialize H-matrix operator (hm_mvp) using precomputed file
%   gamb.hmat.savefn must exist and point to a stored H-matrix
[gamb.hmat.id, nnz] = hm_mvp('init', gamb.hmat.savefn, 4);

%% ------------------------------------------------------------------------
%  MESH GEOMETRY AND COORDINATE ROTATION
% -------------------------------------------------------------------------
% Recompute patch geometry from mesh (if not already in workspace)
patch_stuff = make_triangular_patch_stuff(el, nd);

centroids = patch_stuff.centroids_faces;   % patch centroids (x,y,z)
strikes   = patch_stuff.strike_faces;      % patch strike angles (degrees)

% Rotate and shift interface so that average strike aligns with the y-axis.
%   This simplifies defining asperity nodes along "strike" (Dpts).
str = mean(strikes)*pi/180;                % average strike in radians
R   = [cos(str) -sin(str); ...
       sin(str)  cos(str)];

centroids_rot      = (R * centroids(:,1:2)')';  % rotate x-y
centroids_rot(:,3) = centroids(:,3);            % keep depth as-is

% Shift along rotated strike direction so minimum y is at zero
centroids_rot(:,2) = centroids_rot(:,2) - min(centroids_rot(:,2));

%% ------------------------------------------------------------------------
%  DEFINE ASPERITY NODE LOCATIONS (ALONG STRIKE / DEPTH)
% -------------------------------------------------------------------------
% Nodes along rotated strike direction
Dpts = linspace(0, max(centroids_rot(:,2)), num_Dpts);

% Build grid of asperity nodes (D along strike, Z in depth)
[ptsD, ptsZ] = meshgrid(Dpts, Zpts);
ptsDZ = [ptsD(:) ptsZ(:)];   % [num_nodes x 2] strike-depth nodes

% Find nearest mesh centroid for each asperity node
%   A is in (strike, depth) coordinates: depth is negative z
A    = [centroids_rot(:,2), -centroids_rot(:,3)];
Mdlc = KDTreeSearcher(A);          % kd-tree for nearest neighbors
[idx, dist] = knnsearch(Mdlc, ptsDZ);
asperity_centroids = centroids(idx,:);  % physical centroids tied to nodes

%% ------------------------------------------------------------------------
%  STRESSING RATE SETUP AND DATA UNCERTAINTIES
% -------------------------------------------------------------------------
% Plate interface convergence rates "rates" are in mm/yr; convert to m/yr
srate = rates/1000;

nel = size(el,1);   % number of mesh elements

% Enforce minimum uncertainty for baseline data
% (avoid over-weighting baselines with unrealistically small errors)
Sigbase(Sigbase < 1) = 1;

%% ------------------------------------------------------------------------
%  BUILD GREEN’S FUNCTION MATRICES FOR CHOSEN DATA TYPE
% -------------------------------------------------------------------------
% G   : maps slip / backslip to data (unweighted)
% GG  : maps slip / backslip to normalized data (G ./ sigma)
% dd  : observed data normalized by data uncertainties
if invert_vel
    % Invert horizontal velocities (east + north)
    G  = [Ge; Gn];
    GG = [Ge./repmat(Sige,1,size(Ge,2)); ...
          Gn./repmat(Sign,1,size(Gn,2))];
    dd = [Ve./Sige; Vn./Sign];
else
    % Invert baseline elongation rates
    G  = Gbase;
    GG = Gbase./repmat(Sigbase,1,size(Gbase,2));
    dd = Vbase./Sigbase;
end

% Append vertical GPS data if available
if ~isempty(vert_filename)
    G  = [G;  Gu];
    GG = [GG; Gu./repmat(Sigu,1,size(Gu,2))];
    dd = [dd; Vu./Sigu];
end

num = 1:nel;   % element index vector (for convenience)

%% ------------------------------------------------------------------------
%  INITIAL VALUES: CONTINUE PREVIOUS CHAIN OR START FROM OTHER FOLDER
% -------------------------------------------------------------------------
if continuing
    % Load last state from current folder
    data      = load(fullfile(folder_name, 'M_radii.txt'));
    radii     = data(end,:)';
    data      = load(fullfile(folder_name, 'M_ring_tau.txt'));
    ring_taus = data(end,:)';
end

if startval
    % Load last state from a different run (starting_folder)
    data      = load(fullfile(starting_folder, 'M_radii.txt'));
    radii     = data(end,:)';
    data      = load(fullfile(starting_folder, 'M_ring_tau.txt'));
    ring_taus = data(end,:)';
end

% Ensure radii are non-negative
radii(radii < 0) = 0;

% Vector of unknown parameters in stacked form:
%   X = [radii; ring_taus]
X        = [radii; ring_taus];
stepsize = [stepsize_radii; stepsize_ring_tau];

%% ------------------------------------------------------------------------
%  PRECOMPUTATIONS FOR GEOMETRY / RING CONSTRUCTION
% -------------------------------------------------------------------------
% Angles to parameterize asperity ring boundary
thetas   = linspace(0, 2*pi, 50);
cos_t    = cos(thetas(:));
sin_t    = sin(thetas(:));

% Precompute tangent of dip for each face
tan_dips = tan(patch_stuff.dip_faces*pi/180);

% KD-tree in rotated strike-depth space for ring identification
KDTree = KDTreeSearcher([centroids_rot(:,2), -centroids_rot(:,3)]);

% Find all patches within two radii of each asperity centroid
Mdl = KDTreeSearcher(centroids);   % kd-tree in original 3D centroid space
rbig   = 250;    % large radius of influence (km)
rsmall = 150;    % smaller radius of influence (km)

[idx_rbig,   dists] = rangesearch(Mdl, asperity_centroids, rbig);
[idx_rsmall, dists] = rangesearch(Mdl, asperity_centroids, rsmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Forward model (initial calculation for starting state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute locked patches and ring stressing for current radii/ring_taus
[i_locked, ring_tau] = get_locked_and_ring_tau(...
    ptsDZ, radii, ring_taus, D, centroids_rot, minrad, ...
    cos_t, sin_t, tan_dips, KDTree);

% Compute backslip rate on patches (mm/yr) using H-matrix solver
%   (global update: all elements at once)
bs_rate = compute_creeprate( ...
    gamb.hmat.id, 1:nel, 1:nel, nel, i_locked, ring_tau, srate, srate);

% Creep rate (mm/yr) on interface is the remainder after backslip
creep_rate = rates - bs_rate;

% Forward model predictions for the geodetic data
dhat  = G * bs_rate;
ddhat = GG * bs_rate;

% Log-likelihood based on normalized data misfit
logrho = -0.5 * (dd - ddhat)' * (dd - ddhat);

% Store current state as "previous"
Xprev            = X;
logrhoprev       = logrho;
dhatprev         = dhat;
i_locked_prev    = i_locked;
radii_prev       = radii;
ring_tau_prev    = ring_tau;
creep_rate_prev  = creep_rate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% File setup for saving the MCMC chain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If not resuming, create new (empty) output files
if ~continuing
    fid = fopen(fullfile(folder_name,'M_radii.txt'),'w');         fclose(fid);
    fid = fopen(fullfile(folder_name,'M_ring_tau.txt'),'w');      fclose(fid);
    fid = fopen(fullfile(folder_name,'logrho.txt'),'w');          fclose(fid);
    fid = fopen(fullfile(folder_name,'dhat.txt'),'w');            fclose(fid);
    fid = fopen(fullfile(folder_name,'locked_index.txt'),'w');    fclose(fid);
    fid = fopen(fullfile(folder_name,'creep_rates.txt'),'w');     fclose(fid);
end

% Open files in append mode for writing chain samples
fidM_radii   = fopen(fullfile(folder_name,'M_radii.txt'),'a');
fidM_tau     = fopen(fullfile(folder_name,'M_ring_tau.txt'),'a');
fidrho       = fopen(fullfile(folder_name,'logrho.txt'),'a');
fiddhat      = fopen(fullfile(folder_name,'dhat.txt'),'a');
fid_locked   = fopen(fullfile(folder_name,'locked_index.txt'),'a');
fid_creeping = fopen(fullfile(folder_name,'creep_rates.txt'),'a');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Metropolis Sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numparams = length(X);   % total number of parameters (radii + ring_taus)
numaccept = 0;           % counter for accepted proposals

rng('shuffle');          % initialize random number generator

Logrhos = zeros(1,1e6);  % store log-likelihood history

for iter = 1:1e6

    %-------------------------------------------------------
    % Propose parameter update (random single-parameter step)
    %-------------------------------------------------------
    count = mod(iter, numparams) + 1;   % parameter index to update

    % Random step direction and magnitude
    r = (-1)^round(rand(1)) * rand(1);  % +/- random
    r = r * stepsize(count);            % scale by parameter-specific step
    X(count) = X(count) + r;            % apply proposal

    % Split parameter vector into radii and ring_taus
    radii     = X(1:length(radii));
    ring_taus = X(1+length(radii):end);

    % Keep track of which asperity node is being modified
    % (used to determine local region for efficient creep updates)
    if mod(iter,numparams) == 1 || mod(iter,numparams/2) == 1
        % reset asperity index after cycling through half/full parameter set
        asp_num = 1;
    else
        asp_num = asp_num + 1;
    end

    %-------------------------------------------------------
    % Check parameter bounds (e.g., ring_taus >= 0 & <= max)
    %-------------------------------------------------------
    if all(ring_taus >= 0 & ring_taus <= 7e5)

        %---------------------------------------------------
        % Forward model calculation for proposed state
        %---------------------------------------------------
        [i_locked, ring_tau] = get_locked_and_ring_tau(...
            ptsDZ, radii, ring_taus, D, centroids_rot, minrad, ...
            cos_t, sin_t, tan_dips, KDTree);

        % Local update of backslip rate using H-matrix:
        %   idx_rbig{asp_num} / idx_rsmall{asp_num} specify elements
        %   within large/small distance of current asperity
        bs_rate = compute_creeprate( ...
            gamb.hmat.id, idx_rbig{asp_num}, idx_rsmall{asp_num}, ...
            nel, i_locked, ring_tau, srate, bs_rate);

        % Updated creep rate (mm/yr) along interface
        creep_rate = rates - bs_rate;

        % Predictions for geodetic data (unweighted and normalized)
        dhat  = G * bs_rate;
        ddhat = GG * bs_rate;

        % New log-likelihood for proposed state
        logrho2 = -0.5 * (dd - ddhat)' * (dd - ddhat);

        %---------------------------------------------------
        % Metropolis acceptance step
        %---------------------------------------------------
        accept = metropolis_log(1, logrho, 1, logrho2);

    else
        % Proposal violates bounds: automatically reject
        accept = 0;
    end

    %-------------------------------------------------------
    % Update chain state based on acceptance decision
    %-------------------------------------------------------
    if accept == 1
        % Accept proposal
        logrho          = logrho2;
        Xprev           = X;
        dhatprev        = dhat;
        logrhoprev      = logrho;
        i_locked_prev   = i_locked;
        radii_prev      = radii;
        ring_tau_prev   = ring_taus;
        creep_rate_prev = creep_rate;
        bs_rate_prev    = bs_rate;
        numaccept       = numaccept + 1;
    else
        % Reject proposal: revert to previous state
        X          = Xprev;
        dhat       = dhatprev;
        logrho     = logrhoprev;
        i_locked   = i_locked_prev;
        radii      = radii_prev;
        ring_taus  = ring_tau_prev;
        creep_rate = creep_rate_prev;
        bs_rate    = bs_rate_prev;
    end

    %-------------------------------------------------------
    % File writing: record a sample once per parameter cycle
    %-------------------------------------------------------
    if count == 1
        % Append current model and likelihood to chain files
        fprintf(fidM_radii,'\n');   fprintf(fidM_radii,'%6.8f\t',radii');
        fprintf(fidrho,'\n');       fprintf(fidrho,'%6.8f\t',logrho);
        fprintf(fiddhat,'\n');      fprintf(fiddhat,'%6.8f\t',dhat');
        fprintf(fidM_tau,'\n');     fprintf(fidM_tau,'%6.8f\t',ring_taus');
        fprintf(fid_locked,'\n');   fprintf(fid_locked,'%6.8f\t',i_locked');
        fprintf(fid_creeping,'\n'); fprintf(fid_creeping,'%6.8f\t',creep_rate');

        % Simple progress reporting to command window
        disp(['Completed sample number ' num2str(iter)])
        disp(['Acceptance rate: ' num2str(numaccept/iter*100) '%'])
        disp(['Log likelihood: ' num2str(logrho)])
    end

    %-------------------------------------------------------
    % Store likelihood history and diagnostic plot
    %-------------------------------------------------------
    Logrhos(iter) = logrho;

    % Quick diagnostic plot of log-likelihood every 100 iterations
    if mod(iter,100) == 0
        figure(100); clf
        plot(Logrhos(1:iter))
        xlabel('Iteration')
        ylabel('Log-likelihood')
        title('MCMC log-likelihood trace')
        drawnow
    end
end

%----------------------------------------------------------------------
% Close all open chain files at the end of the run
%----------------------------------------------------------------------
fclose(fidM_radii);
fclose(fidM_tau);
fclose(fidrho);
fclose(fiddhat);
fclose(fid_locked);
fclose(fid_creeping);



