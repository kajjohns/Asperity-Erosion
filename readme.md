# Eroding Asperity Inversion Workflow

This repository contains MATLAB codes to invert geodetic data for **eroding asperities** on a subduction megathrust. The codes accompany the submitted manuscript: Johnson, K.M., and E.M. Sherrill (2025, submitted), Non-stationary Locked Asperities and Subduction Zone Earthquake Potential. 

The workflow:

1. Builds a **triangular mesh** of the plate interface.
2. Computes **elastic half-space Green’s functions** for geodetic observables.
3. Runs a **Metropolis MCMC inversion** for asperity radii and process-zone stressing rates.
4. Produces **summary plots** and **fancy map-based visualizations** of coupling, locking probability, stressing rates, and data fit.
5. Computes **moment rates** and equivalent **500-year Mw** statistics from the posterior.

The codes are currently set up for **Cascadia**, but the workflow is generic and can be adapted to other subduction zones.

Also included in this repository are codes for reproducing and plotting the inversion results  for the three subduction zones in the manuscript. These scripts are in folders called Run_Nankai_Inversion, Run_Cascadia_Inversion, and Run_Hikurangi_Inversion. In each folder run the scripts in the following order: (1)setup_, (2) mcmc_inverions_, (3) plot_MCMC_inversion_.

Note: The only files not provided are the very large hmat binary files that are generated in the build_GFs step (or setup step).  These will need to be built by setting compute_hmat = true;

Note: You may need to compile and Mex hm_mvp codes (heirarchical matrix, matrix-vector product). Instructions for this are in the readme in the hmmvp0.16 folder. 
 
 
---

## Repository Structure

A typical layout (names as used in the scripts):

```text
.
├── build_mesh.m
├── build_GFs.m
├── mcmc_inversion.m
├── plot_mcmc_inversion_results.m
├── plot_fancy_mcmc_inversion_results.m
├── setup_Hmat.m              % H-matrix setup script (stress kernel)
├── make.m                    % builds hmmvp MEX files
├── eroding_asperity_simulation.m  % main driver (if used in your setup)
├── hmmvp0.16/                % H-matrix / fast kernel implementation
├── tools/
│   ├── llh2local.m
│   ├── local2llh.m
│   ├── make_tri_mesh.m
│   ├── make_triangular_patch_stuff.m
│   ├── make_dispG_triangular.m
│   ├── make_baseline_rate_changes.m
│   ├── Get_Gs_baselines.m
│   ├── get_locked_and_ring_tau.m
│   ├── get_locked_radii_indices.m
│   ├── get_ring_tau.m
│   ├── metropolis_log.m
│   ├── slanCM/               % colormap utilities
│   └── dem/                  % DEM/ETOPO1 plotting utilities
├── ne_10m_coastline/         % coastline shapefiles or mat files
├── Cascadia_obs/             % optional observation-specific helpers
├── data/
│   ├── Cascadia_contours.txt
│   ├── Cascadia_rake_rate.txt
│   ├── cascadia_horizontal_gps_*.txt
│   ├── Cascadia_vertical_homogenized.txt
│   ├── origin.txt
│   └── (any other region-specific data files)
├── hmat/                     % H-matrix storage (created by setup_Hmat)
└── cmap.mat                  % custom red–blue colormap
```

---

## Dependencies

- **MATLAB** (R2018a or later recommended).
- **Statistics & Machine Learning Toolbox** (for KDTreeSearcher).
- **hmmvp** (v0.16) H-matrix library (included in `hmmvp0.16/`).
- **Nikkhoo triangular dislocation code** (used via `tools/Nikkhoo/`).
- **ETOPO1 DEM** access or pre-downloaded tiles (used by `tools/dem/` and `plot_global_etopo1`).
- Coastline data in `ne_10m_coastline/`.

Make sure all subfolders are on your MATLAB path (although all scripts do `addpath` for what they need, but you can also add the repo root recursively).

---

## Data Formats

### 1. Interface Contours

`data/Cascadia_contours.txt`  
Columns:
- `lon` (deg)
- `lat` (deg)
- `depth` (km, **negative** down)

### 2. Convergence Rake & Rate

`data/Cascadia_rake_rate.txt`  
Columns:
- `lon` (deg)
- `lat` (deg)
- `rake` (deg; direction of overriding plate wrt subducting plate)
- `rate` (mm/yr; plate convergence speed)

### 3. Reference Origin

`data/origin.txt`  
Single row:
- `[lon0 lat0 elev0]` used by `llh2local` / `local2llh` to define the local Cartesian frame (km).

### 4. Horizontal GPS

`data/cascadia_horizontal_gps_*.txt`  
Columns:
- `lon` (deg)
- `lat` (deg)
- `Ve` (east velocity; mm/yr or m/yr)
- `Vn` (north velocity; mm/yr or m/yr)
- `Sige` (1σ east; same units as Ve)
- `Sign` (1σ north; same units as Vn)

Units are controlled with `is_hz_meter` in `build_GFs.m`.

### 5. Vertical GPS

`data/Cascadia_vertical_homogenized.txt`  
Columns:
- `lon` (deg)
- `lat` (deg)
- `Vu` (up velocity; mm/yr or m/yr)
- `Sigu` (1σ vertical; same units as Vu)

Units controlled with `is_vert_meter` in `build_GFs.m`.  
If no vertical data are used, set `vert_filename = ''`.

---

## Step 0 — Build H-matrix MEX Files (once per platform)

Before running any inversion, you must compile the hmmvp MEX files.

1. Edit `make.m` if needed (e.g., `isunix`, `use_omp`, compiler flags).
2. From MATLAB, in the repo root (or in `hmmvp0.16` if that’s where `make.m` lives):

   ```matlab
   make
   ```

This builds the required MEX binaries (e.g., for `hm_mvp`) for your OS/architecture.

---

## Step 1 — Build Triangular Mesh: `build_mesh.m`

**Purpose:** Construct a triangular mesh of the subduction interface and interpolate convergence rake/rate to patch centroids.

Key parameters inside `build_mesh.m`:

```matlab
contr_file     = 'Cascadia_contours.txt';
conv_rate_file = 'Cascadia_rake_rate.txt';
max_depth      = 46;             % km
intv           = 8;              % nominal patch size (km)
bbox           = [-128 38; -120 54];  % plotting bounds
```

Workflow:

1. Load depth contours and convert `(lon,lat,depth)` → local `(x,y,depth)` via `llh2local`.
2. Discard contours below `max_depth`.
3. Build a triangular mesh:

   ```matlab
   [el, nd] = make_tri_mesh(contrkm, intv);
   patch_stuff = make_triangular_patch_stuff(el, nd);
   ```

4. Interpolate rake and convergence rate to patch centroids.
5. Produce diagnostic plots:
   - Contours in local coordinates.
   - Triangular mesh.
   - Rake and rate fields.
   - Strike & dip vectors.

**Outputs:**

```matlab
save mesh_file el nd patch_stuff rakes rates
```

This `mesh_file.mat` is used by downstream scripts.

---

## Step 2 — Build Green’s Functions & H-matrix: `build_GFs.m`

**Purpose:**  
Load geodetic data, compute displacement Green’s functions (GFs) for slip/backslip on the mesh, and (optionally) construct the H-matrix stress kernel.

Key switches:

```matlab
invert_vel   = false;   % true: invert velocities; false: invert baselines
compute_hmat = false;   % true: (re)build H-matrix; false: use existing
compute_disp = true;    % true: compute GFs; false: load from disp_filename
disp_filename = 'test_disp';

horizontal_filename = 'cascadia_horizontal_gps_240412.txt';
is_hz_meter = true;

vert_filename = 'Cascadia_vertical_homogenized.txt';
is_vert_meter = false;

bbox = [-128 38; -120 54];
save_filename = 'setup_cascadia';
```

Workflow:

1. Add paths and load:
   - `mesh_file.mat` (from Step 1),
   - `origin.txt`,
   - GPS/vertical data.
2. Convert station locations to local `(x,y)` coordinates.
3. If `invert_vel == false`, construct **baseline observables** using:
   - `make_baseline_rate_changes`, which detects triangles that intersect the surface and builds baseline geometries.
4. Plot observed velocities and/or baseline elongation rates.
5. If `compute_disp` is true:
   - Compute GFs at GPS sites with `make_dispG_triangular`:
     ```matlab
     [Ge, Gn, Gu] = make_dispG_triangular(el, nd, [xysites z=0]', rakes, []);
     ```
   - Optionally compute baseline GFs with `Get_Gs_baselines`.
   - Adjust sign conventions for rake (GFs are flipped).
6. Run `setup_Hmat` (script) to build or load the H-matrix stress kernel (stored in `gamb.hmat.savefn` and typically written into `hmat/`).

**Output:**

- A `.mat` file with name `save_filename` containing:
  - mesh (`el`, `nd`, `patch_stuff`, `rakes`, `rates`),
  - data and GFs (`Ge`, `Gn`, `Gu`, `Gbase`, etc.),
  - baseline info if used,
  - H-matrix metadata (`gamb.hmat.*`).

---

## Step 3 — MCMC Inversion: `mcmc_inversion.m`

**Purpose:**  
Run a Metropolis MCMC for:

- **Asperity radii** (at fixed strike–depth node locations),
- **Ring stressing rates** around asperities (process-zone stress).

These parameters control the locked patches and associated backslip rates, which are forward modeled against the geodetic data using the H-matrix + GFs.

Key setup:

```matlab
num_Dpts = 30;
Zpts = [4 7 10 13 15 20 25 30];  % asperity node depths (km)

radii = 20 * ones(num_Dpts * length(Zpts), 1);
stepsize_radii = 10 * ones(size(radii));

minrad = 10;  % km minimum radius (smaller -> effectively off)

ring_taus = 0.15e5 * ones(size(radii));    % Pa/yr
stepsize_ring_tau = 0.15e5 * ones(size(ring_taus));

D = 8;  % km, ring width

continuing      = true;
startval        = false;
starting_folder = 'Cascadia_outputs';
folder_name     = 'Cascadia_outputs';
```

Workflow:

1. Add paths and ensure you have loaded the setup (`setup_cascadia.mat`).
2. Initialize H-matrix via:

   ```matlab
   [gamb.hmat.id, nnz] = hm_mvp('init', gamb.hmat.savefn, 4);
   ```

3. Rotate mesh so average strike aligns with y-axis; build a grid of asperity nodes in (strike, depth) space; map them to nearest mesh centroids.
4. Build KD-trees and precompute geometric factors (for efficient ring construction and patch queries).
5. Set up data vectors:
   - `G`, `GG`, `dd` (unweighted and normalized GFs and data),
   - incorporate baselines and/or verticals depending on switches.
6. If `continuing` or `startval` is `true`, load the last sample as the initial point.
7. Compute initial locked patches and creep/backslip distribution using:
   - `get_locked_and_ring_tau`,
   - `compute_creeprate` (H-matrix).
8. Initialize output files in `folder_name/`:
   - `M_radii.txt`, `M_ring_tau.txt`, `logrho.txt`, `dhat.txt`, `locked_index.txt`, `creep_rates.txt`.
9. Run Metropolis loop (default `for iter = 1:1e6`):
   - Propose a perturbation to a single parameter.
   - Recompute local creep/backslip using H-matrix (local updates with `idx_rbig`, `idx_rsmall`).
   - Compute predicted data, log-likelihood, and accept/reject.
   - Append parameters, log-likelihood, predictions, locked indices, and creep rates to text files every full pass through the parameter set.
   - Plot log-likelihood trace as a quick diagnostic.

You can reduce the iteration count for testing.

---

## Step 4 — Basic Posterior Maps: `plot_mcmc_inversion_results.m`

**Purpose:**  
Read MCMC outputs and generate basic maps of:

- Probability of locking,
- Mean process-zone stressing rate,
- Mean creep rate,
- Coupling ratio,
- Slip deficit rate,
- Data fit (baselines and verticals).

Key inputs inside the script:

```matlab
folder_name = 'Cascadia_outputs';
discard     = 65;   % burn-in to discard
```

Workflow:

1. Load `logrho.txt` and plot to visually choose an appropriate burn-in.
2. Load chain files:
   - `M_radii.txt`, `M_ring_tau.txt`, `dhat.txt`, `locked_index.txt`, `creep_rates.txt`.
3. Discard the first `discard` rows from all chains.
4. Reconstruct per-patch ring stress (`Ring_Taus`) for each sample using:
   - `get_locked_radii_indices`,
   - `get_ring_tau`.
5. Convert patch nodes to `(lon,lat)` via `local2llh`.
6. Plot:
   - `mean(locked_index)` → probability of locking.
   - `mean(Ring_Taus,2)` → mean process-zone stress.
   - `mean(creep_rates)` → mean creep rate.
   - `1 - mean(creep_rates)'/(srate*1000)` → coupling ratio.
   - `(srate*1000) - mean(creep_rates)'` → slip deficit rate.
7. Compare observed vs predicted:
   - Baseline elongation rates (if `invert_vel == false`).
   - Vertical rates (observed, predicted, residuals) if vertical data are present.

These plots provide quick diagnostics and science figures in a simple map projection.

---

## Step 5 — Fancy Topographic Maps: `plot_fancy_mcmc_inversion_results.m`

**Purpose:**  
Generate publication-quality visualizations using shaded-relief ETOPO1:

- Coupling ratio, probability of locking, and mean ring stress over topography.
- Observed/model baseline elongation and vertical rates on DEM.
- Misfit vs distance from trench.
- Moment-rate statistics and Mw histograms with dual axes.

Key pieces:

- Uses `plot_global_etopo1` + `dem` (from `tools/dem/`) to draw shaded-relief maps with a given light azimuth (`LA`).
- Uses the same MCMC outputs as `plot_mcmc_inversion_results.m`, but overlaid on topography.
- Computes trench distance using the deepest mesh nodes and `dsearchn`.
- Computes:
  - Total and locked-only moment rates (`Mo_total`, `Mo_locking`),
  - Equivalent Mw for a 500-yr interval (`Mw_500_total`, `Mw_500_locking`),
  - Variance reduction for baselines and verticals (overall and within 300 km of the trench).
- Plots histograms of moment rates with a **second x-axis** in equivalent Mw.

Edit `lat_range`/`lon_range` and color limits as needed for other regions.

---

## Typical “From-Scratch” Run

In MATLAB:

```matlab
% 0. Build hmmvp MEX files (only once per platform)
make

% 1. Build triangular mesh and interface geometry
build_mesh

% 2. Compute Green’s functions and H-matrix (or load if already done)
build_GFs    % check invert_vel, compute_hmat, compute_disp, filenames

% 3. Run MCMC inversion (this can be expensive)
mcmc_inversion   % check num iterations, folder_name, continuing/startval

% 4. Basic posterior maps
plot_mcmc_inversion_results

% 5. Fancy DEM-based plots
plot_fancy_mcmc_inversion_results
```

---

## Adapting to Other Subduction Zones

To apply the workflow elsewhere:

1. Replace input files in `data/`:
   - Interface contours,
   - Convergence rake/rate grid,
   - GPS horizontal/vertical data,
   - Origin.
2. Update:
   - `bbox`, `lat_range`, `lon_range` (plots),
   - `max_depth`, `intv` (mesh),
   - Station/vertical filenames and unit flags in `build_GFs.m`.
3. Adjust:
   - MCMC hyperparameters (`num_Dpts`, `Zpts`, priors, step sizes),
   - Burn-in length (`discard`),
   - Plot ranges and color scales.

---

## Citation / Acknowledgments

If you use these codes or figures produced by them in a publication, please:

- Cite the associated **eroding asperity** paper(s) and any code-specific reference the PI provides.
- Acknowledge third-party libraries where appropriate:
  - hmmvp H-matrix library,
  - Nikkhoo triangular dislocation code,
  - ETOPO1 global relief model,
  - Natural Earth coastlines (if used in `ne_10m_coastline`).


---
