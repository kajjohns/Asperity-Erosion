# Cascadia Vertical Velocity Compilation and GIA-Corrected Field

This repository contains MATLAB scripts and supporting data used to
compile, homogenize, and correct vertical land-motion measurements
across the Cascadia subduction zone. The workflow integrates GNSS
velocities, leveling-based vertical velocities, and Glacial Isostatic
Adjustment (GIA) corrections to generate a unified vertical-velocity
field and associated visualizations.

The primary MATLAB script is:

    combine_vertical_velocities.m

This script produces an output file:

    Cascadia_vertical_combined.txt

containing longitude, latitude, vertical velocity (mm/yr), and
uncertainty (mm/yr).

------------------------------------------------------------------------

## 1. Included Data Sources and Required Citations

### **GNSS Velocities (MIDAS)**

Blewitt, G., Kreemer, C., Hammond, W. C., & Gazeaux, J. (2016).\
*MIDAS robust trend estimator for accurate GPS station velocities
without step detection.*\
*JGR Solid Earth, 121*, 2054--2068.\
https://doi.org/10.1002/2015JB012552\
(Data from: http://geodesy.unr.edu)

### **Leveling-Based Vertical Velocities**

#### Burgette et al. (2009)

Burgette, R. J., Weldon II, R. J., & Schmidt, D. A. (2009).\
*Interseismic uplift rates for western Oregon and along-strike variation
in locking on the Cascadia subduction zone.*\
*JGR Solid Earth, 114*, B01408.

#### Newton et al. (2021)

Newton, T. J., Weldon, R., Miller, I. M., Schmidt, D., Mauger, G.,
Morgan, H., & Grossman, E. (2021).\
*An assessment of vertical land movement to support coastal hazards
planning in Washington State.*\
*Water, 13*, 281.

### **GIA / Long-Wavelength Vertical Corrections**

Lau, N., Blewitt, G., & Becker, T. W. (2020).\
*Present-day crustal vertical velocity field for the contiguous United
States.*\
*JGR Solid Earth, 125*, e2020JB020066.

------------------------------------------------------------------------

## 2. Repository Structure

    /
    ├── Burgett_data/
    ├── Newton_Supplemental/
    ├── MIDAS/
    ├── Lau/
    ├── dem/
    ├── tl_2017_us_state/
    ├── tools/
    └── combine_vertical_velocities.m

------------------------------------------------------------------------

## 3. Overview of Workflow

1.  **Load MIDAS GNSS Velocities**\
2.  **Load Leveling Data (Newton and Burgette)**\
3.  **Reference Frame Alignment**\
4.  **Outlier Filtering**\
5.  **Spatial Homogenization**\
6.  **Apply Lau et al. (2020) GIA Corrections**\
7.  **Generate Maps and Profiles**\
8.  **Write Final Dataset**

------------------------------------------------------------------------

## 4. Requirements

-   MATLAB R2020+\
-   netCDF support\
-   DEM and mapping utilities included in `/dem`

------------------------------------------------------------------------

## 5. Contact / Usage Notes

If you use any component of the repository, please ensure that all
required citations (listed above) appear in your work.
