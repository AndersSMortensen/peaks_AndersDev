# peaks AGENTS.md

`peaks` is a Python package for loading, processing, and analysing
Angle-Resolved Photoemission Spectroscopy (ARPES) data. Data are represented
as `xarray.DataArray` objects and most tools are available both as standalone
functions and as accessor methods (e.g. `data.smooth(...)`).

## Conventions

- **Accessor notation**: all tools listed here are registered as methods on
  `xr.DataArray` (and where applicable `xr.Dataset` / `xr.DataTree`). Call
  them as `data.<tool>(...)`.
- **Axis names**: `eV` (energy), `theta_par` / `k_par` / `kx` / `ky` / `kz`
  (angle or momentum), `x1` / `x2` (spatial), `hv` (photon energy).
- **Units**: axes carry `pint` units. Smoothing/broadening FWHMs should be
  given in axis units (or as `pint.Quantity`).
- **Analysis history**: every tool appends a record to `data.history` for
  provenance tracking.

## Tools

### load
- description: Load one or more ARPES data files into an `xr.DataArray` (single file) or `xr.DataTree` (multiple files). Supports `.ibw`, `.nxs`, `.zip`, `.nc`, `.zarr`, and other beamline formats. Metadata is parsed automatically.
- function: peaks.core.fileIO.data_loading:load
- parameters:
  - fpath: str or list (required) — file path(s) to load
  - lazy: bool (optional, default=None) — load lazily with dask; auto-determined by file size if None
  - loc: str (optional, default=None) — beamline/location identifier; auto-detected if None
  - metadata: bool (optional, default=True) — parse metadata into DataArray attributes
  - parallel: bool (optional, default=False) — load multiple files in parallel (h5/nxs only)
  - quiet: bool (optional, default=False) — suppress warnings
- returns: xr.DataArray, xr.Dataset, or xr.DataTree
- when_to_use: Start of every analysis workflow
- family: io
- strategy: file_load
- cost: low

### save
- description: Save a DataArray or Dataset as a NetCDF file, or a DataTree as a Zarr store. Peaks metadata is serialised into the file attributes where possible.
- function: peaks.core.fileIO.data_saving:save
- parameters:
  - data: xr.DataArray, xr.Dataset, or xr.DataTree (required) — data to save
  - fpath: str (required) — output file path (extension optional; .nc for DataArray/Dataset, .zarr for DataTree)
- returns: None
- when_to_use: Saving processed data for later reuse or sharing
- family: io
- strategy: file_save
- cost: low

### bin_data
- description: Bin (coarsen) data along one or more dimensions by integer factors. A thin wrapper around `xr.DataArray.coarsen` that also updates analysis history.
- function: xr.DataArray accessor .bin_data()
- parameters:
  - binning: int (optional, default=None) — apply the same bin factor to all dimensions; takes priority over per-dimension kwargs
  - boundary: str (optional, default='trim') — boundary handling: 'trim', 'exact', or 'pad'
  - **binning_kwargs: int — per-dimension bin factors, e.g. eV=3, theta_par=2
- returns: xr.DataArray with reduced dimensions
- when_to_use: Low SNR data, large datasets needing size reduction, quick-look analysis
- addresses: noise, large_dataset
- effects: reduces_noise, reduces_resolution
- family: denoising
- strategy: binning
- cost: low

### bin_spectra
- description: Shortcut to bin_data that automatically identifies and bins the spectral dimensions (eV and theta_par/k_par) by a uniform factor, leaving spatial dimensions unchanged.
- function: xr.DataArray accessor .bin_spectra()
- parameters:
  - binning: int (optional, default=2) — bin factor for spectral dimensions
  - boundary: str (optional, default='trim') — boundary handling
- returns: xr.DataArray with binned spectral dimensions
- when_to_use: Quick noise reduction when the spectral dimensions are known but spatial ones should be preserved
- addresses: noise
- effects: reduces_noise, reduces_resolution
- family: denoising
- strategy: binning
- cost: low

### smooth
- description: Gaussian smoothing along one or more dimensions. FWHM values are specified in axis units (or as `pint.Quantity`); dimensions not listed are left unchanged.
- function: xr.DataArray accessor .smooth()
- parameters:
  - **smoothing_kwargs: float or pint.Quantity — FWHM per axis, e.g. eV=0.05, theta_par=0.3
- returns: xr.DataArray with smoothed intensity
- when_to_use: Noisy data where binning would lose too many points; always apply before derivative/curvature analysis
- addresses: noise
- effects: reduces_noise, reduces_resolution
- family: denoising
- strategy: gaussian
- cost: low

### norm
- description: Normalise data intensity. Can normalise to unity (max), by the mean, by an integrated distribution curve along a given dimension, or by the integrated intensity of an arbitrary slice/ROI.
- function: xr.DataArray accessor .norm()
- parameters:
  - dim: str (optional, default=None) — 'all' to normalise by global mean; axis name (e.g. 'eV') to normalise by an integrated DC along that direction; None to normalise to unity
  - **kwargs: slice — coordinate slices defining a ROI to normalise by, e.g. eV=slice(-0.5, -0.3)
- returns: xr.DataArray with normalised intensity
- when_to_use: Before comparing spectra across different scans or photon energies; correcting for varying photon flux
- addresses: intensity_variation
- effects: normalises_intensity
- family: normalization
- strategy: region_based
- cost: low

### bgs
- description: Background subtraction. Supports subtraction of a constant, mean, integrated distribution curve along a given dimension, a region-defined background, or a Shirley background (for XPS).
- function: xr.DataArray accessor .bgs()
- parameters:
  - subtraction: float, pint.Quantity, or str (optional, default=None) — value to subtract; 'all' for mean; axis name for integrated DC; 'Shirley' for iterative Shirley background
  - num_avg: int (optional, default=1) — Shirley: number of points to average at data endpoints
  - offset_start: float (optional, default=0) — Shirley: offset subtracted from start value
  - offset_end: float (optional, default=0) — Shirley: offset subtracted from end value
  - max_iterations: int (optional, default=10) — Shirley: maximum convergence iterations
  - **kwargs: slice — coordinate slices defining the background region
- returns: xr.DataArray with background subtracted
- when_to_use: Removing a slowly varying or step-like background before fitting or derivative analysis
- addresses: background
- effects: removes_background
- family: background_subtraction
- strategy: region_based, shirley
- cost: low

### k_convert
- description: Convert angular ARPES data to momentum (k-space) using the free-electron final-state model. Handles 2D dispersions, 3D Fermi maps, and hν scans (producing kz). Reads analyser geometry and photon energy from metadata.
- function: xr.DataArray accessor .k_convert()
- parameters:
  - eV: slice (optional) — binding energy range for the output, e.g. slice(-1.0, 0.1, 0.002)
  - eV_slice: tuple (optional) — (energy, width) to return a single MDC-like energy slice after conversion
  - kx: slice (optional) — kx range for the output
  - ky: slice (optional) — ky range for the output
  - kz: slice (optional) — kz range for hν scans
  - return_kz_scan_in_hv: bool (optional, default=False) — return hν-scan data as (hv, eV, k_∥) instead of (kz, eV, k_∥)
  - quiet: bool (optional, default=False) — suppress progress bar and warnings
- returns: xr.DataArray in k-space (k_par/kx/ky replaces theta_par; kz replaces hv for hν scans)
- when_to_use: Band mapping, Fermi surface mapping, kz dispersion; any analysis requiring momentum coordinates
- addresses: angular_coordinates
- effects: converts_to_momentum
- family: coordinate_transform
- strategy: free_electron_final_state
- cost: medium

### estimate_EF
- description: Quickly estimate the Fermi level position from the peak of the first derivative of the angle-integrated spectrum. For hν scans, fits a polynomial correction. This is an approximate method suitable for seeding fit routines or GUIs, not for publication-quality determination.
- function: xr.DataArray accessor .estimate_EF()
- returns: float (single scan) or numpy.ndarray (hν scan) — estimated Fermi level in eV
- when_to_use: Quick energy calibration check; seed value for fit_gold or EF correction; when Fermi level position is unknown
- addresses: energy_calibration
- effects: identifies_fermi_level
- family: calibration
- strategy: derivative_peak
- cost: low

### fit_gold
- description: Fit a gold reference spectrum to a LinearDosFermiModel (Fermi-Dirac × linear DOS + background, convolved with a Gaussian) to extract EF, temperature, and instrumental energy resolution. For 2D data (gold angle map), also fits and returns a polynomial EF correction as a function of angle.
- function: peaks.core.fitting.fit:fit_gold
- parameters:
  - data: xr.DataArray (required) — 1D or 2D gold reference spectrum
  - EF_correction_type: str (optional, default='poly4') — type of EF correction for 2D data: 'poly4', 'poly3', 'quadratic', 'linear', or 'average'
  - **kwargs — initial parameter overrides (e.g. bg_slope=0)
- returns: xr.Dataset with fit parameters (EF, T, sigma_conv, DOS coefficients, background coefficients), uncertainties, lmfit ModelResult objects, and EF_correction attribute
- when_to_use: Determining instrumental energy resolution; obtaining a precise EF correction before applying to sample data
- addresses: energy_calibration, resolution
- effects: measures_fermi_level, measures_resolution
- family: calibration
- strategy: fermi_dirac_fit
- cost: medium

### fit
- description: Fit an `lmfit.Model` to a DataArray, broadcasting over non-independent dimensions. Sequential fitting (updating initial parameters from the previous slice) is supported for 2D data, which greatly improves convergence along a slowly varying axis.
- function: xr.DataArray accessor .fit()
- parameters:
  - model: lmfit.Model (required) — model to fit
  - params: lmfit.Parameters (required) — initial parameters
  - independent_var: str (optional, default=None) — dimension name for the independent variable; required for 2D+ data
  - sequential: bool (optional, default=True) — update initial params from previous slice when fitting 2D data
  - reverse_sequential_fit_order: bool (optional, default=False) — reverse the order of sequential fitting
- returns: xr.Dataset with best-fit parameter values, uncertainties, and lmfit ModelResult objects
- when_to_use: Peak fitting (MDCs, EDCs), lineshape analysis, extracting band positions and widths
- addresses: peak_extraction
- effects: extracts_parameters
- family: fitting
- strategy: lmfit
- cost: medium

### MDC
- description: Extract one or more Momentum Distribution Curves (MDC) from a dispersion at specified energy value(s), with optional energy integration window.
- function: xr.DataArray accessor .MDC()
- parameters:
  - E: float, list, tuple, or numpy.ndarray (optional, default=0) — energy value(s); tuple format (start, end, step) for a range
  - dE: float (optional, default=0) — total integration window in eV (integrates ±dE/2)
- returns: xr.DataArray — 1D (or stacked 2D) MDC(s)
- when_to_use: Band dispersion extraction, linewidth analysis, Fermi surface extraction from a 3D cube
- addresses: dimensionality_reduction
- effects: extracts_momentum_profile
- family: slicing
- strategy: energy_cut
- cost: low

### EDC
- description: Extract one or more Energy Distribution Curves (EDC) from a dispersion at specified momentum/angle value(s), with optional integration window.
- function: xr.DataArray accessor .EDC()
- parameters:
  - k: float, list, tuple, or numpy.ndarray (optional, default=0) — k or theta_par value(s); tuple format (start, end, step) for a range
  - dk: float (optional, default=0) — total integration window in axis units (integrates ±dk/2)
- returns: xr.DataArray — 1D (or stacked 2D) EDC(s)
- when_to_use: Gap measurement, self-energy analysis, lineshape fitting at fixed momentum
- addresses: dimensionality_reduction
- effects: extracts_energy_profile
- family: slicing
- strategy: momentum_cut
- cost: low

### DC
- description: General distribution curve extraction along any coordinate, at one or more values with an optional integration window. MDC and EDC are convenience wrappers around this function.
- function: xr.DataArray accessor .DC()
- parameters:
  - coord: str (optional, default='eV') — coordinate to cut along
  - val: float, list, tuple, or numpy.ndarray (optional, default=0) — value(s) to extract at; tuple (start, end, step) for a range
  - dval: float (optional, default=0) — total integration window in coordinate units
- returns: xr.DataArray — extracted DC(s)
- when_to_use: Extracting cuts along non-standard axes (e.g. hv, x1, x2, ana_polar)
- addresses: dimensionality_reduction
- effects: extracts_profile
- family: slicing
- strategy: coordinate_cut
- cost: low

### DOS
- description: Integrate over all non-energy dimensions to return the best approximation to the density of states available from the data.
- function: xr.DataArray accessor .DOS()
- returns: xr.DataArray — 1D energy spectrum (DOS)
- when_to_use: Overview of spectral weight; seed for estimate_EF; before Fermi-Dirac fitting
- addresses: dimensionality_reduction
- effects: extracts_dos
- family: slicing
- strategy: full_integration
- cost: low

### tot
- description: Integrate spatial map data over all non-spatial (energy + angle/k) dimensions to produce a real-space intensity image, or over all spatial dimensions to produce a dispersion.
- function: xr.DataArray accessor .tot()
- parameters:
  - spatial_int: bool (optional, default=False) — if True, integrate over spatial (x1, x2) dimensions instead
- returns: xr.DataArray — spatially resolved or spatially integrated map
- when_to_use: Spatial maps (nanoARPES): obtain total intensity image or spatially averaged dispersion
- addresses: dimensionality_reduction
- effects: extracts_spatial_map, extracts_dispersion
- family: slicing
- strategy: full_integration
- cost: low

### deriv
- description: Differentiate data along one or more specified dimensions, in sequence. General-purpose wrapper around `xr.DataArray.differentiate` that preserves metadata and analysis history.
- function: xr.DataArray accessor .deriv()
- parameters:
  - dims: str or list of str (required) — dimension(s) to differentiate along, in order (repeat to obtain higher-order derivatives)
- returns: xr.DataArray with differentiated values
- when_to_use: Custom derivative orders or mixed derivatives not covered by the shortcut functions
- addresses: feature_visibility
- effects: enhances_features, destroys_intensity
- family: visualization
- strategy: finite_difference
- cost: low

### d2E
- description: Shortcut for double differentiation along the energy (eV) axis. Enhances band features and suppresses slowly varying background.
- function: xr.DataArray accessor .d2E()
- returns: xr.DataArray with d²I/dE² values
- when_to_use: Band structure visualisation. Apply AFTER smoothing and LAST in a processing chain — destroys quantitative intensity.
- addresses: feature_visibility
- effects: enhances_features, destroys_intensity
- family: visualization
- strategy: second_derivative
- cost: low

### d2k
- description: Shortcut for double differentiation along the momentum (or angle) dimension. Works whether data is in angle or k-space.
- function: xr.DataArray accessor .d2k()
- returns: xr.DataArray with d²I/dk² values
- when_to_use: Enhancing features along the momentum direction; complement to d2E
- addresses: feature_visibility
- effects: enhances_features, destroys_intensity
- family: visualization
- strategy: second_derivative
- cost: low

### dEdk
- description: Shortcut for sequential differentiation: first along eV, then along the momentum/angle dimension.
- function: xr.DataArray accessor .dEdk()
- returns: xr.DataArray with d²I/dEdk values
- when_to_use: Mixed derivative enhancement; highlights features that vary in both energy and momentum
- addresses: feature_visibility
- effects: enhances_features, destroys_intensity
- family: visualization
- strategy: mixed_derivative
- cost: low

### dkdE
- description: Shortcut for sequential differentiation: first along the momentum/angle dimension, then along eV.
- function: xr.DataArray accessor .dkdE()
- returns: xr.DataArray with d²I/dkdE values
- when_to_use: Mixed derivative enhancement; order-reversed complement to dEdk
- addresses: feature_visibility
- effects: enhances_features, destroys_intensity
- family: visualization
- strategy: mixed_derivative
- cost: low

### curvature
- description: 2D curvature analysis following Zhang et al., Rev. Sci. Instrum. 82, 043712 (2011). Better preserves feature sharpness than the second derivative. Free parameters must be provided for both axes of the 2D data; set a parameter to 0 for 1D curvature along the other axis.
- function: xr.DataArray accessor .curvature()
- parameters:
  - **parameter_kwargs: float — free parameters per axis in format axis=value, e.g. eV=1, theta_par=10; both axes required
- returns: xr.DataArray with curvature values
- when_to_use: Publication-quality band visualisation; preferred alternative to d2E for preserving lineshape; apply after smoothing
- addresses: feature_visibility
- effects: enhances_features, destroys_intensity
- family: visualization
- strategy: curvature
- cost: medium

### min_gradient
- description: Minimum gradient analysis following He et al., Rev. Sci. Instrum. 88, 073903 (2017). Computes the ratio of the data to the Gaussian-filtered gradient magnitude; enhances band features while being less sensitive to noise than derivatives.
- function: xr.DataArray accessor .min_gradient()
- parameters:
  - **smoothing_kwargs: float — Gaussian filter FWHM per axis in axis units, e.g. eV=0.05, theta_par=0.2; axes not listed default to 1 pixel
- returns: xr.DataArray with minimum gradient values
- when_to_use: Band visualisation in noisy data; robust alternative to curvature when derivative methods fail
- addresses: feature_visibility, noise
- effects: enhances_features
- family: visualization
- strategy: minimum_gradient
- cost: medium

### sym
- description: Mirror-symmetrise data about a specified coordinate value along one axis. The result is the sum of the original and its mirror image. Can also return just the flipped copy.
- function: xr.DataArray accessor .sym()
- parameters:
  - flipped: bool (optional, default=False) — return only the mirrored data rather than original + mirror
  - fillna: bool (optional, default=True) — fill NaNs in non-overlapping regions with 0
  - **sym_kwarg: float or pint.Quantity — axis and centre value, e.g. theta_par=1.4 or eV=0; only one axis supported; defaults to eV=0
- returns: xr.DataArray — symmetrised (or flipped) data
- when_to_use: Enforcing time-reversal symmetry in dispersions; symmetrising about normal emission; preparing data for Fermi surface analysis
- addresses: asymmetry
- effects: enforces_symmetry
- family: symmetry
- strategy: mirror
- cost: low

### sym_nfold
- description: Apply n-fold rotational symmetry to 2D or 3D data about a specified centre. Useful for enforcing crystal point-group symmetry in Fermi surface maps.
- function: xr.DataArray accessor .sym_nfold()
- parameters:
  - nfold: int (required) — rotational order (e.g. 4 for C4, 6 for C6)
  - expand: bool (optional, default=True) — expand the coordinate grid to show all symmetrised data
  - fillna: bool (optional, default=True) — fill NaNs from rotation with 0 before summing
  - **centre_kwargs: float — centre of rotation per axis, e.g. kx=0.1, ky=0.0; defaults to (0, 0)
- returns: xr.DataArray — n-fold symmetrised data
- when_to_use: Fermi surface maps with known crystal symmetry; increasing effective statistics by exploiting symmetry
- addresses: asymmetry, statistics
- effects: enforces_symmetry, improves_statistics
- family: symmetry
- strategy: rotational
- cost: medium

### estimate_sym_point
- description: Estimate the symmetry centre (e.g. normal emission, Γ point) of data by phase cross-correlation of the data with its mirror image. Returns a dict of {axis: centre_value}.
- function: xr.DataArray accessor .estimate_sym_point()
- parameters:
  - dims: str, list, or tuple (optional, default=None) — dimensions to estimate the centre along; defaults to all dimensions
  - upsample_factor: int (optional, default=100) — subpixel precision factor for phase cross-correlation
- returns: dict mapping dimension name to estimated centre coordinate value
- when_to_use: Before symmetrisation; finding normal emission in a dispersion; locating the Γ point in a Fermi surface map
- addresses: alignment
- effects: identifies_centre
- family: alignment
- strategy: phase_cross_correlation
- cost: low

### rotate
- description: Rotate 2D (or 3D) data by a specified angle about a given centre of rotation, using bilinear interpolation.
- function: xr.DataArray accessor .rotate()
- parameters:
  - rotation: float (required) — rotation angle in degrees (anticlockwise)
  - **centre_kwargs: float — centre of rotation per axis, e.g. kx=0, ky=0; defaults to data midpoint
- returns: xr.DataArray — rotated data on the original coordinate grid
- when_to_use: Correcting angular misalignment in Fermi surface maps; aligning crystal axes with k-space axes
- addresses: misalignment
- effects: corrects_rotation
- family: alignment
- strategy: rotation
- cost: medium

### degrid
- description: Remove a periodic mesh grid artefact from 2D data by identifying and zeroing high-intensity spikes in the FFT outside a protected central low-frequency region.
- function: xr.DataArray accessor .degrid()
- parameters:
  - width: float (optional, default=0.1) — fraction of total width to protect as the central low-frequency region
  - height: float (optional, default=0.1) — fraction of total height to protect as the central low-frequency region
  - cutoff: float (optional, default=4) — intensity threshold as a multiple of the FFT mean for spike removal
- returns: xr.DataArray — data with grid artefact removed
- when_to_use: nanoARPES or µ-ARPES data from instruments with a physical mesh; before spatial or spectral analysis
- addresses: grid_artefact
- effects: removes_artefact
- family: artefact_removal
- strategy: fft_filtering
- cost: low

### drop_nan_borders
- description: Trim all-NaN rows/columns from the edges of each dimension of a DataArray or Dataset.
- function: xr.DataArray accessor .drop_nan_borders()
- returns: xr.DataArray or xr.Dataset — trimmed data
- when_to_use: After k-conversion or rotation, which introduce NaN borders; before fitting or display
- addresses: nan_borders
- effects: trims_data
- family: cleaning
- strategy: border_trim
- cost: low

### drop_zero_borders
- description: Trim all-zero rows/columns from the edges of each dimension of a DataArray or Dataset.
- function: xr.DataArray accessor .drop_zero_borders()
- returns: xr.DataArray or xr.Dataset — trimmed data
- when_to_use: After operations that pad with zeros; before fitting or display
- addresses: zero_borders
- effects: trims_data
- family: cleaning
- strategy: border_trim
- cost: low

### radial_cuts
- description: Extract radial cuts of a 2D Fermi surface map (or 3D cube) as a function of azimuthal angle about a central point. Returns a 3D DataArray with dimensions (azimuthal_angle, radial_distance, eV).
- function: xr.DataArray accessor .radial_cuts()
- parameters:
  - num_azi: int (optional, default=361) — number of evenly spaced azimuthal angles between 0° and 360°
  - num_points: int (optional, default=200) — number of radial points per cut
  - radius: float (optional, default=2) — maximum radial distance in data axis units
  - **centre_kwargs: float — centre of the radial fan, e.g. kx=0.1, ky=0.0
- returns: xr.DataArray — radial cuts stacked along an azimuthal_angle dimension
- when_to_use: Warping square Fermi surfaces to circular geometry; azimuthal analysis of Fermi surface maps
- addresses: coordinate_transform
- effects: extracts_radial_profile
- family: slicing
- strategy: radial_cut
- cost: medium

### extract_cut
- description: Extract an arbitrary straight-line cut through 2D data between two specified points in coordinate space.
- function: xr.DataArray accessor .extract_cut()
- parameters:
  - start_point: tuple (required) — (coord1, coord2) of the cut start in data coordinate values
  - end_point: tuple (required) — (coord1, coord2) of the cut end in data coordinate values
  - num_points: int (optional, default=None) — number of points along the cut; defaults to the data resolution
- returns: xr.DataArray — 1D cut with a path-length coordinate
- when_to_use: Extracting off-axis cuts in Fermi surface maps; kz–kx or kz–ky cuts through 3D data
- addresses: dimensionality_reduction
- effects: extracts_profile
- family: slicing
- strategy: arbitrary_cut
- cost: low

### mask_data
- description: Mask a region of interest (ROI) defined by a polygon path, returning the masked data and optionally its integrated intensity within the ROI.
- function: xr.DataArray accessor .mask_data()
- parameters:
  - ROI: list of tuples (required) — vertices of the polygonal ROI in coordinate space, e.g. [(kx0,ky0), (kx1,ky1), ...]
  - return_integrated: bool (optional, default=True) — also return the integrated intensity within the ROI
- returns: xr.DataArray (masked data) or tuple (masked_data, integrated_intensity)
- when_to_use: Selecting specific Brillouin zone regions for analysis; excluding spurious regions from fitting
- addresses: region_selection
- effects: selects_region
- family: slicing
- strategy: polygon_mask
- cost: low
