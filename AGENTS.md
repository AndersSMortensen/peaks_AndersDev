# AGENTS.md Specification & Template

This document defines the **AGENTS.md** format — a convention for scientific
Python packages to describe their tools in a way that AI agent systems can
discover and use automatically.

## For Library Authors

To make your package's tools discoverable, add an `AGENTS.md` file to your
package's root directory (next to `__init__.py`). The Science Pipeline (and
any other agent system that adopts this convention) will find it automatically
when your package is installed.

**You know your tools best.** We encourage library authors to write their own
AGENTS.md rather than relying on downstream systems to guess tool capabilities
from docstrings or API surface. Accurate descriptions lead to better AI-driven
tool selection.

### File location

```
your_package/
  __init__.py
  AGENTS.md          <-- here
  core/
  ...
```

## Format Specification

The file uses structured Markdown. Each `### heading` defines one tool.
Fields are specified as `- key: value` lines below the heading.

### Required fields

- `description` — What the tool does (1-2 sentences, plain language)
- `function` — How to call it: `module.path:function_name` or accessor notation

### Optional fields

- `parameters:` — Block of parameter definitions (indented sub-items)
- `returns` — What the function returns
- `when_to_use` — Conditions or scenarios for applying this tool
- `addresses` — Comma-separated list of artefacts/problems this tool corrects
- `effects` — Comma-separated list of outcomes
- `family` — Operation family (e.g. `denoising`, `calibration`, `visualization`)
- `strategy` — Algorithm within the family (e.g. `gaussian`, `binning`)
- `cost` — `low`, `medium`, `high`

### Parameter format

```markdown
- parameters:
  - param_name: type (required|optional, default=X) — description
```

## Example (peaks-arpes)

The following is a **template** showing how the [peaks](https://github.com/phrgab/peaks)
package could describe its tools. Library authors should write and maintain
this file — these descriptions are illustrative, not authoritative.

## Tools

### bin_data
- description: Bin (downsample) data along one or more dimensions by integer factors. Reduces noise at the cost of resolution.
- function: xr.DataArray accessor .bin_data()
- parameters:
  - eV: int (optional, default=1) — bin factor along energy axis
  - theta_par: int (optional, default=1) — bin factor along parallel angle
  - x1: int (optional, default=1) — bin factor along first spatial axis
  - x2: int (optional, default=1) — bin factor along second spatial axis
  - t: int (optional, default=1) — bin factor along time/delay axis
- returns: xr.DataArray with reduced dimensions
- when_to_use: Low SNR data, large datasets needing size reduction, quick-look analysis
- addresses: noise, large_dataset
- effects: reduces_noise, reduces_resolution
- family: denoising
- strategy: binning
- cost: low

### smooth
- description: Gaussian smoothing along energy and/or angle dimensions. FWHM given in axis units.
- function: xr.DataArray accessor .smooth()
- parameters:
  - eV: float (optional) — FWHM of Gaussian kernel in eV
  - theta_par: float (optional) — FWHM in degrees along parallel angle
- returns: xr.DataArray with smoothed intensity
- when_to_use: Noisy data where binning would lose too many points, before derivative analysis
- addresses: noise
- effects: reduces_noise, reduces_resolution
- family: denoising
- strategy: gaussian
- cost: low

### k_convert
- description: Convert angular ARPES data to momentum (k-space). Requires photon energy metadata.
- function: xr.DataArray accessor .k_convert()
- parameters:
  - V0: float (optional, default=10) — inner potential in eV for kz conversion
- returns: xr.DataArray with k_par (or kx, ky) replacing theta_par
- when_to_use: Band mapping, Fermi surface mapping, any analysis requiring momentum coordinates
- addresses: angular_coordinates
- effects: converts_to_momentum
- family: coordinate_transform
- strategy: free_electron
- cost: medium

### estimate_EF
- description: Estimate the Fermi level position by fitting a Fermi-Dirac distribution to the angle-integrated spectrum.
- function: xr.DataArray accessor .estimate_EF()
- returns: float — Fermi level position in eV
- when_to_use: Energy calibration, before gap extraction, when Fermi level position is unknown
- addresses: energy_calibration
- effects: identifies_fermi_level
- family: calibration
- strategy: fermi_dirac_fit
- cost: medium

### d2E
- description: Second derivative along energy axis. Enhances band features and suppresses smooth background.
- function: xr.DataArray accessor .d2E()
- returns: xr.DataArray with second derivative values
- when_to_use: Band structure visualisation, feature enhancement. Apply LAST — destroys quantitative intensity.
- addresses: feature_visibility
- effects: enhances_features, destroys_intensity
- family: visualization
- strategy: second_derivative
- cost: low

### curvature
- description: Curvature analysis for enhanced band visualisation. Better feature preservation than second derivative.
- function: xr.DataArray accessor .curvature()
- parameters:
  - a0_eV: float (optional, default=1.0) — free parameter for energy axis
  - a0_k: float (optional, default=1.0) — free parameter for momentum axis
- returns: xr.DataArray with curvature values
- when_to_use: Publication-quality band visualisation, alternative to second derivative
- addresses: feature_visibility
- effects: enhances_features, destroys_intensity
- family: visualization
- strategy: curvature
- cost: medium

### fermi_dirac_fit
- description: Fit a Fermi-Dirac distribution to the energy spectrum to extract Fermi level position, temperature, and resolution.
- function: peaks.core.fitting.fit:fermi_dirac_fit
- parameters:
  - temperature_guess: float (optional, default=300) — initial temperature guess in K
- returns: FitResult with position, width, temperature, resolution
- when_to_use: Precise energy calibration, resolution determination, temperature measurement from Fermi edge
- addresses: energy_calibration, resolution
- effects: measures_fermi_level, measures_resolution
- family: calibration
- strategy: fermi_dirac_fit
- cost: medium

### symmetrise
- description: Symmetrise ARPES data around a specified k-point or energy. Used to enforce crystal symmetry in Fermi surface maps.
- function: peaks.core.process.tools:symmetrise
- parameters:
  - center: float (optional, default=0) — center point for symmetrisation
  - axis: str (optional, default=k_par) — axis along which to symmetrise
- returns: xr.DataArray with symmetrised intensity
- when_to_use: Fermi surface analysis, when crystal symmetry should be enforced for visualisation
- addresses: asymmetry
- effects: enforces_symmetry
- family: symmetry
- strategy: mirror
- cost: low

### mdc_extract
- description: Extract Momentum Distribution Curve (MDC) at a specified energy. Reduces 2D data to 1D momentum profile.
- function: peaks.core.process.tools:mdc_extract
- parameters:
  - energy: float (required) — energy value in eV for the MDC cut
  - integration_window: float (optional, default=0.02) — integration width in eV around the cut energy
- returns: xr.DataArray — 1D momentum profile
- when_to_use: Peak fitting analysis, band dispersion extraction, linewidth analysis
- addresses: dimensionality_reduction
- effects: extracts_momentum_profile
- family: slicing
- strategy: energy_cut
- cost: low

### edc_extract
- description: Extract Energy Distribution Curve (EDC) at a specified momentum. Reduces 2D data to 1D energy profile.
- function: peaks.core.process.tools:edc_extract
- parameters:
  - momentum: float (required) — momentum value for the EDC cut
  - integration_window: float (optional, default=0.02) — integration width in inverse Angstroms
- returns: xr.DataArray — 1D energy profile
- when_to_use: Gap extraction, lineshape analysis, self-energy analysis
- addresses: dimensionality_reduction
- effects: extracts_energy_profile
- family: slicing
- strategy: momentum_cut
- cost: low

### normalize_by_region
- description: Normalise intensity by the mean or integrated intensity of a specified region (e.g. above Fermi level).
- function: peaks.core.process.tools:normalize_by_region
- parameters:
  - region: str (required) — region specification, e.g. 'above_ef' or 'eV>0.1'
  - method: str (optional, default=mean) — normalisation method: 'mean', 'max', 'sum'
- returns: xr.DataArray with normalised intensity
- when_to_use: Before comparing spectra, correcting for varying photon flux, preparing for publication
- addresses: intensity_variation
- effects: normalises_intensity
- family: normalization
- strategy: region_based
- cost: low
