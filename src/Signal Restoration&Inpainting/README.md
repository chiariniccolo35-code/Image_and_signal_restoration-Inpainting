# Source Code — 1D Signal Restoration & Inpainting

MATLAB code implementing Generalized Tikhonov regularization for 1D signal denoising, denoising+deblurring, and inpainting. Corresponds to **Chapter 2** of the report.

## Files

### Main Driver

**`MAIN_1_EX2.m`**
- Single entry point for all 1D signal experiments
- Loads/generates the corrupted signal via `GENERATE_DATA_1D/`
- Sweeps a range of regularization parameters λ and tracks relative error
- Controlled by three boolean flags at the top of the script:

```matlab
TEST_DENOISING   = 0;   % Signal denoising (AWGN only)
TEST_RESTORATION = 1;   % Signal denoising + deblurring (Gaussian blur + AWGN)
TEST_INPAINTING  = 0;   % Signal inpainting (missing samples)
```

All three experiments use a **Generalized Tikhonov** formulation with a **second-order finite-difference discretization matrix D** (penalizing signal curvature, rather than the first-order gradient used in the 2D image case):

- **Denoising:** `x*_λ = argmin_x { ‖x - b‖₂² + ‖Dx‖₂² }`
- **Denoising + Deblurring:** `x*_λ = argmin_x { ‖Ax - b‖₂² + ‖Dx‖₂² }` (A = blur operator, reflective boundary conditions)
- **Inpainting:** `v*_λ = argmin_v ‖D(Sᵗb + Sᶜᵗv)‖₂²` (recovering unknown samples v from known samples b via a sub-sampling matrix S)

### Data Generation

**`GENERATE_DATA_1D/GENERATE_DATA_1D.m`**
- Loads or synthesizes a 1D test signal
- Applies AWGN (µ = 0, σ = 0.05 in the report) and/or Gaussian blur and/or sub-sampling (for inpainting)
- Produces the corrupted observation `b` consumed by `MAIN_1_EX2.m`

**`GENERATE_DATA_1D/MakeSignal.m`**
- Utility for generating standard synthetic test signals (e.g. piecewise-smooth / Heaviside-like signals) commonly used in signal-processing benchmarks

**`GENERATE_DATA_1D/generate_A_from_b_k_1D.m`**
- Builds the 1D blurring operator matrix A from a given blur kernel `b_k`

**`GENERATE_DATA_1D/generate_b_k_1D.m`**
- Generates the 1D blur kernel itself (e.g. Gaussian kernel) used by `generate_A_from_b_k_1D.m`

**`GENERATE_DATA_1D/baart.m`, `GENERATE_DATA_1D/shaw.m`, `GENERATE_DATA_1D/phillips.m`**
- Classical 1D inverse-problem test-case generators (standard ill-posed problem benchmarks from the regularization literature), available as alternative signal sources

**`GENERATE_DATA_1D/wavread.m`**
- Utility for reading `.wav` audio files as 1D signals (see `input/sp1.wav`)

**`GENERATE_DATA_1D/input/`**
- Pre-generated test signals (`signal_1.mat` … `signal_6.mat`), an ECG recording (`ECG_data.mat`), inpainting test data (`inpainting_data.mat`), and an audio sample (`sp1.wav`)

## Typical Usage

```matlab
cd src/1D
addpath('GENERATE_DATA_1D')

% Edit the flags inside MAIN_1_EX2.m to select the experiment, then run:
MAIN_1_EX2
```

### Parameter Sweep Example

```matlab
lambdas_min = 10^-4;
lambdas_max = 30;
lambdas_n   = 500;
lambdas     = linspace(lambdas_min, lambdas_max, lambdas_n);
```

The script evaluates relative error across the full λ-range and reports the optimal value (see main `README.md` for the exact optimal λ and relative error obtained for each test signal in the report).

## Notes on the Signal Inpainting Formulation

For inpainting, the signal `x ∈ Rᴴ` is split into:
- **Known samples** `b` (K samples), selected via sub-sampling matrix `S ∈ Rᴷˣᴺ` (identity matrix with rows removed)
- **Unknown samples** `v` (N−K samples), recovered by solving the smoothness-regularized least-squares problem above, using `Sᶜ` (the complementary rows of the identity not in S)

The reconstructed signal is then assembled as `x = Sᵗb + Sᶜᵗv`.
