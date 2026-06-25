# Source Code — 2D Image Restoration & Inpainting

MATLAB code implementing Tikhonov (TIK-L2) and Total Variation (TV-L2) regularization for image restoration, plus TIK-L2 for inpainting. Corresponds to **Chapter 1** of the report.

## Files

### Main Driver

**`MAIN_1_EX1.m`**
- Single entry point for all 2D experiments
- Loads/generates the corrupted image via `GENERATE_DATA_2D/`
- Sweeps a range of regularization parameters λ and tracks relative error / ISNR
- Controlled by three boolean flags at the top of the script:

```matlab
TEST_TV_L2_GRAD_DESC    = 1;   % TV-L2 restoration via Gradient Descent
TEST_TIK_L2_U_REST      = 0;   % TIK-L2 restoration (closed-form, FFT-based)
TEST_TIK_L2_U_INP       = 0;   % TIK-L2 inpainting
```

### Core Algorithm

**`TV_L2_U_GD.m`**

```matlab
function [xest, itrs, r] = TV_L2_U_GD(b, b_k, lambda, eps, x0, tau, itrs_max, ch_th, debug_cw, xtr)
```

Solves the TV-L2 restoration problem via **Gradient Descent**. Since the TV term `‖∇x‖₂` is non-differentiable at zero, a smoothed/differentiable surrogate (controlled by `eps`) is used so the gradient descent step is well defined.

| Argument | Meaning |
|---|---|
| `b` | Degraded (blurred + noisy) image |
| `b_k` | Blur kernel |
| `lambda` | Regularization parameter |
| `eps` | Smoothing parameter for differentiability of the TV term |
| `x0` | Initial condition |
| `tau` | Gradient descent step size |
| `itrs_max` | Maximum number of iterations |
| `ch_th` | Change threshold (stopping criterion on residual) |
| `xtr` | Ground-truth (noise-free) image, used to track relative error during iterations |

The blur operator is applied in the Fourier domain (`psf2otf`), since the blur matrix A is Block-Circulant-with-Circulant-Blocks (BCCB) and convolution can be computed efficiently via FFT.

The TIK-L2 restoration problem (the other main branch of `MAIN_1_EX1.m`) is solved in **closed form directly via FFT**, since for BCCB A the normal equations diagonalize in the Fourier domain — no iterative solver is needed in that case.

### Differential Operators

**`gradx.m` / `grady.m`**
- Discrete forward-difference horizontal/vertical gradient operators
- Used to build the discretized `∇x` needed by the TV-L2 regularization term

**`div.m`**
- Discrete divergence operator — the adjoint of the gradient
- Used inside `TV_L2_U_GD.m` to compute the gradient of the regularization term during the descent iterations

### Evaluation Utilities

**`compute_rel_err.m`**
```matlab
rel_err = compute_rel_err(x_est, x_true)
```
Computes `‖x_est - x_true‖ / ‖x_true‖`.

**`compute_snr.m`**
Computes the Signal-to-Noise Ratio / Improvement in SNR (ISNR) between the restored and the original image.

### Data Generation

**`GENERATE_DATA_2D/GENERATE_DATA_2D.m`**
- Loads a test image from `GENERATE_DATA_2D/input/`
- Applies Gaussian blur (configurable band/sigma) and AWGN (µ = 0, σ configurable, e.g. σ = 10 in the report)
- For inpainting experiments, applies a binary mask instead of (or in addition to) blur
- Produces the corrupted observation `b` consumed by `MAIN_1_EX1.m`

**`GENERATE_DATA_2D/compute_snr.m`** — local copy of the SNR utility used during data generation/visualization.

**`GENERATE_DATA_2D/input/`** — test images and masks (cameraman, parrot, QR-codes, satellite, checkboards, squirrel-related masks, etc.) — see `data/README.md` (project root `data/` if you split it out) or the input folder directly for the full list.

## Typical Usage

```matlab
cd src/2D
addpath('GENERATE_DATA_2D')

% Edit the flags inside MAIN_1_EX1.m to select which experiment to run,
% then simply run:
MAIN_1_EX1
```

### Calling the Solvers Directly

```matlab
% TV-L2 restoration via Gradient Descent
[x_est, iters, residual] = TV_L2_U_GD(b, blur_kernel, lambda, eps, x0, tau, max_iters, change_thresh, debug_flag, x_true);

% Relative error of the result
err = compute_rel_err(x_est, x_true);
```

## Notes on the Two Restoration Methods

| | TIK-L2 | TV-L2 |
|---|---|---|
| Regularization term | `‖Lx‖₂²` (quadratic) | `‖Lx‖₂` (non-smooth) |
| Solver | Closed-form via FFT | Gradient Descent (smoothed surrogate) |
| Best suited for | Smooth images | Piecewise-smooth images with sharp edges |
| Computational cost | Very low (single FFT) | Higher (iterative) |

See the main `README.md` for the experimental results obtained with each method on the Cameraman, QR-code, and egg-carton (sinusoidal) test images.
