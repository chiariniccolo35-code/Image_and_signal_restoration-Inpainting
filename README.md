# Restoration and Inpainting in Image and Signal Processing

A MATLAB implementation of **Tikhonov (TIK-L2)** and **Total Variation (TV-L2)** regularization for image/signal **restoration** (denoising + deblurring) and **inpainting**.

## Project Overview

This project implements and compares two classical variational regularization methods — **Tikhonov** and **Total Variation** — applied to two different domains:

1. **Image restoration & inpainting** (2D): denoising + deblurring of images corrupted by Gaussian blur and Additive White Gaussian Noise (AWGN), and reconstruction of images with missing pixels (inpainting)
2. **Signal restoration & inpainting** (1D): denoising, denoising+deblurring, and inpainting of 1D signals

This is **distinct from the companion ADMM/super-resolution project** — here the regularization problems are solved either via **closed-form FFT-based solutions** (TIK-L2, since the blur operator is Block-Circulant-with-Circulant-Blocks) or via **Gradient Descent** (TV-L2, since the TV term is non-differentiable at 0 and requires a smoothed/differentiable approximation).

## Chapter 1 — Image Restoration

Images are corrupted by **Gaussian blur** and **AWGN (µ = 0, σ = 10)**.

### 1.1 — Tikhonov Restoration (TIK-L2)

```
x*_λ = argmin_x { ‖Ax - b‖_2^2 + λ‖Lx‖_2^2 },   λ ∈ R+
```

- **A** — blurring matrix (solved via FFT, since A is Block-Circulant-with-Circulant-Blocks)
- **L** — first-order derivative discretization matrix
- **λ** — regularization parameter (trade-off between fidelity and smoothness)

Applied to three test images: Cameraman, QR-code, and a sinusoidal/"egg-carton" image.

### 1.2 — Total Variation Restoration (TV-L2)

```
x*_λ = argmin_x { ‖Ax - b‖_2^2 + λ‖Lx‖_2 },   λ ∈ R+   (solved via Gradient Descent)
```

Unlike Tikhonov, the TV regularization term is **convex but non-differentiable**; it is solved here via **Gradient Descent** with a smoothed (differentiable) approximation of the TV norm (`TV_L2_U_GD.m`), applied to the same three test images for direct comparison with TIK-L2.

### 1.3 — Image Inpainting

```
x*_λ = argmin_x { ‖Sx - b‖_2^2 + λ‖Lx‖_2^2 },   λ ∈ R+
```

Where **S = diag(s₁₁, ..., sₙₙ)** is the binary inpainting mask (sᵢᵢ = 1 if pixel i is known, 0 if it belongs to the inpainting/missing region).

Two families of experiments are carried out:
- **Object removal** (Squirrel — remove a white thread; Parrot — remove a cage), each tested both noise-free and with AWGN
- **Random-mask reconstruction** on the Cameraman image, comparing three mask topologies: Random Point Sampling, Random Horizontal Strip Sampling, and Random Vertical Strip Sampling

## Chapter 2 — Signal Restoration

All signal experiments use **Generalized Tikhonov** regularization with AWGN (µ = 0, σ = 0.05), and a **second-order derivative** discretization matrix D (penalizing signal curvature rather than the gradient as in the image case).

### 2.1 — Signal Denoising

```
x*_λ = argmin_x { ‖x - b‖_2^2 + ‖Dx‖_2^2 }
```

Applied to four synthetic 1D test signals.

### 2.2 — Signal Denoising + Deblurring

```
x*_λ = argmin_x { ‖Ax - b‖_2^2 + ‖Dx‖_2^2 }
```

Same signals as above, now also blurred (Gaussian blur, reflective boundary conditions on A).

### 2.3 — Signal Inpainting

Given only **K < N** known samples **b** of a signal **x ∈ Rᴴ**, with **S ∈ Rᴷˣᴺ** the (row-reduced identity) sub-sampling matrix selecting the known samples, the signal is written as:

```
x = Sᵗ b + Sᶜᵗ v
```

where **v** (the unknown samples) is recovered by solving:

```
v*_λ = argmin_v ‖D(Sᵗb + Sᶜᵗv)‖_2^2
```

with **D** the second-order finite-difference discretization matrix (smoothness prior on the reconstructed signal).

## Project Structure

```
.
├── README.md                          # This file
├── docs/
│   ├── REPORT_LAB_1_IMAGING.pdf      # Full technical report (Chapters 1-2)
│   └── Lab_1.pdf                      # Original assignment specification
│
├── src/
│   ├── 2D/                            # Image restoration & inpainting (Chapter 1)
│   │   ├── MAIN_1_EX1.m               # Main script: TIK-L2, TV-L2 (restoration & inpainting)
│   │   ├── TV_L2_U_GD.m               # TV-L2 solver via Gradient Descent
│   │   ├── compute_rel_err.m          # Relative error metric
│   │   ├── compute_snr.m              # SNR/ISNR metric
│   │   ├── gradx.m / grady.m          # Discrete horizontal/vertical gradient operators
│   │   ├── div.m                      # Discrete divergence operator (adjoint of gradient)
│   │   └── GENERATE_DATA_2D/          # Data generation: loads images, applies blur+noise+mask
│   │
│   └── 1D/                            # Signal restoration & inpainting (Chapter 2)
│       ├── MAIN_1_EX2.m               # Main script: denoising, deblurring, inpainting
│       └── GENERATE_DATA_1D/          # Data generation: loads/generates 1D signals, applies blur+noise
│
└── results/
    └── IMAGE_RESTORATION_PICTURES/
        ├── TIK_L2/                    # Tikhonov restoration results (Cameraman, Parrot, QR-code, Satellite, Sinusoid)
        ├── TV-RESTORATION/            # Total Variation restoration results (same images)
        └── INPAITING/                 # Inpainting results (Cameraman masks, Parrot, QR-code, Squirrel)
```

## Core MATLAB Files

| File | Role |
|------|------|
| `MAIN_1_EX1.m` | Main driver for **2D** experiments: TIK-L2 restoration (FFT-based), TV-L2 restoration (Gradient Descent), TIK-L2 inpainting |
| `MAIN_1_EX2.m` | Main driver for **1D** experiments: signal denoising, denoising+deblurring, inpainting |
| `TV_L2_U_GD.m` | Implements the **Gradient Descent** solver for the TV-L2 restoration problem |
| `gradx.m`, `grady.m` | Discrete horizontal/vertical gradient operators (forward differences) |
| `div.m` | Discrete divergence operator — adjoint of the gradient, used inside TV-L2 GD |
| `compute_rel_err.m` | Computes the relative reconstruction error ‖x* − x_true‖ / ‖x_true‖ |
| `compute_snr.m` | Computes SNR / ISNR (Improvement in Signal-to-Noise Ratio) |
| `GENERATE_DATA_2D/` | Loads a test image, applies Gaussian blur + AWGN (and/or an inpainting mask) to generate the corrupted input `b` |
| `GENERATE_DATA_1D/` | Loads/synthesizes a 1D test signal, applies blur + AWGN (and/or sub-sampling for inpainting) |

## How to Use

### Prerequisites

- MATLAB R2019a or later
- Image Processing Toolbox (recommended)

### Running the 2D Image Experiments

```matlab
cd src/2D
MAIN_1_EX1
```

Toggle which algorithm to run by editing the flags at the top of the script:

```matlab
TEST_TV_L2_GRAD_DESC    = 1;   % Run TV-L2 restoration via Gradient Descent
TEST_TIK_L2_U_REST      = 0;   % Run TIK-L2 restoration (FFT-based)
TEST_TIK_L2_U_INP       = 0;   % Run TIK-L2 inpainting
```

### Running the 1D Signal Experiments

```matlab
cd src/1D
MAIN_1_EX2
```

Toggle the experiment type:

```matlab
TEST_DENOISING   = 0;
TEST_RESTORATION = 1;   % denoising + deblurring
TEST_INPAINTING  = 0;
```

### Parameter Sweeps

Both main scripts sweep over a range of regularization parameters λ and select the optimum (by relative error / ISNR):

```matlab
lambdas_min = 0.01;
lambdas_max = 0.3;
lambdas_n   = 10;
lambdas     = linspace(lambdas_min, lambdas_max, lambdas_n);
```

## Evaluation Metrics

- **Relative error:** `‖x* − x_true‖ / ‖x_true‖` — lower is better
- **ISNR (Improvement in Signal-to-Noise Ratio):** quantifies the numerical improvement of the restored image/signal over the corrupted input — higher is better
- Both metrics are computed across the λ-sweep, and the optimal λ is the one minimizing relative error / maximizing ISNR

## Author

**Niccolò Chiari**  
Master's degree in Mathematics (Applied Curriculum)  
University of Bologna, Academic Year 2024/2025

Course: *B3093 — Mathematical and Machine Learning Methods in Imaging*

## License

Educational project. Available for academic use.

---

For complete results, figures, and quantitative comparisons, see `results/README.md`. For full derivations and additional examples, consult `docs/REPORT_LAB_1_IMAGING.pdf`. The original assignment specification is available in `docs/Lab_1.pdf`.
