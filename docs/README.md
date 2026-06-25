# Documentation

## REPORT_LAB_1_IMAGING.pdf

The complete technical report, covering two chapters:

### Chapter 1 — Image Restoration (pages 1–20)
1. **1.1 Tikhonov Restoration (TIK-L2)** — closed-form FFT-based solution; results on Cameraman, QR-code, and a sinusoidal/egg-carton image
2. **1.2 Total Variation Restoration (TV-L2)** — Gradient Descent solution; same three test images, compared against TIK-L2
3. **1.3 Image Inpainting** — TIK-L2 formulation with a binary mask S; object-removal examples (squirrel, parrot) and random-mask reconstruction (Cameraman, 3 mask types)

### Chapter 2 — Signal Restoration (pages 21–32)
1. **2.1 Signal Denoising** — Generalized Tikhonov (2nd-order derivative penalty) on 4 test signals
2. **2.2 Signal Denoising and Deblurring** — same signals, now also blurred
3. **2.3 Signal Inpainting** — sub-sampling-matrix formulation for recovering missing signal samples

## Lab_1.pdf

The original assignment specification (course: *B3093 — Mathematical and Machine Learning Methods in Imaging*, University of Bologna, A.Y. 2023-2024), describing the two exercises this project implements:

1. **2D Image processing** — complete `MAIN_1_EX1.m` implementing TIK-L2 restoration (FFT-based, no explicit matrix A), TV-L2 restoration (gradient descent), and TIK-L2 inpainting
2. **1D Signal processing** — complete `MAIN_1_EX2.m` implementing denoising, deblurring, and inpainting for 1D signals

## Reading Guide

| If you want to understand... | Read... |
|---|---|
| The TIK-L2 vs. TV-L2 comparison on images | Report Ch. 1.1–1.2 |
| Why TV underperforms on the QR-code | Report Ch. 1.2 (QR-code subsection) |
| The inpainting mask formulation | Report Ch. 1.3 (intro) |
| Why λ changes so much between noisy/noise-free inpainting | Report Ch. 1.3 (Squirrel/Parrot examples) |
| The 1D signal inpainting sub-sampling-matrix trick | Report Ch. 2.3 |
| What exactly was assigned vs. what was explored further | `Lab_1.pdf` vs. `REPORT_LAB_1_IMAGING.pdf` |

## Code ↔ Report Mapping

| Report Section | Code |
|---|---|
| 1.1 TIK-L2 restoration | `src/2D/MAIN_1_EX1.m` (`TEST_TIK_L2_U_REST` branch) |
| 1.2 TV-L2 restoration | `src/2D/MAIN_1_EX1.m` (`TEST_TV_L2_GRAD_DESC` branch) + `src/2D/TV_L2_U_GD.m` |
| 1.3 Image inpainting | `src/2D/MAIN_1_EX1.m` (`TEST_TIK_L2_U_INP` branch) |
| 2.1–2.3 Signal experiments | `src/1D/MAIN_1_EX2.m` (`TEST_DENOISING`, `TEST_RESTORATION`, `TEST_INPAINTING` branches) |

## Author

**Niccolò Chiari**, Master's degree in Mathematics (Applied Curriculum), University of Bologna, A.Y. 2024/2025.
