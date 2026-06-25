# Results — Image Restoration & Inpainting

This folder contains the figures generated for **Chapter 1 (Image Restoration)** of the report: Tikhonov restoration, Total Variation restoration, and TIK-L2 inpainting, applied to five test images.

## Directory Structure

```
IMAGE_RESTORATION_PICTURES/
├── TIK_L2/                  # Tikhonov (TIK-L2) restoration results
│   ├── Cameraman/std_dvt = 10, 20
│   ├── Parrot/std_dvt = 10, 20
│   ├── Qr code/std_dvt = 10, 20
│   ├── Satellite/std_dvt = 10, 20
│   └── Sinusoid/std_dvt = 10, 200
│
├── TV-RESTORATION/           # Total Variation (TV-L2) restoration results
│   ├── Cameraman/std_dvt = 10, 20
│   ├── Parrot/std_dvt = 10, 20
│   ├── Satellite/std_dvt = 10, 20
│   ├── Sinusoid/std_dvt = 10, 20
│   └── qr_code/std_dvt = 10, 20
│
└── INPAITING/                # TIK-L2 inpainting results
    ├── Cameraman/mask_1, mask_2, mask_3   # Random Point / Horizontal-Strip / Vertical-Strip masks
    ├── Parrot/No noise, std_dvt = 10      # Cage removal (noise-free vs. noisy)
    ├── QR code/mask 1, mask 2, mask 3     # Additional inpainting mask experiments
    └── Squirrel/No noise, std_dvt = 10    # Thread removal (noise-free vs. noisy)
```

> Note: the report's main text discusses results for **σ = 10** (and the Sinusoid/egg-carton example). The `std_dvt = 20` / `std_dvt = 200` subfolders contain additional results at higher noise levels, generated for further comparison but not all individually discussed in the main report text.

## TIK_L2 — Tikhonov Restoration

Images corrupted by Gaussian blur + AWGN, restored by solving the closed-form FFT-based Tikhonov problem `argmin_x ‖Ax-b‖₂² + λ‖Lx‖₂²`.

| Image | Optimal λ | Relative Error |
|---|---|---|
| Cameraman | 0.0467 | 0.10120 |
| QR-code | 0.015 | 0.14162 |
| Sinusoid (egg-carton) | 2.388 | 0.01313 |

Each image folder contains:
- The original image
- The blur kernel and noise realization used
- The corrupted (blurred + noisy) image
- The restored image
- The ISNR-vs-λ plot used to identify the optimal λ

**Best suited for:** smooth images (Sinusoid achieves the lowest relative error by far). Performs worse on images with sharp edges (Cameraman, QR-code), since the quadratic gradient penalty cannot separate true edges from noise.

## TV-RESTORATION — Total Variation Restoration

Same corrupted images as above, restored instead via Gradient Descent on the (smoothed) TV-L2 objective `argmin_x ‖Ax-b‖₂² + λ‖Lx‖₂`.

| Image | Optimal λ | Relative Error | Comparison to TIK-L2 |
|---|---|---|---|
| Cameraman | 0.05 | 0.11190 | **TV better** — preserves sharp edges |
| QR-code | 0.061 | 0.18012 | TIK better — TV underperforms here |
| Sinusoid (egg-carton) | 0.13889 | 0.01884 | TIK slightly better — image is smooth |

**Counter-intuitive QR-code result:** despite being piecewise-constant (in principle ideal for TV), the QR-code's dense fine detail combined with the relatively low noise level (σ=10) means TV-L2's relative error (0.18012) is barely better than the uncorrected corrupted image (0.18346), and worse than TIK-L2 (0.14162).

## INPAITING — TIK-L2 Inpainting

### Object Removal (Squirrel, Parrot)

| Image | Task | Condition | Optimal λ | Relative Error |
|---|---|---|---|---|
| Squirrel | Remove white thread | Noise-free | 8.9×10⁻⁴ | 0.09347 |
| Squirrel | Remove white thread | Noisy (σ=10) | 0.639 | 0.11220 |
| Parrot | Remove cage | Noise-free | 1.03×10⁻² | 0.30233 |
| Parrot | Remove cage | Noisy (σ=10) | 5.57×10⁻¹ | 0.30645 |

Each subfolder (`No noise/`, `std_dvt = 10/`) contains the original image, the mask, the corrupted/masked image, the restored image, and the corresponding ISNR plot.

**Key observation:** optimal λ jumps by 2–3 orders of magnitude between the noise-free and noisy cases — without noise, minimal regularization suffices since the masked-but-otherwise-clean image is already close to the truth; with noise, much stronger smoothing is required.

### Random-Mask Reconstruction (Cameraman)

Three different random inpainting masks applied to the Cameraman image with AWGN (σ=10):

| Mask (`mask_1/2/3`) | Type | Optimal λ | Relative Error |
|---|---|---|---|
| `mask_1` | Random Point Sampling (RPS) | 7.7×10⁻² | 6.9×10⁻² |
| `mask_2` | Random Horizontal Strip Sampling (RHSS) | 8.6×10⁻³ | 6.5×10⁻² |
| `mask_3` | Random Vertical Strip Sampling (RVSS) | 8.4×10⁻³ | 6.9×10⁻² |

**Key observation:** RHSS (`mask_2`) achieves the lowest numerical relative error, but RPS (`mask_1`) gives the best *visual* reconstruction — illustrating that relative error / ISNR don't always align with perceived image quality.

### QR Code Inpainting (`QR code/mask 1, 2, 3`)

Additional inpainting experiments on the QR-code image with different masks, included for completeness alongside the Cameraman mask study above (see the figures in each `mask N/` subfolder).

## Evaluation Metrics Used

- **Relative error:** `‖x_est − x_true‖ / ‖x_true‖` — the primary numerical metric reported throughout
- **ISNR (Improvement in Signal-to-Noise Ratio):** plotted against λ for every experiment to visually identify the optimal regularization parameter; the λ that maximizes ISNR generally coincides with (or is very close to) the λ that minimizes relative error

## Summary — Best Method per Image (Restoration)

| Image | Best Method | Why |
|---|---|---|
| Cameraman | TV-L2 | Sharp edges, TV preserves them better |
| QR-code | TIK-L2 | Dense fine detail + low noise hurts TV here |
| Sinusoid | TIK-L2 (marginally) | Smooth image, ideal case for quadratic penalty |

## References

For the mathematical derivations behind every result in this folder, see:
- Main `README.md` — project overview and full result tables
- `docs/REPORT_LAB_1_IMAGING.pdf` — Chapter 1 (Image Restoration)
- `src/2D/README.md` — code that generated these results
