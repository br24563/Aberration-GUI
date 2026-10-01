# Aberration Explorer GUI

An interactive MATLAB teaching tool for visualizing monochromatic optical
aberrations over circular, annular, or vignetted (off-axis) pupils, using
either a Seidel or a Zernike wave-aberration description.

Launch it, drag the sliders, and watch all of the classic aberration
visualizations update live.

## Features

- **9 classic third-order coefficients** — Seidel W-sums (`W011x`, `W011y`,
  `W111`, `W020`, `W040`, `W131`, `W222`, `W220`, `W311`) or the equivalent
  low-order Zernike circle polynomials (`Z0`–`Z8`), each with a slider **and**
  a numeric-entry box.
- **9 higher-order terms** (4th/5th/6th radial order): trefoil, secondary
  astigmatism, quadrafoil, secondary coma, and secondary spherical aberration
  — toggled on with a checkbox.
- **Non-circular pupils**: annular (central obscuration) or a
  vignetted/decentered aperture (intersection of two offset circles).
- **Field survey mode**: overlays several field heights at once on the
  field-dependent plots (wavefront lineouts, ray intercept, spot diagram, MTF
  lineouts).
- **Polychromatic PSF/MTF**: incoherent sum of monochromatic PSFs over
  wavelengths spanning a user-set fractional bandwidth.
- **All 17 classic visualizations**:
  - Wavefront: 3D surface, 2D color map, x/y lineouts
  - Ray intercept plot and spot diagram
  - MTF (2D map and x/y lineouts) and PSF
  - Single- and double-pass interferograms
  - x/y lateral-shearing interferograms
  - Shack-Hartmann spot pattern
  - Knife-edge and wire shadow tests (horizontal and vertical)

## Requirements

- **Base MATLAB only** — no toolboxes required.
- Tested with MATLAB R2026a (should work on recent releases; the code uses
  standard `uicontrol`/`uipanel` graphics, no App Designer dependencies).

## Usage

```matlab
% Launch the interactive tool
AberrationExplorerGUI

% Run the numerical self-tests (no GUI opens; prints PASS/FAIL)
AberrationExplorerGUI('selftest')
```

The self-test battery checks the gradient formulas against finite
differences, the zero-aberration Strehl ratio, pupil-sampling correctness,
Zernike-fit consistency, and OTF normalization/cutoff. Run it after any edit
to the physics/plotting code.

**Sign / axis convention**: pupil coordinates `(x, y)` are normalized so the
pupil edge sits at radius 1 (`x^2 + y^2 = 1`); the object/image field point
is taken to lie along `+y` at normalized field height `h`.

## Files

| File | Description |
| --- | --- |
| `AberrationExplorerGUI.m` | The complete, self-contained tool (GUI + physics + self-tests) |

## Acknowledgment

The concept for a tool like this — an interactive Seidel/Zernike wavefront
explorer covering wavefront maps, ray-intercept and spot-diagram plots,
MTF/PSF, interferograms, shearing interferograms, Shack-Hartmann spots, and
knife-edge/wire tests — is drawn from a MATLAB teaching tool written by
Prof. Brian Kruschwitz for optics coursework at the University of Rochester
(Institute of Optics).

This repository is a completely new, independent implementation: no source
code, GUI layout, comments, or text from that tool are reused anywhere here,
and it adds a substantial set of features (higher-order aberrations,
non-circular pupils, field surveys, a polychromatic PSF model, and slider
controls). Full credit for the pedagogical concept and the choice of this
particular set of visualizations belongs to Prof. Kruschwitz.

## References

The physics follows standard material covered in any optics textbook:

- V. N. Mahajan, *Optical Imaging and Aberrations*.
- J. W. Goodman, *Introduction to Fourier Optics*.

Specifically: Seidel and Zernike wave-aberration polynomials, and
Fourier-optics PSF/OTF relations via the Wiener–Khinchin theorem.

## License

Licensed under the [MIT License](LICENSE).
