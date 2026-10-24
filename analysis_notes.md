# Lab 5 analysis and review notes

The source notebook is main.tex, preserved unchanged. report.tex is the report entry point. The supplied handout is LabON05-AtomicSpectra-v25.1.0.pdf. Its internal headings say Lab 6, but the user's folder and assignment are Lab 5.

## Experiment 2

The helium D3 offsets are 11.11, 11.18, 11.30, and 11.30 cm. Their unweighted mean is 11.2225 cm, sample SD is 0.09394147114 cm, and standard error is 0.04697073557 cm. This follows the handout's explicit instruction to use SD for a new spectral line position and standard error for the mean D3 position. The original 0.1 cm hydrogen estimates are consistent with that empirical scale. The violet reading had no separately recorded uncertainty.

Use the handout's D3 wavelength of 587.6 nm. The helium grating to scale separation is 45.71 - 9.13 = 36.58 cm. The calculated spacing is about 2003.405117 nm, with uncertainty about 7.663739 nm from the mean D3 position. This is a calculation from a known reference, not a fit.

The hydrogen separation changed to 42.20 - 14.35 = 27.85 cm. Reusing the helium separation would give incorrect wavelengths. Use lambda = d y / sqrt(L² + y²), retaining the full geometry. Propagate the D3 mean uncertainty through d and the helium sample SD through each hydrogen offset. The bench position uncertainties were unrecorded. Both separations are held fixed, so quoted uncertainties exclude that contribution and any common alignment error.

Assign the red, cyan, and blue/violet lines to initial levels 3, 4, and 5, with final level 2. The wavelengths are about 663.790494, 486.389685, and 438.017843 nm, with propagated uncertainties about 6.218146, 6.443619, and 6.498741 nm. Evaluate the predictions from the handout's Rydberg constant, 1.0973731568508e7 per m, before rounding. They are 656.112276, 486.009094, and 433.936691 nm. The individual deviations are 1.23, 0.059, and 0.628 quoted uncertainties. Position placement contributes about 6 nm per line, compared with about 2 nm from the grating spacing. The recorded line width is comparable to the position scatter and is not counted again as a separate independent error.

Fit reciprocal wavelength against 1/ni² using generalized least squares. Both slope and intercept are free. The reciprocal wavelength covariance retains the common grating error. The Rydberg constant is minus the slope. Recover nf = sqrt(-slope/intercept), retaining slope and intercept covariance. The fit gives R about 1.106853077e7 per m, uncertainty about 4.124411775e5 per m, and nf about 2.010925030 ± 0.023333433. Do not shrink these input error scales based on the small residual. Three points give one degree of freedom, with chi² about 0.2904.

The conditional fit p value is 0.589949. A separate comparison with the fixed Balmer prediction uses the full wavelength covariance and gives chi² = 1.789639 for three degrees of freedom, with conditional p = 0.617192. The R and nf estimates differ from their references by 0.229849 and 0.468214 quoted uncertainties. These p values treat the estimated reading scale as a fixed Gaussian standard uncertainty. With four reference placements and one free-fit residual degree of freedom, the tests provide approximate consistency checks, not a precise validation of the error scale. A common multiplicative wavelength error changes R but cancels in the slope-to-intercept ratio that determines nf.

The report gives the added setup2.jpg photograph, geometry sketch, reference spectrum, and fit separate full-width spaces. All report photos are now vertically arranged. Photographs are cropped only within the document layout, preserving the originals. Residual bars are the same propagated input uncertainties as the main fit and the points are not connected. The spectrum does not claim measured line intensities. The records do not specify which side was used for each hydrogen offset, so the geometry sketch labels a generic apparent displacement rather than inventing a second set of measurements.

No hydrogen reflection grating trace or UV fluorescence observation was recorded. The optional helium Bohr fit in Section 2.3 is not required and is omitted.

## Earlier report feedback applied

Reviewed lab0/lab0GRADEDsubmission.pdf, lab1/graded5clLab1.pdf, and the graded 5BL report Example/submission_404800050.pdf. Also read the 5BL capstone guidelines for general analysis expectations without importing its report structure.

The Lab 0 and Lab 1 selected deductions mostly concern setup images, sketches, and measurement uncertainties. Lab 0 also notes that residual plots need error bars. Lab 1's qualitative comparison deduction reinforces comparing like quantities rather than making unsupported claims. These are the selected deductions, not a claim that an analysis penalty was selected.

The 5BL grader asks for important results in the introduction, numbered equations, and discrete residual points, and queries whether residuals were normalized. Its analysis and the full-mark Lab 0 section give examples of propagated errors, fit statistics, and explicit agreement checks.

Experiment 2 now includes numerical results in the introduction, an apparatus geometry sketch, and the hydrogen photograph. The added analysis explains the sample SD versus standard error choice, derives the grating and wavelength uncertainties, and compares calculated wavelengths with Balmer predictions. It reports the free fit's slope and intercept uncertainties, retains their covariance in nf, and interprets chi² and agreement tests. A concise discussion identifies the dominant measured uncertainty and the unrecorded systematic contributions. Equation numbers and raw residual error bars address the presentation feedback.

## Experiment 3

The user confirmed a modified qualitative lab following flooding on the first workday. No quantitative analysis of mystery gas readings is included. The notes record warm colored features for B and a broad range of visible lines for C. Within the handout's candidate set, these favor neon and argon. The notes do not resolve nitrogen's molecular band structure, so C's identification remains less secure. No unobserved absence of violet, exact line count, or measured intensity ratio is claimed.

Reference tables inspected:

- NIST neon: https://physics.nist.gov/PhysRefData/Handbook/Tables/neontable2.htm
- NIST argon: https://physics.nist.gov/PhysRefData/Handbook/Tables/argontable2.htm
- NIST hydrogen: https://physics.nist.gov/PhysRefData/Handbook/Tables/hydrogentable2.htm
- NIST helium: https://physics.nist.gov/PhysRefData/Handbook/Tables/heliumtable2.htm

## Critical issues in Isaac's draft

The helium wavelength column contained 0.1 and 0.15 without wavelength units. These were placeholders, not calculated visible wavelengths. They are replaced with values from the exact geometry and the D3 reference. The helium wavelength table now includes uncertainties propagated by the same method as Experiment 2, about 5 nm for each independently measured additional line. D3 remains the fixed reference. The repeated-measurement table also swapped Isaac's left and right readings relative to main.tex. The swap is corrected and does not change the mean or SD.

The incandescent continuum was attributed to pure white light and the eye's sensitivity. Its physical origin is thermal emission from the hot filament. The revised wording keeps Isaac's phrasing with a brief thermal explanation. The CFL description keeps the observed bright features and overall smearing, with a short statement that they are consistent with a mercury vapor CFL. The detailed phosphor explanation was removed at the user's request.

The copied report title still named Lab 4 and the photoelectric effect. It now names Lab 5 and atomic spectra. Placeholder captions and duplicate labels are corrected. Other minor language cleanup supports presentation while preserving the recorded observations.
