# Logistic regression with misclassified binary responses: a simple majority approach based on repeated classifications

R programs for the paper

> V. B. Quinino, D. G. O. Valadares and R. C. Quinino. Logistic regression with
> misclassified responses: a simple majority approach based on repeated
> classifications. Manuscript submitted for publication.

and its Supplementary Material. If you use these programs, please cite the
paper.

## Programs

Each program is named after the section it reproduces; open it and run it (in
RStudio, **Source**; or `Rscript <program>.R`). Keep all files in the same
folder. Results are written to `results/` and figures to `figures/`.

| Program | Reproduces |
|---|---|
| `Section4_2_analytical.R` | Figure 1, Table 1, Figures S1 and S2 |
| `Section4_3_finite_sample.R` | Table 2, Tables S1 to S3, Figure S3 |
| `Section4_4_number_of_classifications.R` | Table 3, Table S4 |
| `Section4_5_unanimous.R` | Section 4.5, Table S5 |
| `Section5_lens.R` | Tables 6 and 7, Table S7, and the goodness of fit of Section 5 |
| `Supp_S2_binomial_regression.R` | Section S2 |
| `Supp_S3_e_estimators.R` | Section S3, Table S6 (estimators of e1 and e2) |
| `Supp_S4_fla_example.R` | Section S4 (example of `fla_fit`) |
| `sm_fit.R` | the SM, the estimators of e1 and e2 and the diagnostic, for use with other data (instructions at the top of the file) |
| `fla_fit.R` | the FLA, for use with other data (instructions at the top of the file) |
| `functions.R` | functions used by the programs; it does not need to be run |

`Section4_5_unanimous.R` and `Supp_S2_binomial_regression.R` use the
replicates saved by `Section4_3_finite_sample.R`. These replicates are
included in `results/`, so the programs run in a few minutes; to simulate
again (about 40 minutes), set `rerun <- TRUE` in `Section4_3_finite_sample.R`.
`Supp_S3_e_estimators.R` also reads its saved replicates from `results/`;
set `rerun <- TRUE` there to simulate again (about 20 minutes).

## Data

The lens inspection data of Section 5 (Table 5 of the paper: number of lenses
by number of 'perfect' classifications among five, for each combination of
laboratory and region) are entered at the beginning of `Section5_lens.R`.

## Requirements

The results were obtained with R 4.3.3 and are the same on Windows, Linux and
macOS. The package `lpSolve` is required; everything else uses base R. EPS
copies of the figures are created with `pdftops` (Poppler) when it is
available; otherwise only the PDF figures are written.

## License

The programs are distributed under the MIT License (see `LICENSE`).
