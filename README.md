# m6APrediction

**m6APrediction** predicts **N6-methyladenosine (m6A) sites** from RNA sequence
and genomic features using a bundled random forest classifier trained on 48,055
candidate sites. It is a small, dependency-light package built for the BIO215
module at XJTLU.

m6A is the most abundant internal modification of eukaryotic mRNA and is
deposited at the consensus motif DRACH. Experimental mapping (MeRIP-seq,
miCLIP) is expensive, so a sequence-based predictor lets you prioritise
candidate sites for follow-up validation. The bundled model recognises m6A
sites from seven features — GC content, RNA biotype, transcript region, exon
length, distance to the nearest splice junction, evolutionary conservation and
the five-base sequence context — and reaches an out-of-bag AUC of **0.887**.

## Installation

Install the package directly from GitHub:

```r
# option 1
# install.packages("devtools")
devtools::install_github("keyue24/m6APrediction")

# option 2
# install.packages("remotes")
remotes::install_github("keyue24/m6APrediction")
```

## Quick start

The package ships with a fitted model and a 100-row example feature table, so
every example below runs out of the box.

```r
library(m6APrediction)

# 1. Load the bundled model and example feature table
ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
  package = "m6APrediction", mustWork = TRUE))
example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
  package = "m6APrediction", mustWork = TRUE))

# 2. Batch prediction for a whole feature table
head(prediction_multiple(ml_fit, example_df[1:5, ], positive_threshold = 0.6))
#>   gc_content RNA_type RNA_region exon_length distance_to_junction ...
#> 1  0.6368159     mRNA        CDS    10.76708             8.882643
#>   predicted_m6A_prob predicted_m6A_status
#> 1              0.876             Positive
#> 2              0.000             Negative
#> 3              0.498             Negative

# 3. Single-site prediction
prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA",
  RNA_region = "CDS", exon_length = 12, distance_to_junction = 5,
  evolutionary_conservation = 0.8, DNA_5mer = "GGACA",
  positive_threshold = 0.5)
#> predicted_m6A_prob predicted_m6A_status
#>              "0.51"           "Positive"
```

## Input features

| Feature | Type | Description |
|---|---|---|
| `gc_content` | numeric | GC content of the surrounding region |
| `RNA_type` | factor | `mRNA`, `lincRNA`, `lncRNA`, `pseudogene` |
| `RNA_region` | factor | `CDS`, `intron`, `3'UTR`, `5'UTR` |
| `exon_length` | numeric | exon length, already log2-transformed |
| `distance_to_junction` | numeric | distance to the nearest splice junction, already log2-transformed |
| `evolutionary_conservation` | numeric | evolutionary conservation score |
| `DNA_5mer` | character | five-base context, A at the central position (e.g. `GGACA`) |

`prediction_multiple()` returns the input table with two appended columns:
`predicted_m6A_prob` (probability of the Positive class) and
`predicted_m6A_status` (`Positive` / `Negative`). A site is labelled Positive
only when its probability is **strictly greater** than `positive_threshold`
(default `0.5`).

## Model performance

Performance is reported on the **out-of-bag predictions** of all 48,055
training sites (11,248 positives, prevalence 0.234) — that is, every site is
scored by the trees that did *not* see it during training, so the figures are
not inflated by resubstitution.

| Metric | Value |
|---|---|
| OOB AUC (ROC) | **0.887** |
| OOB AUPRC (PRC) | **0.702** |
| OOB classification error (threshold 0.5) | **0.161** |

![ROC curve](man/figures/roc.png)

![Precision-recall curve](man/figures/prc.png)

The ROC curve rises steeply, reaching ~80% sensitivity at 20% false positive
rate. Because positives are a minority class (23.4%), the precision-recall
curve is the more informative view: precision stays above 0.8 until recall
exceeds 0.4, then falls gradually towards the 0.234 baseline.

## Package structure

```text
m6APrediction/
├── R/
│   └── main.R                    # three functions with roxygen2 documentation
├── inst/
│   └── extdata/
│       ├── rf_fit.rds            # fitted random forest (500 trees, 15 MB)
│       └── m6A_input_example.csv # 100-row example feature table
├── man/
│   ├── dna_encoding.Rd           # internal helper (not exported)
│   ├── prediction_multiple.Rd
│   ├── prediction_single.Rd
│   └── figures/                  # ROC and PRC curves used above
├── DESCRIPTION
├── NAMESPACE
└── README.md
```

## Functions

| Function | Exported | Purpose |
|---|---|---|
| `prediction_multiple(ml_fit, feature_df, positive_threshold = 0.5)` | yes | batch prediction for a feature table |
| `prediction_single(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5)` | yes | prediction for one candidate site |
| `dna_encoding(dna_strings)` | no (internal) | splits 5-mers into five factor-encoded nucleotide positions |

## Licence

MIT.
