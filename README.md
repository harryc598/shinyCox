
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rshinycox

The goal of rshinycox is to provide a tool for researchers to easily
create shiny applications to display one or more Cox model-predicted
survival curves. This work was motivated by clinical research to provide
an interactive summary of fitted Cox models, make it easy for
researchers to generate these summaries for their own research, and to
support dissemination by minimizing the use of subject-identifiable
data.

Tables of hazard ratios remain an important summary of fitted Cox
models, but predicted survival curves can provide a different summary
that may be more interpretable to clinicians and patients.

While shiny apps are excellent tools for visualization, development of
these apps can be challenging for unfamiliar users. This package allows
a user to input their final Cox model(s) into the function
`shine_coxph()` to create a shiny app for the user with minimal
additional work required.

The function `shine_coxph()` initially uses the original data for
internal verification purposes, but discards the original data when the
app is made, containing only the bare necessities to make predictions
(predictions are made in the application using the object returned by
`prep_coxfit()`). When disseminating the final app, subject-identifiable
data are not distributed.

## Installation

You can install this package using:

``` r
install.packages("rshinycox")
```

You can install the development version of rshinycox from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("harryc598/rshinycox")
```

## Example

Here is an example using multiple treatment arms, which will highlight
the usefulness of this package. If you move around the inputs, you can
see how for different combinations, one treatment might prevail over the
other.

``` r
library(rshinycox)
library(survival)
split_colon <- split(colon, colon$rx)

colon_arm1 <- split_colon$Obs
colon_arm2 <- split_colon$Lev
colon_arm3 <- split_colon$`Lev+5FU`

colon1ph <- coxph(Surv(time, status) ~sex +  age + obstruct + nodes, colon_arm1,
                  x = TRUE, model = TRUE)
colon2ph <- coxph(Surv(time, status) ~ sex + age + obstruct + nodes, colon_arm2,
                  x = TRUE, model = TRUE)
colon3ph <- coxph(Surv(time, status) ~ sex + age + obstruct + nodes, colon_arm3,
                  x = TRUE, model = TRUE)
# This will write a shiny app into whatever your working directory is
# shine_coxph("arm1" = colon1ph, "arm2" = colon2ph, "arm3" = colon3ph, theme = "dashboard", app.dir = getwd())
```
