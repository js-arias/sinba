# Sinba

Sinba is a R package that implements the sinba (Single birth of a trait) model.

## Installing sinba package from GitHub

To install the sinba package from GitHub,
you probably require a C++ compiler
(such as gcc).
Here is a guide to install gcc in different systems:
<https://gcc.gnu.org/install/>.

Then, use devtools to install sinba:

```{R}
library(devtools)
devtools::install_github(repo = "js-arias/sinba")
```

And that is!

## Basic use of sinba model

The package sinba implements the sinba model.
The main use of the model
is to test the correlation between two traits
using a phylogenetic comparative framework
developed by Pagel (1994).

In the following example,
based on section 7.2.2 of the Revell and Harmon (2022) book,
we first create a random tree and three traits for that data.
Using `phytools` package
we draw the tree and the traits.
This is what has been called "Darwin's dilemma"
and the "unreplicated burst scenario"
(Maddison & FitzJohn, 2015),
which has been a recalcitrant problem for testing correlation.
Then we fit the Pagel model with independent
and correlated traits,
and the same models,
but using sinba.
Under sinba model,
each trait can be born on the tree,
and only after its birth does it behave like a Markov model.
If we compare all the models
under the Akaike information criterion
(Akaike, 1974),
we see that the Pagel model prefers correlation,
whereas sinba prefers the independent model,
and it is the best model of the four examined models.

```{r}
# An example from the book of Revell & Harmon (2022)
# section 7.2.2,
# pp. 173-177.
library(phytools)
library(RColorBrewer)

x <- c(0, 0, rep(1, 12), rep(0, 12))
y <- x
z <- c(0, 0, rep(c(0, 1), 6), rep(0, 12))

plotTree.datamatrix(tree,
  data.frame(
    x = as.factor(setNames(x, LETTERS)),
    y = as.factor(setNames(y, LETTERS)),
    z = as.factor(setNames(z, LETTERS)),
    row.names = LETTERS
  ),
  fsize = 1, space = 0.2, xexp = 1.4,
  palettes = c("YlOrBr", "Greens", "Purples")
)

# --- x and y (Darwin's scenario)

data_xy <- data.frame(LETTERS, x, y)

pagel_xy_ind <- fit_pagel(tree, data_xy, model = new_model("IND"))
pagel_xy_corr <- fit_pagel(tree, data_xy, model = new_model("CORR"))

# Sinba model for the Darwin's scenario with independent evolution
sinba_xy_ind <- fit_sinba(tree, data_xy,
  # as this is the default
  # we can safely ignore both the the pi_x and pi_y parameters
  # we put them here just for completeness
  pi_x = c(0, 1), pi_y = c(0, 1)
)
# Sinba model for the Darwin's scenario with correlated evolution
sinba_xy_corr <- fit_sinba(tree, data_xy,
  model = new_model("CORR"),
  # as this is the default
  # we can safely ignore both the the pi_x and pi_y parameters
  # we put them here just for completeness
  pi_x = c(0, 1), pi_y = c(0, 1)
)

# compare AIC values from Pagel's models
# and Sinba models
aic_xy <- setNames(
  c(
    AIC(pagel_xy_ind),
    AIC(pagel_xy_corr),
    AIC(sinba_xy_ind),
    AIC(sinba_xy_corr)
  ),
  c("Pagel-IND", "Pagel-CORR", "Sinba-IND", "Sinba-CORR")
)
aic_xy

# --- x and z (unreplicated burst scenario)

data_xz <- data.frame(LETTERS, x, z)

pagel_xz_ind <- fit_pagel(tree, data_xz, model = new_model("IND"))
pagel_xz_corr <- fit_pagel(tree, data_xz, model = new_model("CORR"))

# Sinba model for the  unreplicated burst scenario with independent evolution
sinba_xz_ind <- fit_sinba(tree, data_xz,
  # as this is the default
  # we can safely ignore both the the pi_x and pi_y parameters
  # we put them here just for completeness
  pi_x = c(0, 1), pi_y = c(0, 1)
)
# Sinba model for the  unreplicated burst scenario with independent evolution
sinba_xz_corr <- fit_sinba(tree, data_xz,
  model = new_model("CORR"),
  # as this is the default
  # we can safely ignore both the the pi_x and pi_y parameters
  # we put them here just for completeness
  pi_x = c(0, 1), pi_y = c(0, 1)
)

# compare AIC values from Pagel's models
# and Sinba models
aic_xz <- setNames(
  c(
    AIC(pagel_xz_ind),
    AIC(pagel_xz_corr),
    AIC(sinba_xz_ind),
    AIC(sinba_xz_corr)
  ),
  c("Pagel-IND", "Pagel-CORR", "Sinba-IND", "Sinba-CORR")
)
aic_xz
```

See the package vignettes that include a more refined explanation
of this example ("*Basic use of sinba model*"),
a guide to build models ("*How to make custom models*"),
and a tour of most of the functionality of the package
("*A tour of sinba package*").

## Authorship and License

Package sinba is developed by [J. Salvador Arias](http://github.com/js-arias/) and [Sergei Tarasov](https://github.com/sergeitarasov),
with funding from the *Soumen Akatemia* (Finnish Academy of Sciences).
It is an open source project distributed under the [MIT License](LICENSE.md).

## Cited literature

Akaike, H.
(1974).
A new look at the statistical model identification.
*IEEE transactions on automatic control*,
**19**, 716-723.
DOI: [10.1109/TAC.1974.1100705](https://doi.org/10.1109/TAC.1974.1100705).

Maddison, W. P., & FitzJohn, R. G.
(2015).
The unsolved challenge to phylogenetic correlation tests for categorical characters.
*Systematic Biology*,
**64**, 127-136.
DOI: [10.1093/sysbio/syu070](https://doi.org/10.1093/sysbio/syu070).

Pagel, M.
(1994).
Detecting correlated evolution on phylogenies: a general method for the comparative analysis of discrete characters.
*Proceedings of the Royal Society of London. Series B: Biological Sciences*,
**255**, 37-45.
DOI: [10.1098/rspb.1994.0006](https://doi.org/10.1098/rspb.1994.0006).

Revell, L. J., & Harmon, L. J.
(2022).
*Phylogenetic Comparative Methods in R*,
Princeton Univ.,
Princeton (U.S.A.).
