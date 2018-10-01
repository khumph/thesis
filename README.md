# Using reinforcement learning to personalize dosing strategies in a simulated cancer trial with high dimensional data

This repository contains the code and files to generate the write-ups for my master's thesis "Using reinforcement learning to personalize dosing strategies in a simulated cancer trial with high dimensional data".

## Abstract

In a simulation of an advanced stage generic cancer trial, I use Q-learning, a reinforcement learning algorithm, to develop dynamic treatment regimes for a continuous treatment, the dose of a single drug. Selected dynamic treatment regimes are tailored to time-varying patient characteristics and to patient subgroups with differential treatment effects. This approach allows estimation of optimal dynamic treatment regimes without a model of the disease process or a priori hypotheses about subgroup membership. Using observed patient characteristics and outcomes from the simulated trial, I estimate Q-functions based on 1) a single regression tree grown by the Classification And Regression Trees (CART) method, 2) random forests, and 3) a slightly modified version of Multivariate Adaptive Regression Splines (MARS). I then compare the survival times of an independent group of simulated patients under treatment regimes estimated using Q-learning with each of the three methods, 10 constant dose regimes, and the best possible treatment regime chosen using a brute force search over all possible treatment regimes with complete knowledge of disease processes and their effects on survival. I also make these comparisons in scenarios with and without spurious high dimensional covariates and with and without patient subgroups with differential treatment effects. Treatment regimes estimated using Q-learning with MARS and random forests greatly increased survival times when compared to the constant dose regimes, but were still considerably lower than the best possible dose regime. Q-learning with a single regression tree did not outperform the constant dose regimes. These results hold across high dimensional and subgroup scenarios. While the MARS method employed produces much more interpretable models than random forests, and therefore has more promise for patient subgroup identification, I show that it is also more sensitive to variations in training data.

## Reproducing the analyses

### Dependencies

The project was tested on macOS 10.14 Mojave using the following:

- [R](https://www.r-project.org/) version 3.5.1
- The R packages [pacman](http://trinker.github.io/pacman_dev/), [data.table](https://github.com/Rdatatable/data.table/wiki), [tidyverse](https://www.tidyverse.org/), [docopt](https://github.com/docopt/docopt.R), [caret](http://topepo.github.io/caret/index.html), rpart, earth, [ranger](https://github.com/imbs-hl/ranger), and [tictoc](http://collectivemedia.github.io/tictoc/).
- [GNU Make](https://www.gnu.org/software/make/) 3.81 - standard on macOS and Linux
- [GNU Bash](https://www.gnu.org/software/bash/) 3.2.57(1) - default shell in terminal on macOS and many Linux distributions

You can ensure you have the required R packages by running the following in an R console:

```r
# Install pacman, if needed
if (!suppressWarnings(require("pacman", quietly = TRUE))) {
   install.packages("pacman", repo = "https://cran.rstudio.com")
}

# Load docopt, install if not installed
pacman::p_load(docopt)
```

The magic of the `pacman` R package will then install any other required packages as they are needed.

### (Re)making the analyses

First, download the project files (e.g., by hitting the "Clone or download" button above).

Then navigate to the main directory of the project (called thesis by default) in a terminal window. Typing `make help` when in the main project directory displays a short description of all of the options to `make`. For example, the previous two steps might look like:

```bash
cd ~/Downloads/thesis
make help
```

To simulate the training clinical trail data, for instance, type

```bash
make simulate
```

when in the `thesis` directory.

Running `make all` (or simply `make`) will produce the following:

1. A directory called "data" with the simulated data files, and baseline data for testing the models.
2. A directory called "results" with R files for the models, and simulated data obtained by using the candidate treatment regimes.
3. A directory called "cache" with cached results from the .Rnw file.
4. A directory called "figure" with figures produced from .Rnw file for inclusion in the write-up.

To speed up reproducing the analyses, you may want to fit a smaller number of repeated samples to train from and fewer test observations for each regime. This reproduces the main thrust of the results without taking much time and disk space. To change these parameters, alter the final two lines in `config.mk`. For example you could change them to:

```make
N_SAMPLES = 5
N_TEST_OBS = 500
```