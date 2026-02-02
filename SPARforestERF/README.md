# SPAR-Forest-ERF Instructions

Software is provided in 3 separate **.R** files to
implement the 3 different exposure response functions (ERF) outlined in
the paper, which are listed below.

- **Linear ERF** - g[x(Sₖ)] = αx(Sₖ), which is the simplest specification
and allows one to quantify the overall linear association α between the
exposure and the risk of disease. The function to fit this model is called
SPARForestERFlinear, and is contained in the file **SPARForestERFlinear.R**.
- **Non-linear ERF** - g[x(Sₖ)], which allows the association to be a smooth
non-linear function whose size depends on the level of the exposure.
This smooth function is represented by a Bayesian p-spline, and is
fitted as a second order random walk latent field model using the INLA
software. The function to fit this model is called SPARForestERFnonlinear
and is contained in the file **SPARForestERFnonlinear.R**.
- **Measurement error ERF** - g[x(Sₖ)] = αx(Sₖ) and x(Sₖ) ∼ N(x̃(Sₖ), σ²ₓ),
which uses a Berkson measurement error model to relate the true unknown
exposure x(Sₖ) to the available error-prone estimate x̃(Sₖ). The
error variance σ²ₓ
has to be estimated in advance, and an example is
given in the motivating air pollution study in the paper. The function
to fit this model is called SPARForestERFME and is contained in the
file **SPARForestERFME.R**.

Two additional files are provided with this software bundle:

- **R - illustration of fitting the SPAR-Forest-ERF model.R** - Illustrates
how to fit the SPAR-Forest-ERF model (and a generalized
linear mixed model (GLMM) competitor with linear confounderresponse
associations) with a linear ERF to a data set generated in
the simulation study.
- **Illustrative data for SPAR-Forest-ERF.Rdata** - Contains the simulated
data set used for the above illustration.

The remainder of this read-me file contains details of the arguments and
outputs for all 3 functions. Following the notation in the paper, the data relate
to a set of _K_ areal units, which comprises a response vector of observed
disease counts **_Y_**, an offset vector of expected disease counts **_e_**, a matrix of
_q_ confounders **_Z_**, and a single exposure _x_. Additionally, one needs a _K x K_
neighbourhood matrix **_W_** that specifies the spatial adjacencies between the
_K_ areal units. Then the arguments to the **SPARForestERFlinear**, **SPARForestERFnonlinear**
and **SPARForestERFME** functions are as follows.

- _Y_ - A _K x 1_ vector of observed disease counts.
- _e_ - A _K x 1_ vector of expected disease counts for use as a model offset.
- _x_ - A _K x 1_ vector containing the exposure of interest.
- _Z_ - A _K x q_ data frame of confounders.
- _W_ - A _K x K_ symmetric spatial neighbourhood matrix.
- _epsilon_ - The threshold value for the mean absolute difference in the
ERFs between successive iterations to use for stopping the algorithm.
- _max.iter_ - The maximum number of iterations allowed for the algorithm.
- _ntrees_ - The number of trees to use when tuning the random forest.
- _mtry_ - The set of possible values for the mtry tuning parameter used
in the random forest algorithm.
- _minnode_ - The set of possible values for the minnode tuning parameter
used in the random forest algorithm.
- _n.sample_ - The number of samples to draw from the spatial model to
propagate uncertainty into the random forest.
- _ntrees.sample_ - The number of trees to use for each sampled response
in the random forest.
- _sigma2.x_ - The measurement error variance in the Berkson measurement
error model for the exposure. This argument is for the measurement
error ERF model only.

Once the model has run it returns a _list_ object with the following 4 elements.

- _convergence.summary_ - A 2 column matrix containing the convergence
diagnostics for the algorithm. Column 1 contains the iteration
number and column 2 contains the mean absolute difference in the
estimated ERF from successive iterations.
- _ERF.estimate_ - A _K x 4_ matrix of estimates (column 2), and 95%
credible intervals (columns 3 and 4) for the ERF _{g[x(Sₖ)]}_, as well
as the values of the exposure _x_ (column 1).
- _confounder.estimate_ - A _K x 1_ vector of estimated confounderresponse
associations {m̃[z(Sₖ)]}.
- _mod.spatial_ - The final fitted spatial model produced by the INLA
software. Model fit critera, relative risk estimates and other model
quantities can be extracted from this object.
