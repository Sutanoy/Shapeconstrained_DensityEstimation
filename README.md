# Shapeconstrained_DensityEstimation
Codes to perform univariate modality constrained density and conditional density estimation
GenerateData is a file that contains several simulated examples which is called by other functions.
FormGammaFromC is the file that constructs \gamma given a coefficient vector in tangent space representation.
FormMModalLikelihoodFromC is the function that computes the objective function at each step of the fmincon algorithm where Boundary values are assumed to be Zero. FormMModalLikelihoodboundaryFromC is the function used when boundaries are NOT constrained to be zero.
ModalConstrainedBoundaryZero is the main file for density estimation when Boundary values are assumed to be zero, ModalConstrainedBoundaryvariable is used when boundary values are not restricted to be zero, BUT have antimodes at the boundaries.
ModalConstrainedCoditional computes the conditional density by calling Modalconstrainedcde function. Here one has a choice to use Meyer basis functions as well.
Meyer basis generator generates the Meyer basis on the [0,1] interval.

We have also added an example R code that mimics the MATLAB version, for k-modal densities with zero boundaries(i.e. f(0)=f(1)=0), support on unit interval,  and an ad-hoc chosen number of basis elements. The basis elements can be chosen via a penalized version of the likelihood function, like AIC. The support can be generalized by scaling the observed data points to the unit interval as a pre-processing step.
