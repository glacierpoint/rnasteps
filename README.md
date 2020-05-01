rnasteps
========

Fit steps to single molecule traces such as RNA unwinding traces

Single-molecule real time trajectories are embedded in high noise. To extract kinetic or dynamic information of the molecules from these trajectories often requires idealization of the data in steps and dwells. One major premise behind the existing single-molecule data analysis algorithms is the Gaussian ‘white’ noise, which displays no correlation in time and whose amplitude is independent on data sampling frequency. This so-called ‘white’ noise is widely assumed but its validity has not been critically evaluated. This R package quantitatively evaluates the structure of the underlying noise, takes the noise structure into account, and identifies steps and dwells in a single-molecule trajectory. Unlike existing data analysis algorithms, this method uses Generalized Least Squares (GLS) to detect steps and dwells. Under the GLS framework, the optimal number of steps is chosen using model selection criteria such as Bayesian Information Criterion (BIC). Comparison with existing step detection algorithms showed that this GLS method can detect step locations with highest accuracy in the presence of correlated noise. Because this method is automated, and directly works with high bandwidth data without pre-filtering or assumption of Gaussian noise, it may be broadly useful for analysis of single-molecule real time trajectories.

Load the package into R as follows:

library(devtools)

install_github("glacierpoint/rnasteps")
library(rnasteps)
library(ifultools)

Reference:

Arunajadai SG, Cheng W (2013) Step Detection in Single-Molecule Real Time Trajectories Embedded in Correlated Noise. PLoS ONE 8(3): e59279. doi:10.1371/journal.pone.0059279
