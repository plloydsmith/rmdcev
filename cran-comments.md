## Test environments
* Windows 11 (local), R 4.1+
* macOS (on github-actions), R release
* ubuntu (on github-actions), R release and devel

## R CMD check results

0 errors | 0 warnings | 2 notes

❯ checking installed package size ... NOTE
    installed size is  5.1Mb
    sub-directories of 1Mb or more:
      libs   4.3Mb

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

## Changes in this version (1.3.0)

* Moved primary estimation backend from rstan to cmdstanr. The rstan backend
  is retained as an optional alternative (rstan is in Suggests).
* Removed hard NAMESPACE imports from rstan; all rstan usage is now conditional
  on the user requesting backend = "rstan" and rstan being installed.
* Added rng_utils.cpp to supply rmdcev_get_rng() / rmdcev_get_stream(),
  replacing the previous runtime dependency on rstan for simulations.
* Fixed summary(), PrepareSimulationData(), and print() methods to work with
  both the cmdstanr and rstan backends.
* STANC_FLAGS in Makevars.win no longer queries rstan at compile time.
