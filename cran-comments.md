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

## Changes in this version (1.3.4)

* Updated test cases for compatibility with the forthcoming rstan 2.39.0.
  Changes to the Stan RNG engine caused some MLE tests to fail to initialise
  from random starting values; those tests now initialise at 0, and one
  snapshot comparison was replaced with a tolerance-based expectation.
  Tests pass under both the current CRAN rstan and rstan 2.39.0.
  Contributed by Andrew Johnson (@andrjohns).

This is a test-only change; no user-facing code was modified.
