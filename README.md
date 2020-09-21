# emcite
Code for Batch Mode Active Learning for Individual Treatment Effect Estimation - Zoltan Puha, Maurits Kaptein and Aurelie Lemmens, IncrLearn @ ICDM 2020.

This repo contains code to run the B-EMCMITE algorithm presented in the paper. 
There are two main functions: 
- `dgp` creates different datasets
- `af` runs the acquisition function.

A full Active Learning procedure can be seen in `vignettes/emcite.Rmd` (with all the helping functions)
