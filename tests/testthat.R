#!/usr/bin/env Rscript

library(devtools)
#devtools::load_all()
library(testthat)
library(mrplew)

test_check("mrplew")
#test_file("testthat/test_mrplew.R")
