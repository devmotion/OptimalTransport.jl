# OptimalTransport.jl

Computation of optimal transport maps.

[![Build Status](https://travis-ci.com/devmotion/OptimalTransport.jl.svg?branch=master)](https://travis-ci.com/devmotion/OptimalTransport.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/devmotion/OptimalTransport.jl?svg=true)](https://ci.appveyor.com/project/devmotion/OptimalTransport-jl)
[![Codecov](https://codecov.io/gh/devmotion/OptimalTransport.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/devmotion/OptimalTransport.jl)
[![Coveralls](https://coveralls.io/repos/github/devmotion/OptimalTransport.jl/badge.svg?branch=master)](https://coveralls.io/github/devmotion/OptimalTransport.jl?branch=master)

## Overview

This package allows the computation of optimal transport maps in the Julia language.

## Usage

Currently only the network simplex method is implemented. An optimal transport map
`P` and the dual potentials `f` and `g` can be computed with
```julia
P, f, g = earthmover(a, b, C)
```
where `a` and `b` are two histograms and `C` is the cost matrix.

## Bibliography

Peyré, G., & Cuturi, M.. (2018). Computational Optimal Transport. [arXiv:1803.00567](https://arxiv.org/abs/1803.00567).