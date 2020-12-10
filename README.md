[![License](https://img.shields.io/github/license/RETURN-project/BenchmarkRecovery)](https://choosealicense.com/licenses/apache-2.0/)
[![Build Status](https://travis-ci.org/RETURN-project/BenchmarkRecovery.svg?branch=master)](https://travis-ci.org/RETURN-project/BenchmarkRecovery)
[![codecov](https://codecov.io/gh/RETURN-project/BenchmarkRecovery/graph/badge.svg)](https://codecov.io/gh/RETURN-project/BenchmarkRecovery)
[![codecov](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/)

# Benchmarking recovery indicators derived from remote sensing time series

This project simulates Landsat data and evaluates the performance of recovery indicators with respect to data and disturbance characteristics.

## Background

The context of this project is the study of the recovery of tropical forests after an abrupt disturbance (typically a forest fire) using satellite images as a data source.

The speed of recovery after a disturbance is known to be correlated with the concept of resilience. This is true not only for forests, but for many dynamical systems. To put it simply: forests that recover fast are more resilient. Forests that recover slowly may be in danger of permanent disappearance.

The specialized literature proposes different metrics for measuring the recovery speed. The performance of these metrics depends on many factors. Some of them are natural, such as the intensity of the perturbation or the seasonality. Others are technical, such as the sampling frequency or the spatial resolution.

## Purpose


The purpose of this project is to **efficiently** **compare** the  **reliability** of different post-disturbance recovery **metrics**.

## Mechanics

1. **Infers** time series' **parameters** and characteristics from optical satellite image data​
2. Uses those parameters to **create** a large collection of synthetic (but realistic) **time series**
3. **Calculates** several state-of-the-art recovery **metrics**

## Simplified workflow


![Simplified workflow](./img/flow.png)

## Install

You can install the master version from R via:

```r
library(devtools)
install_github("RETURN-project/BenchmarkRecovery")
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

[Apache](https://choosealicense.com/licenses/apache/)
