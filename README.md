Diagnostic and Learning Tools for Time Series
================
Joshua Eagan
December 12, 2020

## Project Summary

The R functions in this directory will be useful in teaching students of
time series analysis who are beginning applied work with ARMA(p,q)
models and for more seasoned time series analysts looking for improved
diagnostic tools.

[RPubs project write up](https://rpubs.com/Jmemq7/705674)

-----

## Functions for Learning Time Series

The `ar1.acf()` and `ma1.acf()`functions are designed to help students
learn how to use a graph of the autocorrelation function of a timeseries
to decide whether MA(1) or AR(1) models are the right choice to use in
forcasting future values.

### ar1.acf

This function takes a time series and computes confidence intervals
around its ACF on the assumption that its autocorrelation suggests the
AR(1) model is appropriate. It uses variances from Bartlettâs
Approximation, shown below.

![Alt
text](C:/Users/Josh/Desktop/TS%20analysis/final%20project/Diagnostic%20and%20Learning%20Tools%20for%20Time%20Series/equasions/1.png)

For the AR(1) process, this simplifies down to:

![Alt
text](C:/Users/Josh/Desktop/TS%20analysis/final%20project/Diagnostic%20and%20Learning%20Tools%20for%20Time%20Series/equasions/2.png)
for = 1,2, â¦ where phi is the sample autocorrelation at lag tau

The actual confidence intervals are calculated using the following
formula:

![Alt
text](C:/Users/Josh/Desktop/TS%20analysis/final%20project/Diagnostic%20and%20Learning%20Tools%20for%20Time%20Series/equasions/3.png)

The confidence intervals can be compared against theoretical
autocorrelation values of an AR(1) parameter provided by the user.
`ar1.acf()` uses `ARMAacf()` to retreive the theoretical
autocorrelations.

The `ar1.acf()` function takes the following inputs:

  - **x** time series (object with the class `ts`) to be tested against
    theoretical autocorrelation function values. Confidence intervals
    will be estimated around the estimated autocorrelation function of
    this series.
  - **true\_ar** autoregressive parameter to generate a theoretical
    autocorrelation function to be compared against the provided time
    series
  - **lag.max=NULL** maximum number of lags to estimate autocorrelation
    (defaults to a number of lags chosen by the ACF function based on
    the length of the series)
  - **pvalue=.05** p-value desired for the confidence intervals around
    the autocorrelation function of **x**

<!-- end list -->

``` r
#loading in the required TSA package
library(TSA)

#simulating an AR(1) time series
set.seed(1406565)
true_ar<-.8
x <- arima.sim(n=100,model=list(ar=true_ar))

#run
ar1.acf(x, true_ar = true_ar, lag.max=10)
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

The dots in the graph are the values from the true autocorrelation of an
AR(1) series with the parameter p=.8. The confidence intervals are fit
around the true autocorrelation of x using u where alpha (the pvalue)
was .05.

Notes:

  - The function will not work correctly unless the `TSA` package is
    loaded.
  - The `pvalue` argument automatically updates the title.
  - You can use this function on data for which the true ar(1) parameter
    is unknown. The true autocorrelation can only be known in the case
    of simulated data, which was done in this case to show how the
    function works.

### ma1.acf

Much like `ar1.acf`, this function takes a time series, finds its ACF,
and computes confidence intervals around the ACF testing the assumption
that the data is an MA(1) series using variances from Bartlettâs
Approximation (shown above.)

For the MA(1) case, Bartlettâs Approximation simplifies to the
following:

![Alt
text](C:/Users/Josh/Desktop/TS%20analysis/final%20project/Diagnostic%20and%20Learning%20Tools%20for%20Time%20Series/equasions/4.png)

As above, the confidence intervals around the sample ACF are computed
using the formula below.

![Alt
text](C:/Users/Josh/Desktop/TS%20analysis/final%20project/Diagnostic%20and%20Learning%20Tools%20for%20Time%20Series/equasions/5.png)

The `ma1.acf()` function takes the following inputs:

  - **y** time series (object with the class `ts`) to be tested against
    theoretical autocorrelation function values. Confidence intervals
    will be estimated around the estimated autocorrelation function of
    this series.
  - **true\_ma** MA parameter to generate a theoretical autocorrelation
    function for comparison against the ACF of **y**
  - **lag.max=NULL** maximum number of lags to estimate autocorrelation
    (defaults to a number of lags chosen by the ACF function based on
    the length of the series)
  - **pvalue=.05** p-value desired for the confidence intervals around
    the autocorrelation function of **y**

<!-- end list -->

``` r
#loading the required TSA package
library(TSA)

#simulating an MA(1) time series
set.seed(240)
true_ma<-.8
y <- arima.sim(n=100,model=list(ma=true_ma))

#run
ma1.acf(y, true_ma=true_ma)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

The dots in the graph are the values from the true autocorrelation of an
MA(1) series with the parameter q=.8. The confidence intervals are fit
around the true autocorrelation of y using u where alpha (the pvalue)
was .05.

Notes:

  - The function will not work correctly unless the `TSA` package is
    loaded.
  - The `pvalue` argument automatically updates the title.
  - You can use this function on data for which the true ma(1) parameter
    is unknown. The true autocorrelation can only be known in the case
    of simulated data, which was done in this case to show how the
    function works.

-----

## Diagnostic Tools

`eacf_pic()` and `ar1.est()` are a pair of functions which will be more
useful in an applied setting where you donât have the ability to know
the true ARMA(p,q) order underlying a timeseries.

### eacf\_pic

The `eacf_pic()` function produces a visualization of the sample
extended autocorrelation function using
[geom\_tile()](https://ggplot2.tidyverse.org/reference/geom_tile.html)
in the `ggplot2` package. The main value added of this in comparison to
the
[eacf()](https://www.rdocumentation.org/packages/TSA/versions/1.3/topics/eacf)
function in the `TSA` package is that rather than a binary indicator
significance, âXâ or âOâ, being displayed in the console, this function
produces a picture of the eacf that clearly shows significance using
colors (red tiles are significant at the 5% level and grey tiles are
insignificant) while also indicating the level of significance using
color intensity.

`eacf_pic()` function takes the following inputs:

  - **x** time series object for which to produce the eacf
  - **ar.max=7** maximum number of autoregressive coefficients to test
  - **ma.max = 13** maximum number of moving average coefficients to
    test
  - **text\_pvals=F** an option to add p-values rounded to the nearest
    hundreth to the tiles

<!-- end list -->

``` r
#loading required libraries
library(TSA)
library(reshape2)
library(tidyverse)

#simulating some time series data
set.seed(243)
x=arima.sim(n = 500, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), sd = sqrt(0.1796))

#run
eacf_pic(x)
```

    ## No id variables; using all as measure variables

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#the output of the eacf function from the TSA package
eacf(x)
```

    ## AR/MA
    ##   0 1 2 3 4 5 6 7 8 9 10 11 12 13
    ## 0 x x x x x o x o o o o  o  o  o 
    ## 1 x x x x x o x o o o o  o  o  o 
    ## 2 x x o o o o o o o o o  o  o  o 
    ## 3 x x x o o o o o o o o  o  o  o 
    ## 4 x x x x o o o o o o o  o  o  o 
    ## 5 x x x x x o o o o o o  o  o  o 
    ## 6 x x x o o o o o o o o  o  o  o 
    ## 7 x x x o o x o o o o o  o  o  o

Notes:

  - The function will not work correctly unless the `TSA`, `reshape2`,
    and `tidyverse` packages are loaded.
  - You will always receive the message `No id variables; using all as
    measure variables`. This is not an error, it is automatically
    printed by the `melt()` function in `reshape2`.
  - A careful reader might notice that there is no color to represent
    the p-value threshhold \(.10 < P </= .20\). This is because there
    are no such values in our example graph, so it was ommitted in the
    legend.
  - The text p-values might be helpful to some, but they are generally
    undesirable due to the neccecity to round to the nearest hundreth.
    This makes p-values below .005 appear as 0.

### AR(1) Diagnostic tool:

Inspired by a [Stack Exchange
post](https://stats.stackexchange.com/questions/310470/lag-wise-confidence-band-for-sample-autocorrelation-function-of-ar1-process#310578),
this function automates an answer given to the problem of computing
lag-wise confidence bands for the sample autocorrelation function of an
AR(1) process. This function is potentially a better solution than
`ar1_acf` in an applied setting where the time series is an AR(1)
process but the true autocorrelation is not known.

The `ar1.est()` function takes only 1 input:

  - **x** a time series suspected to be an AR(1) series.

<!-- end list -->

``` r
#loading required package
library(TSA)

#simulating an AR(1) process
set.seed(25)
x <- arima.sim(n=100,model=list(ar=-.8))

#run
ar1.est(x)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Credit

Each of these functions are built using functions from the
[TSA](https://CRAN.R-project.org/package=TSA) package as well as
functions included in base R. The `eacf_pic()` function relies on a
function in
[ggplot2](https://www.rdocumentation.org/packages/ggplot2/versions/3.3.2)
as well as one in
[reshape2](https://www.rdocumentation.org/packages/reshape2/versions/1.4.4).
I leaned on [this
post](https://stats.stackexchange.com/questions/310470/lag-wise-confidence-band-for-sample-autocorrelation-function-of-ar1-process#310578)
for the creation of `ar1.est` and for some of the syntax for plotting
the `ar1.acf` and `ma1.acf` functions. Special thanks to [Dr.Â Lori
Thombs](http://faculty.missouri.edu/~thombsl/) for her guidance on this
project and for her ideas to write the first three functions in this
document.
