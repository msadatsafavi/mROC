Tutorial for model-based ROC (mROC) analysis
================

This document provides background information and stepwise tutorial for
the R package on model-based ROC (mROC) methodology. This is related to
the manuscript “Model-based ROC (mROC) curve: examining the effect of
case-mix and model calibration on the ROC plot”, currently available
from <https://arxiv.org/abs/2003.00316>

Authors: Mohsen Sadatsafavi, Paramita Saha-Chaudhuri, John Petkau

For inquiries and reporting bugs, please contact Mohsen Sadatsafavi
(<msafavi@mail.ubc.ca>)

## What is mROC?

Imagine you have developed a risk prediction model using some
development dataset. The risk prediction model takes in some predictors
(e.g., sex, age, previous disease history) and returns the risk of an
event (e.g., risk of disease relapse in the next 12 months). You would
like to evaluate the performance of the risk model in a new (external)
validation sample. Among other things, you typically evaluate the
Receiver Operating Characteristic (ROC) curve of the risk prediciton
model in the new sample.

Now, ‘model-based’ ROC (mROC) curve is the ROC curve that should be
observed if the prediction model is calibrated in the external
population. Comparing the empirical ROC and mROC curves in the new
sample can be informative on if the model is calibrated in the new
sample.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("msadatsafavi/mROC")
```

## How the package works

The package provides simple functions for mROC-related methods. It also
comes with exemplary datasets. Below we provide a step-by-step
illustration

## A step-by-step guide.

Imagine the variable y indicates risk of disease recurrence in a unit of
time. We have a prediction model that quantifies this risk given a
patient’s age, disease severity level, sex, and whether the patient has
comorbidity.

The package comes with two exemplary datasets. dev\_data and val\_data.
We use the dev\_data as the development sample and the val\_data as the
external validation sample.

``` r
data(dev_data)
data(val_data)
```

dev\_data has 500 rows. val\_data has 400 rows.

Here are the first few rows of dev\_data:

|      age | severity | sex | comorbidity |   y |
|---------:|---------:|----:|------------:|----:|
| 55.16048 |        0 |   0 |           1 |   1 |
| 51.57261 |        1 |   0 |           0 |   0 |
| 62.80492 |        0 |   0 |           1 |   0 |
| 61.34763 |        1 |   1 |           1 |   1 |
| 58.35566 |        0 |   1 |           0 |   0 |
| 53.52201 |        1 |   0 |           0 |   1 |
| 44.56900 |        0 |   0 |           0 |   0 |

We use the development data to fit a logistic regression model as our
risk prediction model:

``` r
reg<-glm(y~sex+age+severity+comorbidity,data=dev_data,family=binomial(link="logit"))
summary(reg)
```

    ## 
    ## Call:
    ## glm(formula = y ~ sex + age + severity + comorbidity, family = binomial(link = "logit"), 
    ##     data = dev_data)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.2905  -0.8274  -0.6371   1.1583   2.1063  
    ## 
    ## Coefficients:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -1.732851   0.564476  -3.070  0.00214 ** 
    ## sex          0.557104   0.223634   2.491  0.01273 *  
    ## age          0.005254   0.010643   0.494  0.62152    
    ## severity    -0.557294   0.227586  -2.449  0.01434 *  
    ## comorbidity  1.091916   0.209945   5.201 1.98e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 602.15  on 499  degrees of freedom
    ## Residual deviance: 560.40  on 495  degrees of freedom
    ## AIC: 570.4
    ## 
    ## Number of Fisher Scoring iterations: 4

Given this, our risk prediction model can be written as:

$\\bf{ logit(p)=-1.7329+0.5571\*sex+0.0053\*age-0.5573\*severity+1.0919\*comorbidity}$.

First, let’s compare the ROC and mROC in the development data

``` r
pred<-predict.glm(reg, type='response')

library(pROC)
```

    ## Warning: package 'pROC' was built under R version 4.1.1

    ## Type 'citation("pROC")' for a citation.

    ## 
    ## Attaching package: 'pROC'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     cov, smooth, var

``` r
dev_roc<-roc(response=dev_data[,'y'], predictor=pred)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
plot(dev_roc)
title("ROC in the development dataset")
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

And now the much-awaited mROC using these data. Note that we use the
line function to add the mROC on top

``` r
dev_mroc<-mROC(p=pred)
```

``` r
plot(dev_roc)
lines(dev_mroc, col="red")
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Important note: the statistical inference on comparing mROC and ROC cannot be used for internal validation. Such a test is made for external validation.

Now lets calculate the predicted probabilities for each subject in the
validation dataset given the prediction equation.

``` r
pred<-predict.glm(reg,newdata = val_data, type="response")

summary(pred)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.1007  0.1855  0.2794  0.2892  0.4051  0.5715

Using the package pROC, let’s draw the validation ROC curve

``` r
val_roc<-roc(response=val_data[,'y'], predictor=pred)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
plot(val_roc)
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

And now the much-awaited mROC using these data. Note that we use the
line function to add the mROC on top

``` r
val_mroc<-mROC(p=pred)
```

Notice that the mROC function only requires the vector of predicted
probabilities.

To compare the ROC and mROC plots, we juxtapose them next to each other:

``` r
plot(val_roc)
lines(val_mroc, col="red")
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Here, it is obvious that the mROC and ROC curve are not compatible,
indicating that the model is not calibrated.

``` r
res<-mROC_inference(val_data[,'y'],pred)

res
```

    ## Mean calibration statistic (A):3.17053e-315(Obs>Pred) (p:0.01187)
    ## mROC/ROC equality statsitic (B):0.436176 (p:0.00196)
    ## Unified statistic:4.69855 (df:0.02999348,p:0.0004679777)
