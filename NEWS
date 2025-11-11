 Summary of important user-visible changes for statistics 1.7.0:
-------------------------------------------------------------------

 Important Notice: 1) `mad`, `mean`, `median`, `std`, `var` functions
                      shadow core Octave's respective functions prior
                      to Octave v9.1
                   2) incompatibility with the `nan` package

 New functions:
 ==============

 ** fcnnpredict

 ** fcnntrain

 ** fitcdiscr

 ** fitcgam

 ** fitcnet

 ** fitcsvm

 ** loadmodel

 ** signrank

 New classdefs:
 ==============

 ** ClassificationDiscriminant

 ** ClassificationGAM

 ** ClassificationNeuralNetwork

 ** ClassificationPartitionedModel

 ** ClassificationSVM

 ** CompactClassificationDiscriminant

 ** CompactClassificationGAM

 ** CompactClassificationNeuralNetwork

 ** CompactClassificationSVM

 Improvements:
 =============

 ** BinomialDistribution: fixed computations on truncated distribution
                          methods (github issue #128)

 ** ClassificationKNN: new methods and various bug fixes

 ** finv: fix erratic output for very large DF2 (savannah bug #66034)

 ** fitcknn: support for cross-validation

 ** NegativeBinomialDistribution: fixed computations on truncated distribution
                                  methods (github issue #128)

 ** plot method for probability distribution objects updated to support 'cdf'
    and 'probability' PlotType options (github issue #129)

 ** PoissonDistribution: fixed computations on truncated distribution
                         methods (github issue #128)


 Summary of important user-visible changes for statistics 1.7.1:
-------------------------------------------------------------------

 Important Notice: 1) `mad`, `mean`, `median`, `std`, `var` functions
                      shadow core Octave's respective functions prior
                      to Octave v9.1
                   2) incompatibility with the `nan` package

 Improvements:
 =============

 ** ClassificationPartitionedModel: add input validation for cvpartition
                                    object (github issue #160)

 ** ClassificationNeuralNetwork: use named variable for timing training
                                 (github issue #169)

 ** kstest: add support for probability distribution objects, fix previous
            regression related to failing test (github issues #164, #165, #167),
            add more BISTs and a DEMO

 ** RegressionGAM: add iteration limit to the private fitGAM function to avoid
                   infinite loop with Netlib BLAS/LAPACK (github issue #160)


 Summary of important user-visible changes for statistics 1.7.2:
-------------------------------------------------------------------

 Important Notice: 1) `mad`, `mean`, `median`, `std`, `var` functions
                      shadow core Octave's respective functions prior
                      to Octave v9.1
                   2) incompatibility with the `nan` package

 Improvements:
 =============

 ** fcnnpredict, fcnntrain: various bug fixes, allow MacOS users to enable
                            support for OpenMP (github issues #168, #171, #172)

 ** editDistance: use correct index types to avoid compiler warnings
                  (github issue #171)


 Summary of important user-visible changes for statistics 1.7.3:
-------------------------------------------------------------------

 Important Notice: 1) `mad`, `mean`, `median`, `std`, `var` functions
                      shadow core Octave's respective functions prior
                      to Octave v9.1
                   2) incompatibility with the `nan` package

 New functions:
 ==============

 ** glmval

 Improvements:
 =============

 ** glmfit: fixed output and updated functionality and compatibility


 Summary of important user-visible changes for statistics 1.7.4:
-------------------------------------------------------------------

 Important Notice: 1) `mad`, `mean`, `median`, `std`, `var` functions
                      shadow core Octave's respective functions prior
                      to Octave v9.1
                   2) incompatibility with the `nan` package

 Improvements:
 =============

 ** ClassificationDiscriminant, ClassificationKNN, ClassificationNeuralNetwork:
                    fix type of returning labels and ClassNames to be the same
                    as input Y, add custom display and subsref/subsassgn methods

 ** confusionchart: fix empty figures in demos (github issue #178)

 ** ConfusionMatrixChart: properly update cdata in figure (github issue #178)

 ** fitdist: fix occasionally failing BISTs (github issue #174)

 ** kmeans: fix empty figures in demos (github issue #179)


 Summary of important user-visible changes for statistics 1.7.5:
-------------------------------------------------------------------

 Important Notice: 1) `mad`, `mean`, `median`, `std`, `var` functions
                      shadow core Octave's respective functions prior
                      to Octave v9.1
                   2) incompatibility with the `nan` package

 New functions:
 ==============

 ** createns

 ** multiway

 New classdefs:
 ==============

 ** cvpartition

 ** ExhaustiveSearcher

 ** hnswSearcher

 ** KDTreeSearcher

 Improvements:
 =============

 ** chi2pdf, gampdf: fix INF handling (#203)

 ** crosstab: allow row vectors as input, fix numeric input ordering (#184)

 ** cvpartition: old style class has been replaced by classdef implementation.
                 'K-fold' partitioning is now randomized, and 'repartition'
                 method works for both 'k-fold' and 'holdout' partition types.

 ** fitdist: fix NaN handling as missing values (#203)

 ** fpdf, finv: fix output when DoF tend to infinity (#203)

 ** fullfact: fix MATLAB compatibility (#212)

 ** grpstats: fix MATLAB compatibility (#215)

 ** knnsearch, rangesearch: fix slow kd-tree implementation (#151)

 ** vartest, vartest2: fix computation of right tail p-value (#183)

 ** x2fx: fix variable shadowing in nested loops (#204)

 ** ztest: fix second output for one-tailed tests (#199, #200)


 Summary of important user-visible changes for statistics 1.7.6:
-------------------------------------------------------------------

 Important Notice: 1) `mad`, `mean`, `median`, `std`, `var` functions
                      shadow core Octave's respective functions prior
                      to Octave v9.1
                   2) incompatibility with the `nan` package

 New functions:
 ==============

 ** factoran (#257)

 ** nanmean (#239)

 Improvements:
 =============

 ** Update documentation add demos to all distrbution objects. Several minor
    bug fixes.

 ** ClassificationDiscriminant: improve documentation, allow changing 'Cost'
                                and 'Prior' properties in existing objects

 ** ClassificationKNN: 'predict' can use either exhaustive or kdtree method

 ** crosstab: fix NaN input handling, improve efficiency (#214, #219)

 ** crossval: fix MATLAB compatibility (#217, #223, #229)

 ** cvpartition: fix MATLAB compatibility (#227)

 ** multiway: simplified function signature, improved efficiency (#258, #259)

 ** nansum: fix MATLAB compatibility to handle 'vecdim' and 'all' options

 ** pdist2: fix returning the smallest/largest K values for single output

 ** tcdf: improve speed and fix edge case in processing small DF (#232)

 ** ztest2: improve input validation (#226)


 Summary of important user-visible changes for statistics 1.7.7:
-------------------------------------------------------------------

 Important Notice: 1) `mad`, `mean`, `median`, `std`, `var` functions
                      shadow core Octave's respective functions prior
                      to Octave v9.1
                   2) incompatibility with the `nan` package

 Improvements:
 =============

 ** Update documentation to all classification objects. Minor bug fixes (#263)

 ** Update documentation to ConfusionMatrixChart class (#285)

 ** cvpartition: fix error in stratified input vector (#268)
                 fix handling of missing observations (#269)

 ** pca: fix erroneous Hotelling's T^2 statistics when applying weights (#279)
