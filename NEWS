 Summary of important user-visible changes for statistics 1.8.0:
-------------------------------------------------------------------

 Important Notice: 1) Minimum Octave version required is 9.1.0
                   2) Add dependency to the datatypes package
                   3) Incompatibility with the `nan` package
                   4) Update 'libsvm' library to version 3.36
                   5) Using `tablicious` as a drop in replacement
                      for `datatypes` might cause issues.  Future
                      development of statistical functions to enable
                      support for tables and categorical arrays will
                      be based on the implementations in `datatypes`

 New functions:
 ==============

 ** makima

 Improvements:
 =============

 ** fillmissing: add support for the "makima" interpolation method

 ** svmpredict: add support for one class posterior probablity

 Removed functions: (available in core Octave)
 ==================

 ** mad

 ** mean

 ** median

 ** std

 ** var


 Summary of important user-visible changes for statistics 1.8.1:
-------------------------------------------------------------------

 Important Notice: 1) Soon after Octave 11.1 is officially released the
                      'statistics' package will update its dependency to
                      Octave 11, so that further development can rely on
                      advancements to core statistics and data functions
                      that will be available with Octave 11.
                   2) Incompatibility with the `nan` package
                   3) Using `tablicious` as a drop in replacement
                      for `datatypes` might cause issues.

 New functions:
 ==============

 ** dummyvar

 ** parseWilkinsonFormula

 Improvements:
 =============

 All random generator distribution functions now accept empty size vectors
 as third input argument for MATLAB compatibility.

 ** confusionmat: add support for new data types

 ** crosstab: add support for new data types

 ** cvpartition: implement 'summary' method

 ** datasample: improve input validation

 ** fillmissing: char arrays do not have standard missing values

 ** friedman: return ANOVA results in a table

 ** grp2idx: fix MATLAB compatibility, add support for new data types

 ** grpstats: add support for tables and categorical arrays, add ploting
              functionality, fix remaining MATLAB compatiblity issue

 ** hnswSearcher: improve performance in building index and query operations

 ** ismissing: fix MATLAB compatibility, char arrays do not have standard
               missing values

 ** KDTreeSearcher: improve query speed

 ** knnsearch: improve input validation and avoid infinite recursion

 ** pcacov: improve input validation

 ** pdist, pdist2: improve performance, add fast euclidean algorithms

 ** rangesearch: improve input validation and avoid infinite recursion

 ** rmmissing: fix MATLAB compatibility, char arrays do not have standard
               missing values

 ** standardizeMissing: fix MATLAB compatibility, char arrays do not have
                        standard missing values

 ** stepwisefit: fix MATLAB compatibility, replace legacy Draper-Smith
                 algorithm, results may differ from previous implementation

 ** tabulate: add support for tables and categorical arrays

 ** wishrnd: allow for non-integral DF under certain constrains


 Summary of important user-visible changes for statistics 1.8.2:
-------------------------------------------------------------------

 Important Notice: 1) Update dependency to datatypes 1.2.0
                   2) Incompatibility with the `nan` package
                   3) Using `tablicious` as a drop in replacement
                      for `datatypes` might cause issues.


 Improvements:
 =============

 ** anova1: properly handle NaN values as group labels

 ** boxplot: add support for multi grouping variables

 ** grp2idx: improve integration with datatypes class objects

 ** makima: improve MATLAB compatibility

 ** parseWilkinsonFormula: add 'equations' mode

 ** sampsizepwr: fix computations for single-tailed distributions

 ** signtest: properly handle NaN values


 Summary of important user-visible changes for statistics 1.8.3:
-------------------------------------------------------------------

 Important Notice: 1) Minimum Octave version required is 11.1.0
                   2) Update dependency to datatypes 1.2.3
                   3) Incompatibility with the `nan` package
                   4) Using `tablicious` as a drop in replacement
                      for `datatypes` might cause issues.


 Improvements:
 =============

 ** ClusterCriterion: improve speed

 ** cmdscale: add support for second argument

 ** DaviesBouldinEvaluation: fix index

 ** harmmean, trimmean: improve MATLAB compatibility

 ** linkage: improve speed and fix unused variable

 ** pca, princomp: improve speed

 ** randsample: fix edge cases

 ** stepwisefit: improve speed

 ** squareform, tabulate: improve MATLAB compatibility

 ** RegressionGAM.predict: correct ySD computation using residual variance


 Summary of important user-visible changes for statistics 1.8.4:
-------------------------------------------------------------------

 Important Notice: 1) Update dependency to datatypes 1.2.6
                   2) Incompatibility with the `nan` package
                   3) Using `tablicious` as a drop in replacement
                      for `datatypes` might cause issues.


 New functions:
 ==============

 ** LinearModel: object-oriented linear regression model class,
                 returned by `fitlm`

 ** anova: object-oriented interface for analysis of variance,
           built on `anova1`, `anova2`, and `anovan` with a
           MATLAB-compatible object API

 Improvements:
 =============

 ** anova1: add support for categorical arrays as grouping variables

 ** fitlm: now returns a LinearModel object; the previous
           `[TAB, STATS] = fitlm (...)` output form is no longer supported

 ** glmfit: fix a variable error (GitHub PR #440)

 ** parseWilkinsonFormula: fix power-term resolution from a base table column

