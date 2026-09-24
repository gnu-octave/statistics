 Summary of important user-visible changes for statistics 1.9.0:
-------------------------------------------------------------------

 Important Notice: 1) Update dependency to datatypes 1.3.3.
                   2) The 29 probability distribution classes have moved into
                      the `prob` namespace.  The flat class names are gone;
                      this is a rename, not an alias.
                   3) The generalized additive model classes now fit boosted
                      trees by default, as MATLAB does.  Every fitted value,
                      score and loss moves.  The spline engine remains
                      available as 'FitMethod', 'splines'.
                   4) `fitcnet` and `fitrnet` now default to the LBFGS
                      solver, and `fitcnet` to rectified hidden layers and a
                      softmax output.  Every default fit returns different
                      weights.
                   5) Models of the compact classes saved by an earlier
                      version no longer load and must be re-saved from the
                      model they were compacted from.
                   6) Incompatibility with the `nan` package.
                   7) Using `tablicious` as a drop in replacement for
                      `datatypes` might cause issues.

 Breaking changes:
 =================

 This release changes behaviour that working code may depend on.  Each item
 below says what breaks and what to write instead.

 ** anova2:
    NaN values are refused again, as they always were and as MATLAB does.

 ** anovan:
    the returned table has nine columns where it had seven; 'Singular?' is
    inserted at column 4, moving 'Mean Sq.' to column 5, and 'Eta Sq.' and
    'Part. Eta Sq.' are appended at 8 and 9.  A model with any random factor
    keeps the interactions that factor appears in, so `p` and the table are
    longer and the table gains eight further columns; p-values change for
    every random effects model.  An interaction is named 'X1:X2', and the
    random factor marker is no longer carried in the names returned in
    `atab`, `stats.varnames` or the expected mean square expressions.
    Reading the table by name from `atab(1,:)` is unaffected.

 ** boxplot:
    the quartiles are taken from `quantile` at its own default method, the
    one `prctile` and MATLAB use, rather than from core's `statistics`.
    Boxes, whisker fences, the returned statistics matrix and what counts as
    an outlier all move.

 ** ClassificationGAM, RegressionGAM:
    the classes fit boosted trees by default; pass 'FitMethod', 'splines' for
    the previous engine.  Eleven spline options now raise unless that engine
    is selected, and a second output from `RegressionGAM.predict` raises
    under the default engine.  The boosted classifier fits by local scoring,
    so every multi-predictor classification GAM fits to different numbers.  A
    tree leaf now holds at least five observations.  `Interactions` reports
    the K-by-2 matrix of fitted predictor pairs, not the request that was
    passed; use `IntMatrix` to recover the argument.  `ScoreTransform`
    defaults to 'logit' and the raw score is the log-odds pair MATLAB
    reports.  `RegressionGAM.fitBoosted` is no longer a public method.

 ** ClassificationKNN:
    `DistanceWeight` holds the name of the weight rather than a function
    handle, and now applies: nothing read it before, so 'inverse' and
    'squaredinverse' changed no prediction.  The documented spelling
    'squaredinverse' is accepted.

 ** ClassificationSVM:
    `loss` returns different numbers, having scaled every loss by the class
    labels so that only the classification error was right.  `margin` takes Y
    in the form the class labels are in, as documented, where it required the
    +1/-1 coding.  `Alpha` is populated for every kernel and holds unsigned
    magnitudes, `Beta` is the primal coefficient vector for a linear kernel
    only, `SupportVectorLabels` is s-by-1 of +/-1, and `IsSupportVector` is
    logical.  No classification result changes.

 ** cluster, clusterdata:
    'MaxClust' yields exactly N clusters however the merge heights tie, and a
    'Cutoff' on the inconsistency coefficient requires every node below to be
    under the cutoff.  Both returned too few clusters.

 ** ClusterCriterion, CalinskiHarabaszEvaluation, DaviesBouldinEvaluation,
    GapEvaluation, SilhouetteEvaluation:
    the cluster evaluation classes are value classes, not handles, so `e2 =
    e1` copies; assign results back, as `e = addK (e, 4)`.
    `ClusterCriterion` is abstract and can no longer be constructed directly.

 ** CompactClassificationSVM:
    the class drops `IsSupportVector`, which was sized to the training set
    and which MATLAB's compact class does not carry.

 ** cvpartition:
    `test` and `training` return one column per requested set; given a vector
    of set indices they returned one column per set per pass.

 ** evalclusters:
    a matrix of clustering solutions with no `KList` is numbered by the
    clusters each column holds, not by column position, so criterion values
    and `OptimalK` change.  `SilhouetteEvaluation.ClusterSilhouettes` holds
    the mean silhouette of each cluster rather than every observation's
    value.

 ** factoran:
    the function is rewritten and its output signature is now MATLAB's,
    `[lambda, psi, T, stats, F]`.  The third output was the factor scores and
    is now the rotation matrix; scores moved to the fifth, so code taking the
    third output as scores silently receives an m-by-m matrix.  The fit is
    maximum likelihood by default where it was principal axis factoring,
    which remains available as 'Extraction', 'paf'.

 ** fcnntrain, fcnnpredict:
    both take activation names in two arguments, `Activations` and
    `OutputLayerActivation`, rather than numeric codes in a single vector.
    `fcnnpredict` takes the layer weights and biases, not the model structure
    `fcnntrain` returns.  The `Alpha` argument and field are removed.  The
    'elu' activation saturates at -1, as it is defined to, where its negative
    arm carried the leaky rectifier's 0.01 scale.

 ** fitcnet:
    'Solver' defaults to 'lbfgs', MATLAB's solver, where it defaulted to
    'sgd', so every default fit returns different weights.  'LearningRate'
    belongs to the epoch loop and now raises unless a solver is named, while
    'GradientTolerance', 'LossTolerance' and 'StepTolerance' are now accepted
    by default; `TrainingHistory` and `ConvergenceInfo` carry the fields of
    the solver that ran.  `fitrnet` is new in this release.  `fitcnet`
    additionally defaults 'Activations' to 'relu', 'OutputLayerActivation' to
    'softmax' and 'LearningRate' to 0.003.  Initial weights are drawn from
    Octave's generator, so repeated fits no longer return an identical model
    and `rand ('seed', s)` governs a fit.  `ConvergenceInfo` reports the
    value the fit ended at as a scalar and carries the series beside it as
    `History`.  `ModelParameters` reports the fit rather than the network;
    read `Mdl.LayerWeights` and `Mdl.LayerBiases`.

 ** fitdist:
    the 'theta' option defaults to 0 for the generalized Pareto, as MATLAB
    assumes it, where it defaulted to 1; `mle` changes with it.

 ** fpdf:
    the density at x = 0 is 1 at df1 = 2 and Inf below it, where it was zero
    whatever the numerator degrees of freedom.

 ** friedman:
    the second output is a cell array with a header row, as MATLAB's is and
    as `anova1` and `anova2` already were; read `tbl{2,2}` rather than
    `tbl.SS(1)`.  The table is displayed by default.

 ** gpfit:
    the second argument is `alpha`, the confidence level, not `theta`, a
    known location; `paramhat` is 1-by-2 and `paramci` 2-by-2.  Fit a known
    location by shifting the data.  Note `gpfit (x, 0.05)` was a fit at
    location 0.05 and is now a 95% confidence level.

 ** gplike:
    `params` is the two element vector [k, sigma] and `acov` is 2-by-2, where
    it was 3-by-3 with a zero row and column.

 ** hmmviterbi, hmmestimate:
    `hmmviterbi` returns the true maximum probability path, having scored a
    spurious transition out of the last state, and `hmmestimate` no longer
    counts a phantom transition out of an assumed initial state.  Both return
    different values.

 ** LinearModel:
    a weighted fit's `LogLikelihood` carries the 0.5 * sum (log (w)) term and
    counts only nonzero weight observations, following R rather than MATLAB;
    `AIC`, `AICc`, `BIC` and `CAIC` move with it.  Unweighted and robust fits
    are unchanged.

 ** LinearModel:
    `Formula` is a `LinearFormula` object rather than a structure or a
    character vector.  Term order follows the variable order of the data
    rather than the alphabetical order of their names, and no longer depends
    on the order the formula names them, so coefficient names move with it.
    The omitted reference level of a character or string grouping column is
    the level the data presents first, not the alphabetically first one, so
    coefficients change name and meaning; a categorical column is unaffected.
    A model with no intercept and a categorical predictor gives the first
    categorical an indicator for every level by both routes, so its
    coefficients are the group means; this is a deliberate deviation from
    MATLAB, which drops the reference level.

 ** mle:
    the 'normal' and 'lognormal' distributions return the maximum likelihood
    estimate, `std (x, 1)`, where they returned the unbiased `std (x, 0)`.

 ** nbinfit:
    the confidence bounds are no longer clamped to the parameter space, a
    bound outside it being what tells the caller the normal approximation has
    broken down.

 ** negloglik:
    the distribution object method returns the negative log likelihood, as
    documented and as MATLAB returns it; it returned the log likelihood, the
    wrong sign.

 ** pca:
    `tsquared` is 0 where the fit has one degree of freedom or fewer, and NaN
    for dropped rows, returned at the full observation count.

 ** prob distribution objects:
    the 29 classes are now `prob.NormalDistribution` and so on; the flat
    names are gone, so `class (pd)` and `isa (pd, ...)` take the qualified
    name.  `makedist`, `fitdist`, `cdf`, `pdf`, `icdf`, `random` and `mle`
    take the distribution's name and are unaffected.  `DistributionName`
    holds the distribution's name ('Normal'), not the class name.
    `prob.WeibullDistribution` renames `lambda` and `k` to `A` and `B` and
    `prob.RayleighDistribution` renames `sigma` to `B`.  `ParameterValues` is
    a cell for the piecewise linear and multinomial classes, `x` and `Fx` of
    the piecewise linear class are rows, and the 29 classes now share an
    abstract base, `prob.ProbabilityDistribution`, as MATLAB's do.

 ** proflik:
    the distribution object method computes a true profile likelihood, re-
    maximizing over the remaining parameters at each point, where it held
    them fixed at their fitted values and so gave wrong results for every
    multi-parameter distribution.  The default grid is 101 values over the
    98% confidence interval when the selected parameter is the only one
    estimated, and 21 otherwise.

 ** regress_gp:
    the prediction intervals were wrong in both branches and are now correct,
    and both branches return the lower bound in column 1; code reading
    Yint(:,1) as an upper bound must be changed.

 ** signrank:
    two differences agreeing to within eps (x) + eps (y) rank as tied and one
    smaller than that is dropped, as MATLAB does, so `signedrank`, `zval` and
    the p-value move on such data.  `stats` carries `signedrank` where it
    carried `sign`, which belongs to `signtest`, and `zval` is empty for the
    exact test.

 ** signtest:
    the test decision `h` is returned as a logical, as documented and as
    MATLAB returns it, where it was a double.

 ** tabulate:
    the percentage column is count / total * 100 taken literally, so data
    with no observations gives NaN where it gave zero.

 ** The compact classes:
    `CompactClassificationSVM` drops `ModelParameters` and `IsSupportVector`,
    `CompactClassificationNeuralNetwork` drops six properties describing the
    fit, and `CompactClassificationGAM` drops seven, MATLAB's own compact
    classes carrying none of them.  Read any of them off the model that was
    compacted from.  A compact SVM, neural network or GAM model saved by an
    earlier version therefore no longer loads, and `loadmodel` raises
    'invalid model'.  The compact regression and Gaussian process classes are
    new in this release and never carried them.

 ** The distribution functions:
    the 130 *cdf, *pdf and *inv functions validate the class of their numeric
    arguments.  Continuous distributions and every *inv accept double and
    single only; the discrete *cdf and *pdf additionally accept the integer
    types and promote to double, so the result is always a probability.
    Logical and character arrays are rejected.  This fixes silently wrong
    answers: `geocdf (int32 (0:4), 0.3)` returned [0 1 1 1 1], `nbincdf`
    returned zeros and `binopdf` returned Inf.

 ** The random generators:
    the 39 generators and the distribution objects' `random` method treat a
    negative dimension as zero and return an empty array, as core Octave and
    MATLAB do, instead of raising.

 ** The supervised learner classes:
    `ScoreTransform` and `ResponseTransform` hold the name of the transform
    rather than a function handle, so `strcmp (Mdl.ScoreTransform, 'none')`
    works and `predict` applies the transform for you.  `Prior` is a 1-by-K
    row in `ClassNames` order.  `W` weighs each class by its prior rather
    than being uniform.  An observation is dropped only for a missing
    response, so `size (X, 1)`, `numel (W)` and `NumObservations` now agree
    and `RowsUsed` is empty when nothing was dropped.  The `Standardize`
    property is gone from eleven classes; test `! isempty (Mdl.Mu)` instead.
    `Cost` refuses anything that is not a square matrix.  `BinEdges` is a
    cell on the cross-validated classifier, and a column cell on the GAM
    classes.  `RowsUsed` is a logical mask on `ClassificationGAM`,
    `ClassificationSVM` and `RegressionGAM`, where it was a double.  `edge`
    normalizes its 'Weights' within each class to that class's prior.  The
    cross-validated models are read only apart from `Cost`, `Prior` and the
    transforms, either assignment now reaching every fold, and their
    transforms are applied where they used to be inert.
    `CrossValidatedModel` holds the learner's short name, not the class name.
    An observation that no fold held out gets a missing label with NaN
    scores, a deliberate deviation from MATLAB, which names the first class.

 ** tiedrank:
    an array is ranked along its first dimension, as MATLAB does, instead of
    being refused; `tieadj` carries one entry per column.

 ** trimmean:
    a non-scalar percentage is refused, as in MATLAB; it was accepted and
    returned an untrimmed mean behind a warning.

 ** ttest, ttest2:
    `Alpha` is validated; a negative, zero, unit or vector value was
    accepted, a negative one returning a NaN confidence interval together
    with 'do not reject'.

 ** unifit:
    the function returns MATLAB's four outputs [ahat, bhat, aci, bci], where
    it packed both estimates into one vector and both intervals into one
    matrix.  X may now be a matrix, negative data is accepted, and the upper
    endpoint's interval was inverted.

 New functions:
 ==============

 ** andrewsplot:
    Andrews plot of multivariate data, with grouping, four standardizations
    and a quantile mode.

 ** ansaribradley:
    Ansari-Bradley test for equal dispersions, with an exact permutation
    p-value or a tie corrected normal approximation.

 ** biplot:
    biplot of principal component coefficients or factor loadings.

 ** ClassificationDiscriminant.cvshrink:
    cross-validate the discriminant over a grid of `Gamma` and optionally
    `Delta`.

 ** ClassificationDiscriminant.mahal, ClassificationDiscriminant.logp,
    ClassificationDiscriminant.nLinearCoeffs:
    squared Mahalanobis distance to each class mean, the log of the density
    summed over the classes, and the number of predictors kept at a
    regularization threshold; also on the compact class.

 ** ClassificationGAM.addInteractions, RegressionGAM.addInteractions:
    fit interaction terms onto a model that already carries its univariate
    ones.

 ** ClassificationGAM.resume, RegressionGAM.resume:
    add trees to a fitted model, continuing in the phase that ran last.

 ** ClassificationPartitionedLinear, ClassificationPartitionedKernel,
    RegressionPartitionedLinear, RegressionPartitionedKernel:
    the cross-validated counterparts of the linear and kernel learners.

 ** ClassificationPartitionedModel.kfoldLoss,
    ClassificationPartitionedModel.kfoldMargin,
    ClassificationPartitionedModel.kfoldEdge:
    the class had none of the three; `kfoldLoss` takes 'classiferror',
    'classifcost', 'mincost' or a function handle, the margin based losses
    not being offered.

 ** CompactLinearModel:
    a fitted linear model without its training data, returned by
    `LinearModel`'s new `compact` method.

 ** CompactRegressionGAM:
    the compact counterpart of `RegressionGAM`.

 ** copulastat, copulaparam, copulafit:
    complete the copula family, which already shipped `copulacdf`,
    `copulapdf` and `copularnd`.

 ** coxphfit:
    Cox proportional hazards regression, opening survival analysis, with
    stratified models, the counting process form and seven residual types.

 ** dbscan:
    density-based spatial clustering of applications with noise.

 ** discardSupportVectors:
    empty the support vectors of a linear kernel model and leave `Beta` and
    `Bias` to predict, on the four SVM classes; the engine's own copy is
    collapsed, so the memory is freed.

 ** dwtest:
    Durbin-Watson test for autocorrelation in regression residuals, with an
    exact p-value by Imhof's method.

 ** ecdfhist:
    histogram built from the output of `ecdf`.

 ** fitckernel, ClassificationKernel, fitrkernel, RegressionKernel:
    Gaussian kernel models for large data, fitted linearly in a randomized
    feature space, with `resume`.

 ** fitclinear, ClassificationLinear, fitrlinear, RegressionLinear:
    linear models for data with many predictors, with ridge or lasso
    penalties and six solvers; a vector 'Lambda' fits one model per strength
    inside a single object, narrowed afterwards by `selectModels`.

 ** fitcnb, ClassificationNaiveBayes, CompactClassificationNaiveBayes:
    naive Bayes classification with the 'normal', 'kernel', 'mvmn' and 'mn'
    predictor distributions, mixable one per predictor.

 ** fitcox, CoxModel:
    the Cox proportional hazards model as an object, with `survival`,
    `hazardratio`, `coefci`, `linhyptest`, `plotSurvival` and
    `discardResiduals`.

 ** fitglm, GeneralizedLinearModel:
    generalized linear regression as a fitted model object over five
    families, the GLM counterpart of `fitlm` and `LinearModel`.  A binomial
    response given with 'BinomialSize' is the number of successes, as MATLAB
    reads it, and `predict` returns the probability of success.

 ** fitglme, GeneralizedLinearMixedModel:
    generalized linear mixed effects models by penalized quasi-likelihood,
    with four fit methods.

 ** fitlme, fitlmematrix, LinearMixedModel:
    linear mixed effects models from a Wilkinson formula extended with random
    effects terms or from design matrices, by maximum likelihood or REML.

 ** fitnlm, NonLinearModel:
    nonlinear regression by iterative least squares as a fitted model object,
    with an error model and robust fitting.

 ** fitrgp, RegressionGP, CompactRegressionGP:
    Gaussian process regression with ten covariance functions, four explicit
    bases, and `postFitStatistics`.

 ** fitrnet, RegressionNeuralNetwork, CompactRegressionNeuralNetwork:
    neural network regression against the mean squared error with an identity
    output layer.

 ** fitrsvm, RegressionSVM, CompactRegressionSVM:
    support vector machine regression by epsilon-insensitive loss, with
    'SVMtype', 'nu_svr' as an Octave extension.

 ** gamtrain, gampredict, gamboosttrain, gamboostpredict, gamboostpairs,
    gamboostinter, __lbfgs__, __bhtsne__, __knnselect__, __knnbrute__:
    new compiled functions: the spline and boosted-tree GAM engines, the
    limited-memory BFGS solver, the Barnes-Hut t-SNE summation and the
    nearest-neighbour kernels.  They are not meant to be called directly.

 ** glyphplot:
    star glyph plot of multivariate data; Chernoff faces are not supported.

 ** gplotmatrix:
    matrix of grouped scatter plots, with grouped histograms on the diagonal.

 ** hmmdecode:
    posterior state probabilities of a hidden Markov model by a scaled
    forward-backward recursion.

 ** hmmtrain:
    maximum likelihood estimation of hidden Markov model parameters by Baum-
    Welch or Viterbi training.  Baum-Welch is verified bit-for-bit against
    MATLAB; Viterbi training carries two documented deviations.

 ** iforest, IsolationForest:
    isolation forest anomaly detection, with an `isanomaly` method.

 ** invpred:
    inverse prediction from a simple linear regression, with Fieller
    confidence bounds that need not be finite.

 ** jbtest:
    Jarque-Bera hypothesis test of composite normality.

 ** johnsrnd:
    random arrays from the Johnson system of distributions.

 ** KernelDistribution:
    nonparametric kernel smoothing distribution object, created by `fitdist
    (x, 'Kernel')`.

 ** kfoldfun:
    apply a function across the folds of a cross-validated model, on
    `ClassificationPartitionedModel` and `RegressionPartitionedModel`.

 ** kmedoids:
    k-medoids clustering, with the 'pam' and 'small' algorithms.

 ** ksdensity:
    kernel smoothing density estimate, with four kernels, five output
    functions and bounded support; 'Censoring' is not yet implemented.

 ** lasso, lassoglm:
    lasso and elastic net regularized regression by cyclic coordinate
    descent, for least squares and for the five GLM families, with cross-
    validation and the 1-SE rule.

 ** lillietest:
    Lilliefors composite goodness of fit test for the normal, exponential and
    extreme value families.

 ** LinearFormula:
    model formula object, the class of the `Formula` property of a
    `LinearModel` and of a `GeneralizedLinearModel`.

 ** LinearModel.plot, LinearModel.plotInteraction, LinearModel.anova,
    LinearModel.step, LinearModel.compact:
    five new methods.

 ** lof, LocalOutlierFactor:
    local outlier factor anomaly detection, with an `isanomaly` method.

 ** mdscale:
    nonclassical metric and nonmetric multidimensional scaling over five
    stress criteria.

 ** mlecov:
    asymptotic covariance matrix of maximum likelihood estimators, from the
    observed Fisher information.

 ** mnrval:
    predict category probabilities and confidence bounds from a multinomial
    logistic regression model, the companion to `mnrfit`.

 ** mvksdensity:
    multivariate product-kernel smoothing density estimate.

 ** mvregress, mvregresslike:
    multivariate linear regression by maximum likelihood, with three
    algorithms for missing responses.

 ** nanvar, nanstd, nanmedian, nancov:
    complete the nan* family alongside the existing `nanmean`, `nansum`,
    `nanmax` and `nanmin`.

 ** nlinfit, nlparci, nlpredci:
    Levenberg-Marquardt nonlinear least squares with a numeric Jacobian, and
    confidence and prediction intervals for its coefficients and predictions.

 ** nnmf:
    nonnegative matrix factorization, by alternating least squares or
    multiplicative updates.

 ** ocsvm, OneClassSVM:
    one-class support vector machine anomaly detection, with an `isanomaly`
    method.

 ** parallelcoords:
    parallel coordinates plot of multivariate data.

 ** paretotails:
    piecewise distribution object with generalized Pareto tails and an
    empirical middle; the 'kernel' middle segment is not supported.

 ** partialcorr, partialcorri:
    linear and rank partial correlation coefficients, with p-values.

 ** pearsrnd:
    random arrays from the Pearson system of distributions.

 ** perfcurve:
    ROC and other classifier performance curves, returning the curve,
    thresholds, AUC and optimal operating point, with bootstrap bounds.

 ** ppca:
    probabilistic principal component analysis, handling missing data by
    expectation-maximization.

 ** probplot:
    probability plot of a sample against eight reference distributions, with
    censoring and frequency weights.

 ** RegressionPartitionedModel:
    the regression counterpart of `ClassificationPartitionedModel`, returned
    by `crossval` on a `RegressionGAM`, `RegressionGP`,
    `RegressionNeuralNetwork` or `RegressionSVM`.

 ** rica, ReconstructionICA:
    reconstruction independent component analysis for feature extraction,
    with a `transform` method.  `Mu` and `Sigma` are columns and `FitInfo`
    carries the whole minimisation history, as in MATLAB.

 ** robustcov:
    robust multivariate covariance and mean estimate, with the 'fmcd' and
    'ogk' methods; 'olivehawkins' is not implemented.

 ** robustfit:
    robust M-estimator linear regression by iteratively reweighted least
    squares, with nine weight functions.

 ** rocmetrics:
    receiver operating characteristic metrics object, with per-class AUC and
    the `addMetrics`, `average` and `plot` methods.

 ** rotatefactors:
    rotate a factor loading matrix by the orthomax family, promax, or
    procrustes rotation, the last defaulting to 'Type', 'oblique' as in
    MATLAB.

 ** scatterhist:
    scatter plot with marginal histograms or kernel densities.  As in MATLAB,
    the marginal bars point toward the scatter by default and 'Location'
    names the corner the histograms occupy, not the scatter.

 ** sparsefilt, SparseFiltering:
    sparse filtering for feature extraction, with a `transform` method.

 ** spectralcluster:
    spectral clustering using the eigenvectors of a similarity graph
    Laplacian.

 ** StableDistribution, stblpdf, stblcdf, stblinv, stblrnd, stblfit,
    stbllike:
    the stable distribution in the Nolan S0 parameterization, by numerical
    inversion of the characteristic function; `makedist ('Stable')` was
    previously an unsupported stub.

 ** statset, statget:
    the options structure pair used across the package's iterative
    algorithms; seven functions documented an 'Options' argument and nothing
    could build one.

 ** stepwiseglm:
    fit a `GeneralizedLinearModel` by stepwise term selection, with the
    selection trace in its `Steps` property.

 ** stepwiselm:
    stepwise linear regression over five criteria, with term hierarchy
    enforced on both addition and removal.

 ** The classification margin and edge surface:
    fifteen new methods complete it: every full classifier now has `margin`,
    `edge`, `resubPredict`, `resubMargin`, `resubEdge` and `resubLoss`, and
    every compact classifier has `margin` and `edge`.  A compact classifier
    gets no resubstitution method, as in MATLAB.

 ** tsne:
    t-distributed stochastic neighbor embedding, with the 'exact' and
    'barneshut' algorithms.

 Improvements:
 =============

 ** anova:
    the `anova1`, `anova2` and `anovan` backends all run silently, so fitting
    never prints a table or opens a figure.  `varianceComponent` no longer
    returns an interval for a negative variance estimate, and a saturated
    model reports its error term as zero on zero degrees of freedom rather
    than producing no F-ratios at all.

 ** anova1:
    `p` and `F` are NaN for a design with no variation at all and `p` is 0
    where the fit is exact, both having been the wrong way round.  A
    categorical level that no observation uses leaves every field of `stats`,
    so they can be indexed together; this is a deliberate deviation from
    MATLAB, whose own fields disagree in length there.  The `vartype`
    argument is documented as an Octave extension.

 ** anovan:
    the table reports a 'Singular?' flag and the degrees of freedom that can
    actually be estimated for each term.  A design that is not of full column
    rank now takes the minimum-norm least squares solution, where the
    coefficients came from a triangular solve of a singular system and the
    sums of squares were nonsense.  A continuous factor raised to a power is
    squared from the raw predictor rather than the centred one, which could
    put Type II and Type III out by orders of magnitude.  `stats.termcols`
    counts one column per level, and 'sstype' accepts 'H' beside 'h'.

 ** boxplot:
    a single observation, an entirely missing variable and an empty input are
    plotted rather than failing inside the function, integer observations are
    accepted, and a grouped plot with a one-observation group no longer
    fails.

 ** canoncorr:
    rank deficient input was solved against a singular triangular factor, so
    the coefficients ran away and the canonical correlations came back wrong
    with nothing said.  The factorisation now pivots, deficient columns are
    dropped, and a warning is raised as MATLAB raises one.

 ** ClassificationDiscriminant, CompactClassificationDiscriminant:
    all six discriminant types are available where only 'linear' was.
    `DiscrimType`, `Gamma` and `Delta` may be assigned after fitting, re-
    deriving `Sigma`, `LogDetSigma` and `Coeffs` without refitting, and
    `Delta` now eliminates predictors, the new `DeltaPredictor` reporting the
    value at which each drops out.  `predict` applies the `ScoreTransform`
    the model holds, where it ignored it entirely.  Assigning `Prior` or
    `Cost` no longer destroys the model, `LogDetSigma` is accurate on nearly
    collinear predictors, and the class can be loaded at all, `load_model`
    having raised 'invalid model' whatever the model held.

 ** ClassificationGAM, CompactClassificationGAM:
    'Interactions' given as 'all' or as a count built no interaction term at
    all.  Assigning `Cost` was refused on every model.  `ClassNames` keeps
    the type of Y, where it was a cell of character vectors whatever the
    response was.  A numeric response coded 1 and 2 gave a NaN intercept and
    every score with it.  `predict` applies the `ScoreTransform`, which was
    documented, settable and read by nothing.

 ** ClassificationKNN:
    `NumNeighbors`, `Distance`, `DistanceWeight`, `BreakTies`, `IncludeTies`
    and `DistParameter` may be assigned after fitting and are live.
    `predict` applies the score transform once rather than once per query
    row, and its third output honours the cost matrix where it returned 1 -
    posterior.  A model fitted with any metric a kd-tree cannot search could
    not be reloaded; eight of the twelve metrics were affected.  'CacheSize'
    is carried for compatibility.

 ** ClassificationNeuralNetwork:
    the trained parameters are reachable through `LayerWeights` and
    `LayerBiases`; before this a network could be trained with no way to see
    what it had learned.  `loss` accepts every loss MATLAB documents for the
    class and defaults to 'mincost'.  A saved model came back carrying
    another model's numbers, `savemodel` having written 21 of 29 properties
    into a placeholder object.

 ** ClassificationPartitionedModel:
    an explicit `Prior` reached no fold on the discriminant, SVM and neural
    network backings, and a non-default `Cost` reached no fold on the SVM
    backing, which enters it into the box constraint.  A standardized
    `ClassificationKNN` refitted every fold on the raw scale.  `kfoldPredict`
    discarded the predictions of a 'Holdout' partition, raised for a GAM
    backing and refused a cost for a network backing.  The class gains `W`
    and a `disp` of its own, and takes a `Prior` on a naive Bayes backing.

 ** ClassificationSVM, CompactClassificationSVM:
    `predict` and `resubPredict` return the expected cost as a third output,
    and `kfoldPredict` does for an SVM-backed cross-validated model, where it
    raised.  MATLAB returns the cost matrix read transposed; we return the
    row, which is documented.  `resubLoss` returned the wrong number or none,
    and `resubPredict` selected the training rows with a bare mask, so a
    three-predictor model was asked about a single column.  The classes
    report `Prior`, `Cost`, `W`, `CategoricalPredictors` and
    `ExpandedPredictorNames`.

 ** crossval:
    the partition covers only the observations actually trained on; with any
    missing value present each fold trained on rows its own fit then
    discarded.  A function handle over a predictor and a response is no
    longer rejected, the last data variable having been dropped, which is why
    `plsregress`'s 'CV' never worked.

 ** cvpartition:
    'Stratify' and 'GroupingVariables' can be combined, and
    'GroupingVariables' applies to 'Holdout' and 'Leaveout' as well as
    'KFold', where all three accepted it and silently ignored it, giving the
    caller the leakage the option exists to prevent.  A scalar SVAL to
    `repartition` never worked, and seeding is now confined to the call.

 ** evalclusters:
    the returned object has MATLAB's shape and reports MATLAB's names.  The
    silhouette criterion honours the 'Distance' it is given, where it always
    used 'sqeuclidean' and every metric returned the same numbers.

 ** ExhaustiveSearcher, KDTreeSearcher, hnswSearcher, knnsearch, rangesearch:
    the exhaustive search no longer forms the whole distance matrix for the
    four Minkowski-family metrics, and the searches take their K smallest by
    partial selection: a 4000-point self-search at K = 5 falls from 4.5 s to
    1.9 s and searches that could not allocate now run.  The kd-tree search
    is about five times faster, resolving the metric once per query.
    `knnsearch` and `rangesearch` accept a per-call metric.  Asking for more
    neighbours than the data holds is answered with all of them.  Single data
    stays single and every other class converts to double, which fixes
    integer data silently rounding each coordinate difference.  Raising
    `SearchSetSize` on `hnswSearcher` returned the same neighbour for every
    query past a threshold.

 ** fcnntrain, fcnnpredict:
    every network the package trained descended the wrong gradient: the mean
    squared error layer differentiated sum (y - t)^2 as 2y - t rather than 2
    (y - t), so training settled at 2y = t and the network converged on half
    its target.  The backward pass of a hidden activation layer recovered the
    incoming gradient by dividing by its own input, read past the end of the
    array for any widening network, and applied the wrong local derivative
    for six of the eight activations.  The `Accuracy` and `Loss` vectors
    described no network that ever existed and were then stored behind as
    many zeros as the epoch count.  Training visited the samples in a fixed
    order, so with sorted labels the weights swung between the classes
    instead of settling.  Weights are drawn on a range set by the fan-in
    rather than always on [-1, 1].

 ** fitcgam, fitrgam:
    nine new boosted-tree options matching MATLAB in name and default, and an
    interaction asked for by count or 'all' is chosen by testing each
    candidate pair against the predictor phase's residuals and ranking by
    p-value, with 'MaxPValue' setting the cut.  The GAM classes report
    `FitMethod`, `TreeModel`, `ReasonForTermination`, `BinEdges` and
    `PairDetectionBinEdges`, and `ModelParameters` reports MATLAB's thirteen
    fields under the boosted-tree engine.

 ** fitcnet:
    'Prior' and 'Cost' are accepted, having been rejected as unknown
    parameters, so there had been no way to weigh the classes at all.
    'Activations' accepts 'none', which the documentation listed and the code
    rejected.  A model fitted with 'Standardize' set to true was trained on
    the raw predictors while `predict` standardized its own, so
    resubstitution accuracy on a badly scaled problem was 0.275 and is now
    0.925.

 ** fitcsvm:
    a model fitted with 'Standardize' set to true was trained on the raw
    predictors while `predict` standardized its own, so resubstitution
    accuracy on a badly scaled problem was 0.500, chance, and is now 0.912.

 ** glmfit:
    the returned coefficient vector was clobbered by the binomial Anscombe
    residual computation, returning garbage coefficients whenever the `stats`
    output was requested and taking `glmval`'s confidence intervals with it.

 ** glmval:
    accept a numeric 0 or 1 for the 'simultaneous' option and 'BinomialSize'
    as an alias for 'size'.

 ** gpfit:
    the boundary warnings say which boundary was reached and that the
    confidence intervals are withheld.  NaN and Inf observations propagate
    into the estimates rather than being dropped with a warning.

 ** hygecdf, poisstat, wblpdf, ClassificationNeuralNetwork:
    error messages and help text corrected: `hygecdf` named its third
    argument `k` where `hygepdf` names it `M`, `poisstat` named SIGMA,
    `wblpdf`'s three signature lines named `wblinv`, and the `Solver`
    property read 'Gradient Descend'.

 ** LinearModel, fitlm:
    a robust fit reported the wrong `LogLikelihood` and information criteria,
    and a robust fit given weights scored the wrong residuals, so a precisely
    measured observation was less likely to be downweighted.  A variable the
    model does not use is no longer counted as one of its predictors.
    `Steps` is populated by `stepwiselm` and `step`, where it was always
    empty.  `removeTerms` could not name a term on a table model whose
    formula uses only some of the table's variables, `step` refused any model
    holding a power term, and `stepwiselm` silently dropped one from a
    starting model or bound.  Thresholds letting a term enter and leave at
    once are refused, where the search never returned.  `VariableNames` and
    `PredictorNames` are column cell arrays.

 ** mle:
    the custom distribution interface is implemented, having been a stub that
    printed a message and silently fitted a normal.  The 'bernoulli' and
    'unif' families could not be fitted at all, and an empty 'Frequency'
    raised instead of meaning no frequencies.

 ** mnrfit:
    rewritten for MATLAB compatibility: nominal models with any number of
    categories, having errored for more than two, plus ordinal and
    hierarchical models, four links, the full `stats` structure and the
    `EstDisp` option.

 ** ModelParameters:
    the property reaches its final shape on every learner that has one.
    `ClassificationKNN` and `ClassificationDiscriminant` gain it, having been
    the only learners without one; the neural network pair reports the fit as
    it was asked for rather than the network that came out of it;
    `RegressionGP` reports its starting values; and the six partitioned
    classes report the learner's parameters under their own tags.

 ** optimalleaforder:
    a custom 'Transformation' is reachable, having been documented and
    implemented but refused before the code that uses it ran, and a two-leaf
    tree no longer crashes.

 ** paramci, proflik:
    `paramci` reports one column per parameter on `BinomialDistribution`,
    holding the fixed N at its own value, which also unblocks `proflik`
    there, and a fixed parameter is addressable by name again on every
    distribution.  `proflik` returns a third output, defaults `pnum` to the
    first free parameter, keeps the profiled-out parameters inside their own
    range, and refuses a parameter with no confidence interval.  On
    `BurrDistribution` it profiled a censored likelihood, the frequencies
    having landed in the censoring slot.

 ** pca:
    three defects in the handling of missing data and weights, contributed by
    Swayam Shah (GitHub PR #451).  With 'Rows' set to 'all' the guard that
    refuses NaN input fired only when every column held a NaN.

 ** pdist, pdist2:
    no longer return a negative distance: 'cosine', 'correlation' and
    'spearman' put a row at -2.2e-16 from itself.

 ** prob.PiecewiseLinearDistribution:
    `mean`, `var` and `std` of a truncated distribution are computed exactly
    segment by segment, where quadrature integrated across the density's
    jumps and lost four to five digits.

 ** randsample:
    added to the package INDEX, so the function is listed in `pkg describe`
    and the online documentation; it shipped but was undiscoverable there.

 ** regress:
    the coefficient standard errors and confidence intervals are computed
    from the QR factor directly rather than by inverting R' * R, which
    squared the condition number of the design.  On the Longley benchmark
    they now match the NIST certified values to 6e-13 relative, previously
    2e-8, with no spurious singularity warning.

 ** regression_ftest, regression_ttest:
    `regression_ttest` had no working call form, `cov` of two vectors
    returning a matrix, and its fitted values were evaluated at mean (x)
    instead of at x, so the residual sum of squares was really the total.
    `regression_ftest` refused the commonest test of all, dropping a single
    predictor.  Positive tests have been added to both.

 ** RegressionGAM:
    `predict` computes ySD and yInt from the model's own terms rather than
    the stored predictors, so a model built with 'Interactions' or a
    'Formula' no longer reports a standard deviation several times too large,
    or raises.  A single Inf passed validation, and a scalar 'Knots', 'Order'
    or 'DoF' was never expanded.  `savemodel` wrote a text file where every
    other class writes a binary one.

 ** RegressionGAM, ClassificationGAM:
    the spline fitting is compiled and each predictor's design is factorised
    once for the whole fit: fitting a classifier is four times faster and
    cross-validating one between five and nine times.  A rank deficient
    design no longer divides by a singular value that is rounding noise.

 ** sampsizepwr:
    the 'r' test type for the correlation coefficient was unreachable by any
    route, and the sample size came from Fisher's transformation of the
    difference between the two correlations, which is right only when the
    null correlation is zero, the single value the function refused.

 ** slicesample:
    a log density supplied through 'logpdf' was logged a second time, so the
    chain sampled an entirely different distribution without complaint;
    drawing from a standard normal returned values with a mean around -8000.
    Three missing input checks are added.

 ** svmpredict:
    scoring a single observation corrupted the heap and aborted the
    interpreter; two or more rows were fine.  The result pointers referred to
    a temporary that was freed before the prediction loop wrote to it.

 ** The copula family:
    `copulapdf` gains the Gaussian and Student's t families and `copularnd`
    the Frank and Gumbel-Hougaard families, so every family MATLAB supports
    now works in every copula function; only Clayton previously worked in all
    five.  The Ali-Mikhail-Haq and Farlie-Gumbel-Morgenstern Octave
    extensions work throughout as well.

 ** The distribution name wrappers:
    `cdf`, `pdf`, `icdf`, `random`, `makedist`, `fitdist` and `mle` accept
    the same set of distribution names, matched ignoring case, spaces and
    hyphens; they previously accepted disjoint spellings, so `makedist
    ('Extreme Value')` and `cdf ('ExtremeValue', ...)` were both errors while
    their opposites worked.

 ** The distribution objects:
    `prob.KernelDistribution` hides `NumParameters`, `ParameterNames` and
    `ParameterDescription`, MATLAB carrying none of the three, and
    `prob.StableDistribution` no longer lists `ParameterCI`; all remain
    readable by name.  The classes list `Truncation` and `IsTruncated` where
    MATLAB lists them.  `MultinomialDistribution.random` and
    `LoguniformDistribution.random` were broken in every call form, and
    `KernelDistribution`'s `pdf`, `cdf` and `icdf` errored on an empty query.

 ** The learner classes:
    each reaches the property and method surface MATLAB gives it.
    `Intercept`, `CategoricalPredictors`, `ExpandedPredictorNames`, `W`,
    `BetweenSigma`, `KernelParameters`, `BoxConstraints`, `OutlierFraction`,
    `Nu` and `BinEdges` are added where each was missing, and the ten classes
    with a MATLAB counterpart declare `HyperparameterOptimizationResults`,
    read-only and always empty.  Every classification learner accepts a
    response given as a character matrix, which the documentation always
    listed and which had never worked.  Every classifier with a settable
    `Cost` validates the matrix as MATLAB does and accepts a struct of
    `ClassNames` and `ClassificationCosts`; `Prior` accepts the matching
    struct form.  The seventeen supervised classes enforce access with
    property attributes rather than hand-written `subsref` and `subsasgn`
    whitelists.  Setting `ScoreTransform` to 'none' or 'identity' made
    `predict` raise on four classes, and assigning it raised on five others,
    so the property could not be set at all.  'ismax' and 'symmetricismax'
    mark the largest score of each observation, where the maximum was taken
    down the column.

 ** The paired argument parser:
    `crossval`, `grpstats`, `rmmissing`, `stepwisefit` and `scatterhist`
    parse Name-Value arguments through `parsePairedArguments` from the
    `datatypes` package, and the statistics package's own private
    `pairedArgs` is gone.

 ** tiedrank:
    accepts MATLAB's fourth argument, a tolerance, which defaults to exact
    comparison so every existing call is unaffected.

 ** violin:
    a row vector is read as one variable, a cell of row vectors works,
    integer and logical observations are accepted, and a single observation,
    an empty input and a call with no arguments report what is wrong instead
    of failing inside the function.

 ** wblplot:
    a non-finite observation is rejected; `wblplot ([1, Inf, 2, 3])` never
    returned, creating nine graphics objects per iteration for as long as it
    was left alone.

 Removed functions:
 ==================

 ** mvtcdfqmc:
    the quasi-Monte-Carlo integrator behind `mvtcdf` and `mvncdf` is now the
    private helper `__mvtcdfqmc__` and can no longer be called directly.
    MATLAB ships no such function, so nothing portable depended on it.  Call
    `mvtcdf` or `mvncdf` instead.

 Summary of important user-visible changes for statistics 1.9.1:
-------------------------------------------------------------------

 Improvements:
 =============

 ** ClassificationKNN:
    a custom distance function handle works again.  The constructor probed
    the handle against the response instead of the predictors and sized the
    result against the number of response rows, so a correctly written
    `distfun` was rejected and no custom metric could be used at all.  A
    handle assigned to `Distance` after fitting was accepted without being
    probed and failed later inside `predict`, and `disp` could not print
    one.  The `D2 = distfun (ZI, ZJ)` signature is now documented in the
    property's help text.

 ** documentation:
    `doc/statistics.qch` is rebuilt with pkg-octave-doc 0.7.7.  Double-
    clicking the package in the GUI Documentation tab opens an overview
    listing every documented name under its category, rather than whichever
    category came first, and the one formula that reached the page as raw
    TeX now renders as the plain text form of its docstring.  The online
    pages for `NonLinearModel` and `gmdistribution` were published as plain
    function pages, listing none of their properties or methods, because a
    member documented outside texinfo made the class renderer fail and fall
    back silently; both are rendered as classes again, restoring 105 and 72
    documented members.

 Summary of important user-visible changes for statistics 1.9.2:
-------------------------------------------------------------------

 New functions and methods:
 ==========================

 ** rica, ReconstructionICA:
    a `'NonGaussianityIndicator'` name-value pair, and the read-only property
    of the same name, giving one sign per learned feature: `+1` seeks a
    super-Gaussian feature and `-1` a sub-Gaussian one.  It sets the sign that
    feature's contrast term carries in the objective, so it moves the fit.  The
    default is all `+1`, which is the fit as it stood.

 Improvements:
 =============

 ** binoinv:
    an exactly attained probability is reached rather than stepped past, so
    the median of a symmetric binomial is right: `binoinv (0.5, 5, 0.5)`
    returned 3 and `binoinv (0.5, 2001, 0.5)` returned 1001.  A `p` of 1
    returns `n` rather than wherever the CDF first saturates, which was 206
    for 500 trials at 0.25.  An answer above 500 no longer errors on a row
    vector, on a vector `n`, or with a NaN among the probabilities.

 ** prob.BinomialDistribution:
    `median` and `icdf` at 0.5 follow the `binoinv` fix, so
    `median (prob.BinomialDistribution (5, 0.5))` is 2 rather than 3 on
    platforms whose `betainc` is less accurate.

 ** nbinpdf:
    a density is returned wherever one exists.  The three factors of the
    formula left the range of a double before their product did, so
    `nbinpdf (1000, 1001, 0.5)` was `Inf` where the density peaks at 0.0089,
    and `NaN` past the mode.

 ** nbincdf:
    follows the `nbinpdf` fix.  The lower tail sums the density, so it read 1
    from 308 failures on for 1001 successes, where the value is 4.2e-39, `NaN`
    past the mode, and was not monotone.

 ** nbininv:
    an exactly attained probability is reached rather than stepped past.  For
    `r` successes at a probability of 0.5 the median is exactly `r-1`, and
    `nbininv (0.5, r, 0.5)` returned `r` for 188 of the 399 values of `r` from
    2 to 400, `nbininv (0.5, 101, 0.5)` among them.

 ** prob.NegativeBinomialDistribution:
    `median` and `icdf` follow the `nbininv` fix, so
    `median (prob.NegativeBinomialDistribution (101, 0.5))` is 100 rather
    than 101.

 ** documentation:
    `disp` and `display` are hidden on every class that defines one, which was
    already true of 141 of the 147 and is now true of all of them.  They no
    longer appear in `methods (CLASS)` or in the online pages for `anova`,
    `paretotails`, `NonLinearModel`, `ConfusionMatrixChart`, `gmdistribution`
    and `rocmetrics`.  Calling them is unaffected, and `NonLinearModel`'s
    `subsref` is hidden with them.

    `gmdistribution`'s methods are documented: `cdf`, `cluster`, `mahal`,
    `pdf`, `posterior`, `random` and the static `fit` carried no help text at
    all and now answer `help gmdistribution.METHOD`.

    `help` headings name the package rather than a category that does not
    exist: `copularnd`, `createns` and `normplot` were headed `Function File`,
    and `loadmodel` was headed `ClassificationSVM`, a class it is not part of.

    `NonLinearModel` and `LinearMixedModel` document their properties: all 23
    and all 17 carried no help text at all and now answer
    `help CLASS.PROPERTY`.  `NonLinearModel`'s `ErrorModelInfo` is hidden,
    MATLAB not carrying it in its public interface, so
    `properties ('NonLinearModel')` lists 23 names rather than 24; reading
    `MDL.ErrorModelInfo` is unaffected.

    `CoxModel`, `GeneralizedLinearMixedModel` and `gmdistribution` document
    their properties: all 16, all 14 and all 14 carried no help text and now
    answer `help CLASS.PROPERTY`, and `gmdistribution`'s constructor is
    documented.  `gmdistribution`'s `NlogL` is hidden, MATLAB hiding the same
    alias, so `properties ('gmdistribution')` lists 14 names rather than 15;
    reading `OBJ.NlogL` is unaffected and `NegativeLogLikelihood` carries the
    same value.

    `ConfusionMatrixChart`'s fifteen property set methods move into a hidden
    block, so `set.XLabel` and its fourteen siblings no longer appear as
    methods in the online pages.  Assigning the properties is unaffected.

    `ReconstructionICA` and `SparseFiltering` document their properties: the
    same eight on each carried no help text and now answer
    `help CLASS.PROPERTY`.  `FitInfo` states that the fitting trajectory is
    this implementation's own, so its length and its iteration counts differ
    from MATLAB's under either solver.

    Six more property set methods move into hidden blocks, in
    `ClassificationKernel`, `ClassificationLinear`, `RegressionKernel`,
    `RegressionLinear`, `RegressionGP` and `CompactRegressionGP`, so
    `set.ScoreTransform` and `set.ResponseTransform` no longer appear as
    methods in the online pages.  Assigning the properties is unaffected.

    `RegressionGP` and `RegressionGAM` document their constructors, and every
    compact class answers `help CLASS.CLASS` with real text where eight of the
    ten had returned a one-word source comment.  All ten hide the constructor,
    as MATLAB does, so `methods` no longer lists it; `compact` is named in
    every one of their class descriptions as the way to create one.

    Every demo that draws random numbers seeds each generator stream it uses,
    so `demo` gives the same result twice and a rebuild of the online pages
    diffs only where something really changed.  Octave seeds `rand`/`randn`,
    `randg`, `rande` and `randp` separately, and the legacy `('seed', N)` form
    selects a different generator, so many demos that looked seeded were not.

    The 29 `+prob` distribution class constructors answer `help`: each carried
    no help text at all and now gives its call signature, names its parameters
    and states their defaults.  They stay hidden, `makedist` being the usual
    way to create a distribution object.

 ** johnsrnd:
    `quantiles` may be given as a 2-by-4 matrix, the first row four strictly
    increasing, evenly spaced standard normal quantiles and the second row the
    data quantiles at them; a four-element vector still means the fixed
    `[-1.5, -0.5, 0.5, 1.5]` points.  The `SL` family returned a curve of the
    wrong sign whenever `delta` was negative, `lambda` having been fixed at 1
    rather than taken as `sign (delta)`.  A `NaN` among the quantiles is
    refused instead of being carried into the fit.

 ** fitclinear, fitrlinear:
    a fit that stops because the line search cannot improve the objective says
    so, reporting `TerminationCode` -11 and `'Unable to find a step decreasing
    the objective.'` as MATLAB does.  Every exit but three was reported as
    `'Iteration limit exceeded.'`, so a fit that gave up after 59 iterations of
    a limit of 1000 claimed to have run out of iterations.

 ** partialcorr:
    a coefficient is NaN where the controlling variables explain either of the
    two variables completely, as MATLAB returns.  The residual variance was
    tested against exact zero, which floating-point arithmetic does not
    produce, so the coefficient was computed from rounding error and was not
    reproducible between builds: `partialcorr ([1 2; 3 4; 5 5], [2; 4; 6])`
    gave an off-diagonal of -0.13 here and -0.089 elsewhere, where MATLAB gives
    NaN.

 ** lasso:
    `'Standardize'` is forced to false, with a warning, when `'Intercept'` is
    false.  The two were accepted together and the predictors were standardised
    anyway.

 ** mixed effects models:
    `fitlme`, `fitlmematrix`, `fitglme`, `LinearMixedModel` and
    `GeneralizedLinearMixedModel` fit markedly more accurately.  Both profiled
    deviances are now optimised with their closed form gradients rather than by
    finite differences, every matrix is symmetrised before it is factored, and
    the fitting path forms no explicit inverse.  On a balanced one-way fixture,
    where REML must equal the ANOVA estimate exactly, the variance component
    was 6.9e-07 out and is now 1.4e-09.

 Summary of important user-visible changes for statistics 1.9.3:
-------------------------------------------------------------------

 Important Notice: Update dependency to datatypes 1.4.2.

 New functions and methods:
 ==========================

 ** fitctree, ClassificationTree, CompactClassificationTree:
    a CART binary decision tree for classification, over the new compiled
    `treetrain` and `treepredict` engine.  It grows on a numeric predictor
    matrix with a numeric, logical, character or cellstring response, honours
    `Prior`, `Cost` and `Weights`, merges leaves and estimates the cost
    complexity pruning sequence, and reports the property surface MATLAB
    reports.  Its methods are `predict`, `prune`, `predictorImportance`,
    `nodeVariableRange`, `view`, `margin`, `edge`, `loss`, the four `resub`
    forms, `compact`, `crossval`, `cvloss` and `savemodel`.  `compact` returns
    a `CompactClassificationTree`, the tree without the observations it was
    fitted on, carrying the 33 properties MATLAB reports.  `crossval` returns
    a `ClassificationPartitionedModel` over compact folds, taking `KFold`,
    `Holdout`, `Leaveout` or `CVPartition`, each fold grown with the parent's
    class names, prior, cost and observation weights.  `cvloss` returns the
    cross-validated loss of the tree and of its subtrees, with the standard
    error, the leaf count of each subtree and the pruning level chosen by
    `'TreeSize'`; its standard error is the standard error over the folds,
    which is not the quantity MATLAB reports under that name.

    Categorical predictors are taken with `'CategoricalPredictors'`,
    `'MaxNumCategories'` and `'AlgorithmForCategorical'`.  Such a predictor is
    split into two sets of levels, reported in `CutType`, `CutCategories` and
    `CategoricalSplit` and shown by `view` and `nodeVariableRange`; a level a
    node did not see stops an observation there.  Two classes order the
    levels, more classes search every partition up to `MaxNumCategories`
    levels and otherwise take the best of the PullLeft, PCA and OVAbyClass
    heuristics.  Splits match MATLAB R2024a but for ties and, on nodes with
    more than `MaxNumCategories` levels, the heuristic splits and their sides.

    Not implemented, and refused rather than ignored: surrogate splits,
    predictor subsampling, the `'twoing'` split criterion, the `'impurity'`
    pruning criterion, and binning.

 ** fitrtree, RegressionTree, CompactRegressionTree:
    a CART binary decision tree for regression, over the same compiled engine.
    It grows on a numeric predictor matrix and a numeric response, honours
    `Weights` and `QuadraticErrorTolerance`, merges leaves and estimates the
    cost complexity pruning sequence, and reports the 36 properties MATLAB
    reports.  Its methods are `predict`, `prune`, `predictorImportance`,
    `nodeVariableRange`, `view`, `loss`, `resubPredict`, `resubLoss`,
    `compact`, `crossval`, `cvloss` and `savemodel`.  `compact` returns a
    `CompactRegressionTree`, carrying the 28 properties MATLAB reports;
    `crossval` returns a `RegressionPartitionedModel` over compact folds, each
    grown with the parent's growth parameters and a slice of its observation
    weights; `cvloss` returns the cross-validated mean squared error of the
    tree and of its subtrees, with the same standard error caveat as the
    classification tree.

    Categorical predictors are taken with `'CategoricalPredictors'` and
    `'MaxNumCategories'`, the levels being ordered as for two classes.  Splits
    match MATLAB R2024a but for ties and the sides a split's levels take.

    Not implemented, and refused rather than ignored: surrogate splits,
    predictor subsampling and binning.

 ** fitcensemble, ClassificationEnsemble, ClassificationBaggedEnsemble,
    CompactClassificationEnsemble:
    ensembles of decision trees for classification, boosted with AdaBoostM1,
    AdaBoostM2, GentleBoost, LogitBoost, LPBoost, TotalBoost or RUSBoost, or
    bagged.  `predict`, `loss`, `edge`, `margin` with `Learners`, `Mode` and
    `UseObsForLearner`, the `resub*` methods, `compact`, `resume`,
    `removeLearners` and `predictorImportance`; `oobPredict`, `oobLoss`,
    `oobEdge`, `oobMargin` and `oobPermutedPredictorImportance` on the bagged
    ensemble.  The `Subspace` method fits nearest neighbour or discriminant
    learners on random predictor subsets, or on every combination with
    `AllPredictorCombinations`, and averages their scores; weights that are
    not uniform are refused, the learners taking none.  `RUSBoost` grows each
    tree on a weighted, class-balanced sample sized by `RatioToSmallest`.
    A boosting method resamples its observations when `'Resample'` is `'on'`
    or `'FResample'` or `'Replace'` is given, returning the bagged ensemble,
    which records the rows each learner drew and estimates out-of-bag error.
    `'CategoricalPredictors'` is passed on to every tree, which splits those
    predictors into sets of levels, and is carried to compact models and
    cross-validation folds; the `Subspace` method refuses them, its learners
    taking none here, where MATLAB passes them to nearest neighbour learners
    and drops them from discriminant ones.

    `LPBoost` grows each tree on the dual weights of a linear program over the
    trees, whose solution gives the learner weights, and stops when the
    smallest edge is within `'MarginPrecision'` of the program's value.
    `TotalBoost` takes one quadratic step towards the least relative entropy
    under edge constraints, stops when no tree can improve the margin, and
    maximises the smallest margin with GLPK.  Both follow MATLAB R2024a; where
    a program has several solutions the trees grown next, or the weights
    chosen, can differ.  MATLAB R2024a stops as if `'MarginPrecision'` were at
    least 0.01, where a smaller value is taken as given here.  A perfect weak
    learner is kept and ends the fit, where MATLAB drops it; GentleBoost
    applies `LearnRate` to the fit; `resume` continues AdaBoostM2 exactly.

 ** fitrensemble, RegressionEnsemble, RegressionBaggedEnsemble,
    CompactRegressionEnsemble:
    ensembles of regression trees, boosted with LSBoost or bagged.  `predict`
    and `loss` with `Learners`, `Mode` and `UseObsForLearner`, `resubPredict`,
    `resubLoss`, `compact`, `resume`, `removeLearners` and
    `predictorImportance`; `oobPredict`, `oobLoss` and
    `oobPermutedPredictorImportance` on the bagged ensemble, which a boosting
    method returns when `'Resample'` is `'on'` or `'FResample'` or `'Replace'`
    is given.  `'CategoricalPredictors'` is passed on to every tree and
    carried to compact models and folds.

    `regularize` finds non-negative lasso weights for the trees over a path of
    `Lambda` values, `shrink` returns a compact ensemble of the trees a column
    of those weights keeps, and `cvshrink` cross-validates that choice; all
    three are on the bagged ensemble too.  `regularize` returns the minimum of
    the lasso objective, where MATLAB's solver can stop short of it at small
    `Lambda`.

 ** ClassificationPartitionedEnsemble, RegressionPartitionedEnsemble:
    cross-validated ensembles, from `fitcensemble` or `fitrensemble` given
    `CrossVal`, `KFold`, `Holdout`, `Leaveout` or `CVPartition`, or from the
    `crossval` method of an ensemble.  `kfoldPredict`, `kfoldLoss` with
    `Folds` and the modes `average`, `individual` and `cumulative`, `kfoldfun`
    and `resume`, and for classification `kfoldEdge` and `kfoldMargin`.  A
    regression `ResponseTransform` is applied once, where MATLAB applies it
    twice.

 ** fitcecoc, ClassificationECOC, CompactClassificationECOC:
    multiclass classification from binary learners.  A coding matrix says
    which classes each binary learner calls +1, which it calls -1 and which
    sit it out; a new observation goes to every learner and takes the class
    whose row of that matrix its scores match most closely.  Seven learners
    are accepted, `svm`, `tree`, `knn`, `naivebayes`, `discriminant`, `linear`
    and `kernel`, either by name or as a template carrying the options to fit
    them with.  Seven coding designs, seven binary losses and both decoding
    schemes.  A linear or kernel learner keeps no training data, so the fit
    returns the compact model for those two, as MATLAB does.
    `discardSupportVectors` empties the support vectors of every binary
    learner that is a support vector machine on a linear kernel, whose linear
    model stands in for them exactly, and warns when none of them is;
    `selectModels` narrows every binary learner to a subset of its fitted
    regularization strengths.  `'CategoricalPredictors'` is passed as given to
    every binary learner, which codes the predictors its own way; the model
    reports the indices.

 ** ClassificationPartitionedECOC:
    cross-validated ECOC models, from `crossval` or from `fitcecoc` given
    `KFold`, `Holdout`, `Leaveout` or `CVPartition`.  It carries the coding
    matrix, the binary loss and what each observation was to each learner,
    which the general cross-validated class has nowhere to put and without
    which a fold's scores cannot be decoded.  `kfoldPredict`, `kfoldMargin`,
    `kfoldEdge`, `kfoldLoss` and `kfoldfun`.

 ** TreeBagger, CompactTreeBagger:
    random forests of bagged decision trees, for classification and
    regression.  Each tree is grown on a bootstrap sample and chooses every
    split from a random subset of the predictors.  `predict`, `error`,
    `margin` and `meanMargin`, their out-of-bag counterparts `oobPredict`,
    `oobError`, `oobMargin` and `oobMeanMargin`, and `compact`, `growTrees`,
    `append`, `combine` and `setDefaultYfit`.  Out-of-bag predictor importance
    by permutation, the split criterion contributions and split counts of each
    predictor, quantile regression with `quantilePredict`, `quantileError`,
    `oobQuantilePredict` and `oobQuantileError`, and the proximity matrix with
    `fillprox`, `mdsprox`, `proximity` and `outlierMeasure`.
    `'CategoricalPredictors'`, `'MaxNumCategories'` and
    `'AlgorithmForCategorical'` are passed to its trees.

    Labels come back in the type of the response, where MATLAB returns a cell
    array of character vectors.  `oobQuantileError` in `individual` mode
    judges each tree on its own out-of-bag observations, where MATLAB R2024a
    fails, and `growTrees` empties a filled proximity matrix, where MATLAB
    keeps the one computed for the earlier trees.

 ** templateSVM, templateTree, templateKNN, templateNaiveBayes,
    templateDiscriminant, templateLinear, templateKernel, templateEnsemble:
    learner templates, which name a learner and the options it is to be fitted
    with without fitting anything.  `fitcecoc` takes one as its `Learners`;
    `'Learners', 'ensemble'` asks for a LogitBoost ensemble of 100 trees.  The
    default binary loss follows the method, as in MATLAB.  A RUSBoost
    learner's scores, which MATLAB can read with no binary loss and so cannot
    predict with, are divided by the total weight of its trees, as R's ebmc
    and adabag and scikit-learn scale a boosted ensemble's votes, and read
    with `quadratic`.

 ** designecoc:
    the coding design matrix on its own, over seven designs.

 ** treetrain, treepredict:
    the compiled engine for decision trees, shared by the classification and
    the regression learner.  `treetrain` grows a binary tree by recursive
    partitioning, merges the leaf pairs that lower no misclassification and
    computes the cost complexity pruning sequence; `treepredict` evaluates
    one, and both take the level sets of a categorical predictor.  A row
    missing the split predictor descends to neither child and is stopped at
    that node when predicting.  Neither is meant to be called directly.

 ** swtest:
    the Shapiro-Wilk test of normality, by Royston's algorithm AS R94, for
    3 to 5000 values; `'Method', 'shapiro-francia'` gives the Shapiro-Francia
    test instead, for 5 to 5000 values.  MATLAB has no such test.

 ** corr:
    correlation coefficients with their p-values, as MATLAB computes them:
    `'Type'` Pearson, Kendall or Spearman, `'Rows'`, `'Tail'` and
    `'Weights'`.  It shadows core Octave's `corr`, which offers the Pearson
    coefficient alone and no p-value.  Two deviations: a character array is
    refused, where MATLAB correlates the character codes, and an array of
    more than two dimensions is refused, where MATLAB and core Octave answer
    `NaN` for some of them.

 ** fitcknn, ClassificationKNN:
    `'CategoricalPredictors'`, empty or `'all'` as in MATLAB: with every
    predictor categorical the default distance is `'hamming'`.

 ** fitcsvm, fitrsvm, fitcnet, fitrnet (and their classes):
    `'CategoricalPredictors'`, as indices, a logical vector or `'all'`.  Each
    categorical predictor is dummy coded in place, one column per level seen
    in training, named as in `'x1 == 2'` in `ExpandedPredictorNames`, and left
    unstandardized; `X` keeps the predictors as given, and compact models,
    saved models and cross-validation folds carry the coding.  A level not
    seen in training gives no score or prediction.

 Bug fixes:
 ==========

 ** kmeans:
    an empty `'Replicates'` or `'MaxIter'` raised "invalid number of
    replicates" and "invalid MaxIter", where MATLAB takes an empty value as a
    request for the default.  Both now use the default, 1 replicate and 100
    iterations, and an empty `'Replicates'` beside a three-dimensional
    `'Start'` takes the replicate count from its third dimension.

 ** ClassificationGAM, RegressionGAM, fitcgam, fitrgam:
    the bin grids of a boosted-tree fit, `BinEdges` for a predictor with more
    distinct values than cuts and `PairDetectionBinEdges`, placed some cuts
    between different observations from MATLAB's whenever the bin position
    was fractional, and handled repeated values differently.  Each cut is now
    the midpoint of the two observations MATLAB uses, moved up to the next
    distinct value where they are equal, so such fits now agree with MATLAB.

 ** ClassificationGAM, RegressionGAM, fitcgam, fitrgam:
    in an interaction term a missing value added nothing, and a training row
    missing one of the pair's predictors took no part in its trees.  As in
    MATLAB, such a row now stops at the node splitting on its missing
    predictor and takes that node's value, in training and in prediction.
    One case is known to differ: a regression fitted with missing values can
    start its second interaction tree with a different split from MATLAB's.

 ** ClassificationGAM, RegressionGAM, fitcgam, fitrgam:
    interaction terms were fitted on the coarse grid used to detect them, so
    their surfaces differed from MATLAB's.  Pair trees are now fitted to the
    rows, as in MATLAB: a node is cut between two values it holds, a leaf
    keeps at least five rows, trees grow a layer at a time within
    `'MaxNumSplitsPerInteraction'`, and each pair's curvature-weighted mean
    moves into the intercept.  `PairDetectionBinEdges` reports only the
    detection grid.  With `'MaxNumSplitsPerInteraction'` 1 no interaction is
    selected, and a warning says so.

 ** ClassificationGAM, RegressionGAM, fitcgam, fitrgam:
    a training row missing a predictor was left out of the boosted-tree fit.
    It is now fitted, as in MATLAB: it takes no part in that predictor's
    splits and is fitted with the value of the tree's root, and a predictor
    with missing values is held to zero at a missing value, so the intercept
    is the one MATLAB reports.  The spline fit still uses complete rows.

 ** ClassificationGAM, RegressionGAM, fitcgam, fitrgam:
    `'Prior'` did not reach the boosted-tree fit, so it changed only the
    reported prior, and neither learner took a fit-time `'Weights'`.  The fit
    now weighs each class by its prior, spread over its observations in
    proportion to `'Weights'`, as MATLAB does; an empirical prior is each
    class's share of the weight, and `W` reports the weights the fit used.
    Cross-validation folds keep the weights of the rows they hold.
    `'Weights'` is refused with `'FitMethod'` `'splines'`.

 ** stepwisefit:
    data with no complete observation, empty or all missing, failed with
    "operator +: nonconformant arguments"; it now returns the result MATLAB
    does, with `NaN` coefficients, `stats.df0` -1, `stats.dfe` 0 and an empty
    `history`.  X and Y with different numbers of rows now get their own
    message (#484).

 ** ClassificationGAM, RegressionGAM and their compact classes:
    `predict` dropped every row holding a missing value, so it returned fewer
    rows than it was given, and `loss`, `margin` and `edge` failed on such
    data.  Every row is now kept.  Under boosted trees a missing value adds
    nothing from its term; under splines the row has no score (`NaN`), and a
    classifier gives it the class of largest `Prior`.

 ** fitcecoc, ClassificationECOC:
    a model whose binary learners take no observation weights (SVM, KNN, naive
    Bayes, discriminant) failed with "invalid parameter name" on classes of
    unequal size, and its cross-validation failed whenever a fold's classes
    were unequal: the class weights, equal or not, were passed to learners
    that cannot take them.  Such a learner is now given the prior its two
    sides hold, as MATLAB's learners report it; `'Weights'` that vary within a
    class are refused for them with a message that says why.

 ** ClassificationSVM, CompactClassificationSVM:
    `margin`, `loss`, `edge` and `resubLoss` scored the data without the
    standardization the model was fitted with, so on a model fitted with
    `'Standardize'` they gave wrong margins and losses.

 ** ClassificationNeuralNetwork, RegressionNeuralNetwork and their compact
    classes: a row missing a predictor was given a score or a prediction with
    the default ReLU activation, the missing value being rectified to zero; its
    scores are now `NaN`.

 ** ClassificationLinear, ClassificationKernel, ClassificationSVM,
    ClassificationNeuralNetwork, ClassificationKNN, ClassificationECOC,
    their compact classes, and the matching regression models:
    a row missing a predictor was labelled with the first or second class, or
    given finite KNN scores and zero ECOC losses, or a `NaN` prediction.  As in
    MATLAB, a classifier now gives it `NaN` scores (and `NaN` cost) and the
    class of largest `Prior`, and a regression model predicts the weighted lower
    median of the training response.  `fitrgp` keeps such rows in `X` and
    `NumObservations`, as MATLAB does, but leaves them out of the fit.

 ** fitclinear, fitckernel, fitrlinear, fitrkernel, fitrgp:
    `'CategoricalPredictors'` was stored but ignored, so a categorical
    predictor was fitted as a number.  Each categorical predictor is now
    dummy coded in place, one column per level seen in training, named as in
    `'x1 == 2'` in `ExpandedPredictorNames`, and left unstandardized; new data
    is coded the same way, a level not seen in training giving `NaN`.  The
    option also takes `'all'`, and indices beyond the predictors are refused.

 ** ClassificationECOC, ClassificationEnsemble, ClassificationKernel,
    ClassificationLinear, ClassificationSVM, and the fitc functions that
    build them:
    classes given as text were laid out in the order they first appear in
    `Y` instead of sorted, so `ClassNames`, the score columns and per-class
    options such as `Prior`, `Cost` and `RatioToSmallest` disagreed with
    MATLAB and with the other classifiers.  They are now sorted.

 ** ClassificationDiscriminant, ClassificationECOC, ClassificationEnsemble,
    ClassificationGAM, ClassificationKernel, ClassificationKNN,
    ClassificationLinear, ClassificationNeuralNetwork, ClassificationSVM,
    ClassificationTree, and the fitc functions that build them:
    `'ClassNames'` chose which classes to fit but not their order; the model
    now keeps the order given, as MATLAB and ClassificationNaiveBayes do.

 ** ClassificationDiscriminant, ClassificationECOC, ClassificationEnsemble,
    ClassificationGAM, ClassificationKernel, ClassificationKNN,
    ClassificationLinear, ClassificationNaiveBayes,
    ClassificationNeuralNetwork, ClassificationSVM, ClassificationTree,
    TreeBagger, and the fitc functions that build them:
    `'ClassNames'` whose type differed from the labels' failed with an
    internal error or was matched wrongly.  Numeric names for text labels
    now raise "not all 'ClassNames' are present in Y." in every classifier;
    text names match numeric labels by their text, in naive Bayes and
    TreeBagger too, and logical labels as `'true'` and `'false'`; the
    ensembles, ECOC, kernel and linear classifiers now refuse a name absent
    from the labels like the others.

 ** ClassificationTree, RegressionTree, fitctree, fitrtree:
    when `MaxNumSplits` ran out partway through a layer of the tree, the nodes
    left unsplit were the last ones of that layer, whatever their splits were
    worth.  They are now the ones whose splits gain the least, as MATLAB
    documents, so such a tree grows the splits MATLAB grows:
    `fitrtree (meas(:,2:4), meas(:,1), 'MaxNumSplits', 10, 'MinLeafSize', 5)`
    has a resubstitution MSE of 0.0791 where it had 0.0837.

 ** partialcorr:
    a controlling variable shifted by a large constant or scaled by a small
    one, or two variables of very different scales, gave a wrong coefficient
    or `NaN`: `partialcorr (X, [1e8 + t, t.^2])` differed from
    `partialcorr (X, [t, t.^2])` by half.  The data are now centered and
    rescaled before the regression, so the coefficient no longer changes.
    MATLAB's does.  Patch by Dmitri Sergatskov (#483).

 ** ClassificationSVM, fitcsvm:
    a two-class problem whose labels are the numbers -1 and +1 was fitted with
    every class the other way round, so `predict` returned the wrong label for
    every observation and `resubLoss` was 1 where it should have been 0.  Any
    other labels were unaffected, numeric or textual, which is why nothing
    caught it.  The solver is now given the class indices rather than the
    labels themselves.

 ** anova (LinearMixedModel), DFMethod satterthwaite:
    a model whose random-effects terms are confounded, such as
    `(1|g) + (1|g:g2)` where `g2` is a function of `g`, reported DF2 of 0 and a
    meaningless p-value.  It now reports the degrees of freedom of the
    equivalent reduced model.

 ** nanmean, nansum, nanvar, nanstd, nanmedian:
    an empty input that is not 0 by 0 keeps its shape instead of collapsing to a
    scalar.  `nansum (zeros (0,3))` gives `[0 0 0]`, `nanmean ([], 2)` a 0 by 1
    empty and `nanvar (ones (2,0,3,2), 0, 2)` a 2 by 1 by 3 by 2, as MATLAB
    returns; each gave a scalar.  A 0 by 0 input with no dimension given still
    gives the scalar.

 ** nanmax, nanmin:
    an empty input gives an empty result.  The all-NaN mask was assigned into an
    empty array, which grew it, so `nanmax ([])` answered `NaN` and
    `nanmax (zeros (0,3))` a 1 by 3 of `NaN`, values that were never in the
    data.  A dimension vector applied to an empty input no longer errors.

 ** nancov:
    an N-D input is refused with `Inputs must be 2-D.` instead of failing on an
    out-of-bound index.  An explicitly empty `Y`, as in `nancov (x, [])`, is
    taken as given and its length checked, where it was discarded and the
    variance of `X` returned.  `nancov (zeros (0,3), zeros (0,3))` gives a 2 by
    2 of `NaN`, the two arguments being paired into two columns, not a 3 by 3.
    `nancov ([])` gives `NaN`, not an empty.

 ** ranksum:
    `stats` always carries a `zval` field, empty for the exact method, so
    reading it no longer errors on the method that computes no z-statistic.
    Its fields are ordered `ranksum` then `zval`, as MATLAB orders them; `zval`
    came first on every method.

 ** lillietest, jbtest:
    the Monte-Carlo p-value requested through `MCTol` is no longer smoothed by
    an added one, so it reaches exactly 0 and exactly 1 as MATLAB's does.  The
    smallest value either could return was `1 / (reps + 1)`.

 ** jbtest:
    a sample with zero variance no longer aborts with an internal indexing
    error.  Skewness and kurtosis are undefined there, so `jbstat` is NaN, and
    `p` is 0 and `h` is 1.

 ** datasample:
    `k` may be 0, returning an empty of the sampled orientation, as MATLAB
    returns; it was refused as not positive.

 ** makima:
    sample points where `x` or `y` is NaN are dropped in pairs before the fit,
    with a warning, where a single NaN made every returned value NaN.  A NaN in
    one series of a matrix or N-D `y` drops that sample point from all of them,
    as MATLAB does.

 ** ztest:
    a matrix or N-D array is tested along the operating dimension instead of
    being collapsed to a single test, and a scalar is tested rather than
    answering an empty.  `ztest ([1 2; 3 4; 5 6], 3, 1)` gives `[0 0]` and
    `ztest (5, 0, 1)` gives 1, where the first gave the scalar 1 and the second
    a 1 by 0 empty; with `dim` given and a NaN present the call raised
    `nonconformant arguments`.  An empty input keeps its shape and a group with
    no data gives NaN, where every empty gave the scalar 0.

 ** grpstats:
    an empty X keeps its columns and has no groups, where every empty gave a
    0 by 0 whatever the statistic asked for: `grpstats (zeros (0,3), zeros
    (0,1))` gives a 0 by 3, `'meanci'` a 0 by 3 by 2 and `'gname'` an empty
    cell column, as MATLAB returns.  A grouping variable whose values are all
    missing now leaves no groups instead of raising `'group_mean' undefined`,
    and a row missing in any of several grouping variables is dropped rather
    than indexing the group names with NaN.  `'gname'` with no grouping
    variable returns `{'1'}` where it raised `'g_names' undefined`.

 ** grpstats:
    `'sem'` on a matrix divides each column by its own count instead of
    solving a least squares system, which collapsed every column onto one
    value: `grpstats ([1 10; 2 20; 3 30], [1 1 1], 'sem')` gives
    `[0.5774 5.7735]` where it gave `[3.1754 3.1754]`.  Columns of equal
    spread hid it.  The table route was already right.

 ** grpstats:
    `'VarNames'` renames the variables of the output table instead of raising
    `table.renamevars: too few input arguments.`

 ** grpstats:
    a function handle in WHICHSTATS works on a matrix, where it raised
    `binary operator '>' not implemented for 'function handle' by 'scalar'`
    before any statistic was computed.  It is applied to each column of each
    group on the values that are not NaN, as MATLAB applies it, and a table
    names the new variable after the function rather than after `@mean`; the
    table route was passing the NaNs through, so `@mean` answered NaN for a
    group holding one and `@numel` counted them.

 ** dwtest, runstest:
    an empty input gives a result instead of raising.  `dwtest` returns `p = 0`
    and `d = NaN`, `runstest` returns `h = 0`, `pval = 1` and a `stats` of
    `nruns` and `z` both NaN with `n1` and `n0` both 0, as MATLAB returns.
    `dwtest (zeros (0,1), zeros (0,2))` died with `A(0,_): subscripts must be
    either integers 1 to (2^63)-1 or logicals`, a zero-length vector clearing
    the vector check and reaching the body; the other empty shapes raised a
    validation error.  Two deliberate differences from MATLAB: `d` is NaN for
    every empty input, where MATLAB gives 0 for some of them and NaN for
    others, and a character array is refused, where MATLAB accepts one and
    returns the empty result.

 ** ClassificationKNN, fitcknn:
    `predict` on a model fitted to a categorical response raised
    `cat: categorical/vertcat method failed` instead of returning labels.  The
    labels now come back as a categorical array, in the type of the response,
    as they already did for every other response type.

 ** fitcdiscr, fitcgam, fitcknn, fitcnb, fitcsvm:
    a `Prior` given as a structure was refused for a response held in a
    character matrix or a categorical array, the first with "'ClassNames' in
    the 'Prior' structure must name every class of the model" and the second
    with an error from `categorical.eq`.  The structure's class names are now
    matched by text for both, and may themselves be a character matrix.

 ** crossval, kfoldPredict:
    a classifier fitted to a categorical or string response could not be
    cross-validated: `crossval` passed the model's own class names back to the
    learner, whose constructor accepted 'ClassNames' only as a cell array of
    character vectors, a logical or numeric vector or a character array.  All
    eight classifier constructors now accept categorical and string
    'ClassNames', and `kfoldPredict` returns labels in the type of the
    response, missing ones as `<undefined>` or `<missing>`.

 ** savemodel, loadmodel (classification models):
    a classifier fitted to a categorical or string response did not survive a
    save.  `savemodel` warned that classdef objects cannot be saved and wrote
    its class names and response as structures, after which the SVM, linear
    and kernel classifiers loaded with `ClassNames` a structure and predicted
    structures, and the KNN, naive Bayes, discriminant, GAM and neural network
    classifiers refused the file.  Both types now save without a warning and
    load back as they were, and files saved by 1.9.2 load too.

 ** fitcsvm, fitcgam, fitckernel, fitclinear, fitcecoc:
    a categorical response carrying a category with no observations, such as a
    subset of a categorical array, counted that category as a class.  The four
    binary learners refused the fit and `fitcecoc` built a coding matrix for a
    class with no data.  Unused categories are no longer classes, as in the
    other classifiers.

 ** fitcnb, ClassificationNaiveBayes:
    'ClassNames' given as a cell array of character vectors was refused for a
    categorical or string response with "no observation belongs to the named
    classes".  Class names and labels of any textual type are now matched by
    their text.

 ** grp2idx:
    a categorical grouping variable was numbered over the values it held,
    sorted by name, while the names returned listed every category, so the
    numbers and the names disagreed whenever a category was unused or the
    categories were not in alphabetical order, and the levels lost their
    category order and ordinality.  It is now numbered by category code, over
    every category and in category order, as MATLAB does.

 ** anova:
    the fitted values of a one-way model grouped by a categorical array with an
    unused category, or with categories out of alphabetical order, were taken
    from the wrong group's mean.  Each observation now takes the mean of its
    own group.

 ** grpstats:
    grouped by a categorical array with an unused category, or with categories
    out of alphabetical order, the statistics were reported against the wrong
    group names.  Each group now carries its own statistics, one per category
    in category order, as MATLAB reports.

 ** crosstab:
    the rows and columns of a categorical variable are laid out in category
    order, as MATLAB lays them out, where they were sorted by name.

 Improvements:
 =============

 ** fitcdiscr, ClassificationDiscriminant:
    `'CategoricalPredictors'` is refused with a message saying that
    discriminant analysis cannot use categorical predictors, instead of as an
    invalid parameter name.

 ** @seealso, package wide:
    a class member named in an `@seealso` is now written qualified,
    `Class.member`, where 122 of them were bare.  `help margin` finds nothing
    while `help ClassificationTree.margin` prints the method, so a bare name
    told the reader to type something that does not work; the online pages
    link a qualified member to its own page or to the class page and the
    collapsible holding it.  A reference to a function keeps its bare name.

 ** ClassificationPartitionedModel:
    accepts a `ClassificationTree`, reporting `CrossValidatedModel` as `Tree`
    and storing a `CompactClassificationTree` per fold.  `Cost` cannot be
    assigned on one, as it cannot on a cross-validated `ClassificationSVM`,
    the cost being consumed while the tree is fitted.

 ** RegressionPartitionedModel:
    accepts a `RegressionTree`, reporting `CrossValidatedModel` as `Tree` and
    storing a `CompactRegressionTree` per fold.

 ** fitlme, fitlmematrix, fitglme:
    the fits no longer form a matrix whose size is the number of observations.
    At 12800 observations and 20 random effects `fitlme` takes 0.10 seconds and
    0.17 GiB where it took 114 seconds and 3.8 GiB; at 6400 observations
    `fitglme` takes 0.15 seconds and 0.18 GiB where it took 134 seconds and
    1.4 GiB.  Peak memory no longer grows with the number of observations, so
    fits above about 25000 observations are possible where they previously were
    not.

 ** cmdscale:
    points that are all coincident give a single column of zeros rather than an
    empty matrix, whose shape varied with how the input was given: a 4 by 0 for
    `cmdscale (zeros (4,4))` and a 0 by 2 for `cmdscale (0)`.  With `p` given,
    the single column is kept and `e` keeps its length of `p`.

 ** grpstats:
    `grpstats (tbl)` summarises a table with no grouping variable, every
    variable being a data variable and the single group named `All`, as
    MATLAB does; it raised `table.removevars: too few input arguments.`

 ** grpstats:
    `[means, sem, counts, gname] = grpstats (x, group)` returns the four
    default statistics without naming any, as MATLAB does; up to four outputs
    were an error unless WHICHSTATS listed exactly that many.

 ** dwtest:
    a design matrix of full rank leaves no residual degrees of freedom and the
    statistic has no null distribution.  The p-value was 0.5 by the exact
    method, from an empty coefficient set, and the approximate method reported
    0 for an undefined value.  It is now 0, 1 or 0 by tail for the exact method
    and NaN for the approximate one.

 ** regress:
    a row vector `y` is taken as a column, as in `robustfit`, `lasso`,
    `fitlm`, `stepwiselm` and `stepwisefit`; it was refused.

 ** fitcgam, fitrgam, ClassificationGAM, RegressionGAM:
    accept `'CategoricalPredictors'` under the boosted-tree engine, as indices,
    a logical vector or `'all'`.  A tree splits a categorical predictor into two
    sets of levels, a level not seen in training predicts as a missing value,
    and its `BinEdges` and `PairDetectionBinEdges` are empty, as in MATLAB.
    `gamboosttrain`, `gamboostinter` and `gamboostpairs` take a matching
    logical argument flagging columns of level codes.

 ** fitcgam, fitrgam, ClassificationGAM, RegressionGAM:
    under the boosted-tree engine an interaction tree that splits on only one
    of its pair's predictors adds nothing, and a pair left with no tree is
    dropped with a warning, as in MATLAB.  Such trees were kept, so
    `'Interactions', 'all'` fits moved away from MATLAB's after the first
    round.

 ** fitclinear, fitrlinear, fitckernel, fitrkernel:
    the 'bfgs' solver keeps the full inverse Hessian and both quasi-Newton
    solvers take MATLAB's weak Wolfe line search, so fits reach MATLAB's
    iterations and coefficients; they differed in the fourth significant digit.
    `ModelParameters.LineSearch` reports 'weakwolfe'.  The coefficient
    tolerance is tested before the gradient tolerance, as in MATLAB.

 ** fitclinear, fitrlinear (lasso, 'sparsa'):
    the solver tests the gradient tolerance, a tolerance of 0 switches its test
    off, and a line search that finds no decrease reports code -11.  With both
    tolerances at 0 it reported the coefficient tolerance satisfied.

 ** fitcdiscr, ClassificationDiscriminant:
    accept `'Weights'`.  The class means and covariances are weighted, an
    empirical prior sums the weights per class, a row of zero weight is left
    out, and `W` holds the weights over their sum, as in MATLAB.
    Cross-validation folds keep the weights of their rows.

 ** fitcknn, ClassificationKNN:
    accept `'Weights'`.  Each class carries its prior spread over its rows in
    proportion to their weights, and a neighbour votes with that weight; an
    empirical prior sums the weights, standardization is weighted, and a row
    of zero weight is left out, as in MATLAB.  Scores under a prior other than
    the empirical one now follow it: a neighbour always voted as one.

 ** knnsearch, rangesearch, KDTreeSearcher, fitcknn:
    the kd-tree search returns neighbours at equal distance in row order and
    keeps the earlier row as the K-th neighbour, as MATLAB does.  It returned
    them in the order it met them and could return a later row instead, so a
    KNN classifier could vote with a different neighbour.

 ** isoutlier:
    the `'percentiles'` method works on a matrix.  It took both thresholds from
    one `quantile` call and indexed the result linearly, which only holds for a
    vector, so every column was tested against the first column's pair and the
    outlier mask was wrong without saying so.  An option that needs a value now
    says which one is missing; `'movmean'`, `'movmedian'`, `'percentiles'`,
    `'SamplePoints'`, `'ThresholdFactor'` and `'MaxNumOutliers'` all failed on
    an internal subscript error instead.

 ** isoutlier:
    the `'quartiles'` method works.  It passed the operating dimension to
    `quantile` where the probability belongs, so every call raised
    `Q(3): out of bound 2`.  `'ThresholdFactor'` applies to it, and the help
    text no longer says otherwise.

 ** dummyvar:
    a `NaN` in a numeric grouping variable is missing data, filling that
    variable's own block of columns with `NaN` and leaving the others alone, as
    MATLAB does; it raised.  A column that is entirely `NaN` has no levels, so
    it produces no columns and a warning names it.

 ** anova1, kruskalwallis, vartestn, levene_test, bartlett_test:
    a vector `X` without `GROUP` is a single group, whichever way it lies, as
    MATLAB reads it; all five raised instead.  With one group there is nothing
    to compare, so `p` is `NaN`; `kruskalwallis` and `vartestn` returned 0.
    An empty `X` gives `anova1` and `kruskalwallis` an answer with no degrees
    of freedom below zero, where MATLAB reports negative ones, and the three
    variance tests refuse it by name, as MATLAB's `vartestn` does.

 ** knnsearch:
    `'K'` above the number of observations returns all of them, as MATLAB does;
    it raised a nonconformant error on `'kdtree'` and named a range on
    `'exhaustive'`.  This also settles `'K', Inf`, which ran out of memory, and
    a single-row `X`, which returned that one point K times.  An empty `X` or
    `Y` returns empty rather than raising.

 ** LIBSVM:
    the bundled `svm.cpp` and `svm.h` move to the instance-weights variant
    3.37.  Nothing in the two files changed from 3.36 but the version
    constant.

 ** fitcsvm, ClassificationSVM, svmtrain:
    accept observation weights.  `svmtrain` takes a fourth argument of
    instance weights, each scaling its instance's box constraint, from LIBSVM's
    instance-weights variant.  `fitcsvm` gives each observation the box
    constraint n * C * w, w its weight scaled to the cost-adjusted prior; an
    empirical prior sums the weights, standardization is weighted, and a row
    of zero weight is left out, as in MATLAB.  Prior and Cost now enter the fit
    the same way with or without weights, and `loss` scales the weights of
    each class to its prior, as MATLAB does; it counted every row alike.

 ** fitcsvm, ClassificationSVM, CompactClassificationSVM (one-class):
    `predict` labels every row with the model's one class and returns one
    score column, as MATLAB does; it failed with an index error whenever a
    row was an outlier.  A one-class fit ignores observation weights, and its
    classification error is 0, as in MATLAB.
    Its bias puts the least supported support vector on the boundary, or with
    `'OutlierFraction'` leaves that fraction of the training rows below it
    while keeping `Nu`, as MATLAB does; it was LIBSVM's average over the free
    support vectors, and `'OutlierFraction'` refitted the model with that `Nu`.

 ** crosstab:
    a categorical group gives an unused category its zero row or column, as
    MATLAB does, so the table lines up with `labels`; the table dropped the
    row while `labels` kept its name, so every later row was mislabelled.
    The chi-square statistic and p-value are computed over the rows and
    columns holding observations, where MATLAB returns NaN when any is empty.

 ** Classifiers with a categorical response:
    `ClassNames`, and the labels `predict` returns, carry only the categories
    the response holds, as MATLAB reports them; an unused category stayed in
    their category list.

 ** ClassificationKNN, ClassificationDiscriminant,
    CompactClassificationDiscriminant:
    `loss`, `margin` and `edge` match the true labels to the classes by name,
    whatever their type, as every other classifier does and as MATLAB does;
    cellstr, string or categorical labels of another type than `ClassNames`
    were refused.

 ** ClassificationSVM, CompactClassificationSVM (one-class):
    `margin`, and so `edge` and `resubEdge`, return NaN, as MATLAB does; they
    returned twice the score.  `resubLoss` with `'classifcost'` returns 0, as
    MATLAB does; it failed with an index error when a training row was an
    outlier.

 ** fitlme, fitlmematrix, fitglme, LinearMixedModel,
    GeneralizedLinearMixedModel:
    the random-effects design is held as a sparse matrix, and the fit factorises
    its covariance in sparse form unless crossed terms make it fill in, so large
    fits with many groups use far less memory and run much faster;
    `designMatrix` returns it sparse for `'Random'`, as MATLAB does.

 ** LinearMixedModel:
    `residuals` with `'ResidualType', 'Standardized'` no longer builds
    n-by-n matrices; it took 35 s and 7.7 GB at 12000 observations and could
    not run on large data sets, and now takes milliseconds.  `anova` with
    `'DFMethod', 'Satterthwaite'` runs 10 to 250 times faster with many groups.

 Summary of important user-visible changes for statistics 1.9.4:
-------------------------------------------------------------------

 Important Notice: Observation weights given to the machine learning models
                   must be single or double; logical and integer weights are
                   refused.

 Backwards incompatible improvements:
 ====================================

 ** Observation weights:
    `'Weights'` must be single or double at fitting and in `loss`, `edge`,
    `resubLoss`, `resume` and the `TreeBagger` error methods, as in MATLAB.
    Integer weights fitted a different model or failed; logical weights
    worked only in `fitcecoc`, `fitcensemble` and `fitrensemble`.

 New functions and methods:
 ==========================

 ** boxchart, stats.chart.BoxChart:
    box charts from vectors or a table, with notches, either orientation,
    jittered outliers and colour groups.  The returned object redraws when
    its properties are set.

 ** gumbelstat:
    the mean and variance of the Gumbel distribution.

 ** lime:
    local interpretable model-agnostic explanations for any model that
    predicts and for a function handle, with a linear or a tree simple model,
    and a `plot` method.  `'Distance'` also takes a function handle, which
    MATLAB does not.

 ** partialDependence, plotPartialDependence:
    partial dependence of a fitted model on one or two predictors, drawn as a
    line, a surface or individual conditional curves.

 ** PredictiveModel:
    an abstract superclass of every model with a `predict` method, so that
    `isa (Mdl, 'PredictiveModel')` tells whether a model predicts.

 ** shapley:
    Shapley values for any model that predicts and for a function handle,
    exact or estimated, interventional or conditional, and exact at any size
    for linear models and trees.  `plot`, `boxchart`, `swarmchart` and
    `plotDependence` draw them.

 ** swarmchart:
    swarm charts.  `XData` holds where each point was drawn, where MATLAB
    leaves it unspread.

 ** Table input:
    every learner takes a table at fitting, and so do `predict`, `loss`,
    `margin`, `edge` and the other methods that take observations, matched
    to the predictors by name.  A variable holding levels is a categorical
    predictor.  `shapley` and `lime` take tables too.

 Bug fixes:
 ==========

 ** anova1:
    the table for `'unequal'` variances carried a sum of squares column
    headed `F`; it now has the five columns the function prints.

 ** anova2:
    an empty `X` failed on an internal variable; it is refused with a clear
    error, where MATLAB returns NaN p-values.

 ** burrinv:
    every quantile of a Burr distribution with a scale other than 1 was
    wrong, which reached `icdf` and `median` of `BurrDistribution`.

 ** ClassificationDiscriminant:
    `mahal` with `'ClassLabels'` failed for a categorical response and
    refused categorical and string labels.

 ** ClassificationDiscriminant, ClassificationKNN:
    `resubLoss` failed when given both `'LossFun'` and `'Weights'`.

 ** ClassificationGAM, RegressionGAM:
    with `'FitMethod', 'splines'`, `'Knots'` and `'DoF'` given together kept
    only the second.

 ** ClassificationKNN:
    `partialDependence` failed for class labels given as text.

 ** ClassificationSVM:
    `'Nu'` on two classes trained a one-class model, `'OutlierFraction'`
    overwrote a given `'Nu'`, and `loss` and `margin` failed for a
    categorical or string response.

 ** ClassificationTree, ClassificationNaiveBayes, ClassificationLinear,
    ClassificationKernel:
    `edge` on data missing a class was scaled by the prior of the classes
    present.

 ** combnk:
    `combnk (DATA, 0)` returned 0-by-0 instead of the 1-by-0 empty
    combination.

 ** datasample:
    NaN and all-zero weights were sampled from, and an `Inf` weight drew the
    last element, repeatedly even without replacement.

 ** dummyvar:
    a grouping variable holding `Inf` ran out of memory, and an empty one
    with several columns failed.

 ** ff2n:
    `ff2n (0)` was refused; it returns 1-by-0.  The error message printed
    texinfo markup.

 ** fitcecoc:
    weights varying within a class were refused for SVM, KNN, naive Bayes and
    discriminant learners.

 ** fitcknn, fitcsvm:
    `'Leaveout', 'on'` was refused.

 ** fitcnet:
    a prior other than the empirical one did not affect training, and `W`
    was misaligned when a row missed a predictor.

 ** fitcsvm, fitcknn, fitcnet:
    standardization scaled a constant predictor by about 1e16.

 ** fitrgp:
    `W` held the weights as given instead of summing to one.

 ** gpcdf, gppdf, gpinv, gprnd, gpfit:
    the help text gave the wrong condition for the Pareto special case.

 ** grp2idx:
    empty grouping variables returned wrongly shaped outputs or failed, and a
    character array without rows counted as one observation.

 ** iforest:
    fewer than 3 observations were refused with an error about an option the
    caller had not given.

 ** libsvmwrite:
    column indices were written in the wrong format on Windows.

 ** LinearModel:
    `anova` reported a rounding residue, or an infinite F, for a term that
    adds no rank, and failed on a model holding only an intercept.

 ** mvncdf, mvtcdf:
    an options structure missing fields, as `statset` makes, was refused.

 ** randsample:
    `randsample (Inf, 1)` returned a huge number and an infinite `k` ran out
    of memory.

 ** ranksum, vartestn:
    the help text described a rule and fields the functions do not have.

 ** RegressionSVM:
    `IsSupportVector` marked the wrong observations when a row missed a
    predictor.

 ** rmmissing:
    empty input accepted invalid options and returned a wrongly shaped `TF`;
    `'MinNumMissing'` may be 0, and on a vector it counts per element.

 ** Single observation weights:
    `fitclinear` and `fitrlinear` fitted a different model from single
    weights.  `W` keeps the class of the weights; `Prior` and the predictions
    stay double where MATLAB returns single.

 ** standardizeMissing:
    an empty or matrix indicator was refused.

 ** svmtrain, svmpredict:
    an options string of 2048 characters or more was read past its end.

 Improvements:
 =============

 ** Observation weights:
    `fitcnb`, `fitrsvm`, `fitcnet` and `fitrnet` take `'Weights'`, and their
    `resubLoss` and cross-validation weigh by them.

 ** CategoricalPredictors:
    every learner that codes categorical predictors takes them by name as
    well as by index.

 ** TreeBagger, CompactTreeBagger:
    a read-only `CategoricalPredictors` property.

 ** pdist:
    inputs of 1000 observations or more are computed much faster, 0.02 s
    instead of 46 s for `'correlation'` on 1200 rows.

 ** Machine_Learning:
    `inst/Supervised_Learning` is renamed `inst/Machine_Learning`, and its
    `INDEX` group "Machine Learning".

 ** Compiler warnings:
    the compiled sources build without warnings under `-Wall -Wextra`.
