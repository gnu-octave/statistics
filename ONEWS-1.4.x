Summary of important user-visible changes for statistics 1.4.3:
-------------------------------------------------------------------

 New functions:
 ==============

 ** anova1            (patch #10127)
    kruskalwallis

 ** cluster            (patch #10009)

 ** clusterdata        (patch #10012)

 ** confusionchart     (patch #9985)

 ** confusionmat       (patch #9971)

 ** cophenet           (patch #10040)

 ** datasample         (patch #10050)

 ** evalclusters       (patch #10052)

 ** expfit             (patch #10092)
    explike

 ** gscatter           (patch #10043)

 ** ismissing          (patch #10102)

 ** inconsistent       (patch #10008)

 ** mhsample.m         (patch #10016)

 ** ncx2pdf            (patch #9711)

 ** optimalleaforder.m (patch #10034)

 ** pca                (patch #10104)

 ** rmmissing          (patch #10102)

 ** silhouette         (patch #9743)

 ** slicesample        (patch #10019)

 ** wblplot            (patch #8579)

 Improvements:
 =============

 ** anovan.m: use double instead of toascii (bug #60514)

 ** binocdf: new option "upper" (bug #43721)

 ** boxplot: better Matlab compatibility; several Matlab-compatible
    plot options added (OutlierTags, Sample_IDs, BoxWidth, Widths,
    BoxStyle, Positions, Labels, Colors) and an Octave-specific one
    (CapWidhts); demos added; texinfo improved (patch #9930)

 ** auto MPG (carbig) sample dataset added from
    https://archive.ics.uci.edu/ml/datasets/Auto+MPG (patch #10045)

 ** crosstab.m: make n-dimensional (patch #10014)

 ** dendrogram.m: many improvements (patch #10036)

 ** fitgmdist.m: fix typo in ComponentProportion (bug #59386)

 ** gevfit: change orientation of results for Matlab compatibility (bug #47369)

 ** hygepdf: avoid overflow for certain inputs (bug #35827)

 ** kmeans: efficiency and compatibility tweaks (patch #10042)

 ** pdist: option for squared Euclidean distance (patch #10051)

 ** stepwisefit.m: give another option to select predictors (patch #8584)

 ** tricdf, triinv: fixes (bug #60113)


Summary of important user-visible changes for statistics 1.4.2:
-------------------------------------------------------------------

 ** canoncorr: allow more variables than observations

 ** fitgmdist: return fitgmdist parameters (Bug #57917)

 ** gamfit: invert parameter per docs (Bug #57849)

 ** geoXXX: update docs 'number of failures (X-1)' => 'number of failures (X)' (Bug #57606)

 ** kolmogorov_smirnov_test.m: update function handle usage from octave6+ (Bug #57351)

 ** linkage.m: fix octave6+ parse error (Bug #57348)

 ** unifrnd: changed unifrnd(a,a) to return a 0 rather than NaN (Bug #56342)

 ** updates for usage of depreciated octave functions

Summary of important user-visible changes for statistics 1.4.1:
-------------------------------------------------------------------
 ** update install scripts for octave 5.0 depreciated functions

 ** bug fixes to the following functions:
      pdist2.m: use max in distEucSq (Bug #50377)
      normpdf: use eps tolerance in tests (Bug #51963)
      fitgmdist: fix an output bug in fitgmdist
      t_test: Set tolerance on t_test BISTS (Bug #54557)
      gpXXXXX: change order of inputs to match matlab (Bug #54009)
      bartlett_test: df = k-1 (Bug #45894)
      gppdf: apply scale factor (Bug #54009)
      gmdistribution: updates for bug #54278, ##54279
      wishrnd: Bug #55860

Summary of important user-visible changes for statistics 1.4.0:
-------------------------------------------------------------------

 ** The following functions are new:

      canoncorr
      fitgmdist
      gmdistribution
      sigma_pts

 ** The following functions have been moved from the statistics package but are
    conditionally installed:

      mad

 ** The following functions have been moved from octave to be conditionally
    installed:

    BASE
      cloglog
      logit
      prctile
      probit
      qqplot
      table  (renamed to crosstab)

    DISTRIBUTIONS
      betacdf
      betainv
      betapdf
      betarnd
      binocdf
      binoinv
      binopdf
      binornd
      cauchy_cdf
      cauchy_inv
      cauchy_pdf
      cauchy_rnd
      chi2cdf
      chi2inv
      chi2pdf
      chi2rnd
      expcdf
      expinv
      exppdf
      exprnd
      fcdf
      finv
      fpdf
      frnd
      gamcdf
      gaminv
      gampdf
      gamrnd
      geocdf
      geoinv
      geopdf
      geornd
      hygecdf
      hygeinv
      hygepdf
      hygernd
      kolmogorov_smirnov_cdf
      laplace_cdf
      laplace_inv
      laplace_pdf
      laplace_rnd
      logistic_cdf
      logistic_inv
      logistic_pdf
      logistic_rnd
      logncdf
      logninv
      lognpdf
      lognrnd
      nbincdf
      nbininv
      nbinpdf
      nbinrnd
      normcdf
      norminv
      normpdf
      normrnd
      poisscdf
      poissinv
      poisspdf
      poissrnd
      stdnormal_cdf
      stdnormal_inv
      stdnormal_pdf
      stdnormal_rnd
      tcdf
      tinv
      tpdf
      trnd
      unidcdf
      unidinv
      unidpdf
      unidrnd
      unifcdf
      unifinv
      unifpdf
      unifrnd
      wblcdf
      wblinv
      wblpdf
      wblrnd
      wienrnd

    MODELS
      logistic_regression

    TESTS
      anova
      bartlett_test
      chisquare_test_homogeneity
      chisquare_test_independence
      cor_test
      f_test_regression
      hotelling_test
      hotelling_test_2
      kolmogorov_smirnov_test
      kolmogorov_smirnov_test_2
      kruskal_wallis_test
      manova
      mcnemar_test
      prop_test_2
      run_test
      sign_test
      t_test
      t_test_2
      t_test_regression
      u_test
      var_test
      welch_test
      wilcoxon_test
      z_test
      z_test_2

 ** Functions marked with known test failures:
      grp2idx: bug #51928
      gevfir_lmom: bug #31070

 ** Other functions that have been changed for smaller bugfixes, increased
    Matlab compatibility, or performance:

      dcov: returned dcov instead of dcor. added demo.
      violin: can be used with subplots. violin quality improved.
      princomp: Fix expected values of tsquare in unit tests
      fitgmdist: test number inputs to function
      hist3: fix removal of rows with NaN values

 ** added the packages test data to install
