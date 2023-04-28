Summary of important user-visible changes for statistics 1.2.4:
-------------------------------------------------------------------

 ** Made princomp work with nargout < 2.

 ** Renamed dendogram to dendrogram.

 ** Added isempty check to kmeans.

 ** Transposed output of hist3.

 ** Converted calculation in hmmviterbi to log space.

 ** Bug fixes for
    stepwisefit   wishrnd.

 ** Rewrite of cmdscale for improved compatibility.

 ** Fix in squareform for improved compatibility.

 ** New cvpartition class, with methods:

    display repartition test training

 ** New sample data file fisheriris.txt for tests

 ** The following functions are new:

    cdf crossval dcov pdist2 qrandn randsample signtest ttest ttest2
    vartest vartest2 ztest


Summary of important user-visible changes for statistics 1.2.3:
-------------------------------------------------------------------

 ** Made sure that output of nanstd is real.

 ** Fixed second output of nanmax and nanmin.

 ** Corrected handle for outliers in boxplot.

 ** Bug fix and enhanced functionality for mvnrnd.

 ** The following functions are new:

    wishrnd iwishrnd wishpdf iwishpdf cmdscale

Summary of important user-visible changes for statistics 1.2.2:
-------------------------------------------------------------------

 ** Fixed documentation of dendogram and hist3 to work with TexInfo 5.

Summary of important user-visible changes for statistics 1.2.1:
-------------------------------------------------------------------

 ** The following functions are new:

      pcares  pcacov  runstest  stepwisefit hist3

 ** dendogram now returns the leaf node numbers and order that the nodes were displayed in.

 ** New faster implementation of princomp.

Summary of important user-visible changes for statistics 1.2.0:
-------------------------------------------------------------------

 ** The following functions are new:

      regress_gp  dendogram   plsregress

 ** New functions for the generalized extreme value (GEV) distribution:

      gevcdf gevfit gevfit_lmom gevinv gevlike gevpdf gevrnd gevstat

 ** The interface of the following functions has been modified:

      mvnrnd

 ** `kmeans' has been fixed to deal with clusters that contain only
    one element.

 ** `normplot' has been fixed to avoid use of functions that have been
    removed from Octave core. Also, the plot produced should now display some
    aesthetic elements and appropriate legends.

 ** The help text of `mvtrnd' has been improved.

 ** Package is no longer autoloaded.
