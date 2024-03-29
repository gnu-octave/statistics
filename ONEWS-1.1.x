Summary of important user-visible changes for statistics 1.1.3:
-------------------------------------------------------------------

 ** The following functions are new in 1.1.3:

      copularnd mvtrnd

 ** The functions mnpdf and mnrnd are now also usable for greater numbers
    of categories for which the rows do not exactly sum to 1.

Summary of important user-visible changes for statistics 1.1.2:
-------------------------------------------------------------------

 ** The following functions are new in 1.1.2:

      mnpdf mnrnd

 ** The package is now dependent on the io package (version 1.0.18 or
    later) since the functions that it depended of from miscellaneous
    package have been moved to io.

 ** The function `kmeans' now accepts the 'emptyaction' property with
    the 'singleton' value. This allows for the kmeans algorithm to handle
    empty cluster better. It also throws an error if the user does not
    request an empty cluster handling, and there is an empty cluster.
    Plus, the returned items are now a closer match to Matlab.

Summary of important user-visible changes for statistics 1.1.1:
-------------------------------------------------------------------

 ** The following functions are new in 1.1.1:

      monotone_smooth kmeans jackknife

 ** Bug fixes on the functions:

      normalise_distribution  combnk
      repanova

 ** The following functions were removed since equivalents are now
    part of GNU octave core:

      zscore

 ** boxplot.m now returns a structure with handles to the plot elemenets.

Summary of important user-visible changes for statistics 1.1.0:
-------------------------------------------------------------------

 ** IMPORTANT note about `fstat' shadowing core library function:

    GNU octave's 3.2 release added a new function `fstat' to return
    information of a file. Statistics' `fstat' computes F mean and
    variance. Since MatLab's `fstat' is the equivalent to statistics'
    `fstat' (not to core's `fstat'), and to avoid problems with the
    statistics package, `fstat' has been deprecated in octave 3.4
    and will be removed in Octave 3.8. In the mean time, please
    ignore this warning when installing the package.

 ** The following functions are new in 1.1.0:

      normalise_distribution  repanova  combnk

 ** The following functions were removed since equivalents are now
    part of GNU octave core:

      prctile

 ** The __tbl_delim__ function is now private.

 ** The function `boxplot' now accepts named arguments.

 ** Bug fixes on the functions:

      harmmean  nanmax  nanmin  regress

 ** Small improvements on help text.
