classdef UniformDistribution
  ## -*- texinfo -*-
  ## @deftypefn {statistics} UniformDistribution
  ##
  ## Continuous uniform probability distribution object.
  ##
  ## A @code{UniformDistribution} object consists of parameters, a model
  ## description, and sample data for a uniform probability distribution.
  ##
  ## The uniform distribution uses the following parameters.
  ##
  ## @multitable @columnfractions 0.25 0.48 0.27
  ## @headitem @var{Parameter} @tab @var{Description} @tab @var{Support}
  ##
  ## @item @qcode{Lower} @tab Lower limit @tab @math{-Inf < Lower < Upper}
  ## @item @qcode{Upper} @tab Upper limit @tab @math{Lower < Upper < Inf}
  ## @end multitable
  ##
  ## There are several ways to create a @code{UniformDistribution} object.
  ##
  ## @itemize
  ## @item Create a distribution with specified parameter values using the
  ## @code{makedist} function.
  ## @item Use the constructor @qcode{UniformDistribution (@var{Lower},
  ## @var{Upper})} to create a uniform distribution with specified parameter
  ## values.
  ## @end itemize
  ##
  ## It is highly recommended to use @code{makedist} function to create
  ## probability distribution objects, instead of the constructor.
  ##
  ## A @code{UniformDistribution} object contains the following
  ## properties, which can be accessed using dot notation.
  ##
  ## @multitable @columnfractions 0.25 0.25 0.25 0.25
  ## @item @qcode{DistributionName} @tab @qcode{DistributionCode} @tab
  ## @qcode{NumParameters} @tab @qcode{ParameterNames}
  ## @item @qcode{ParameterDescription} @tab @qcode{ParameterValues} @tab
  ## @qcode{Truncation} @tab @qcode{IsTruncated}
  ## @end multitable
  ##
  ## A @code{UniformDistribution} object contains the following methods:
  ## @code{cdf}, @code{icdf}, @code{iqr}, @code{mean}, @code{median},
  ## @code{pdf}, @code{plot}, @code{random}, @code{std}, @code{truncate},
  ## @code{var}.
  ##
  ## Further information about the continuous uniform distribution can be found
  ## at @url{https://en.wikipedia.org/wiki/Continuous_uniform_distribution}
  ##
  ## @seealso{makedist, unifcdf, unifinv, unifpdf, unifrnd, unifit, unifstat}
  ## @end deftypefn

  properties (Dependent = true)
    Lower
    Upper
  endproperties

  properties (GetAccess = public, Constant = true)
    DistributionName = "UniformDistribution";
    DistributionCode = "unif";
    NumParameters = 2;
    ParameterNames = {"Lower", "Upper"};
    ParameterDescription = {"Lower limit", "Upper limit"};
  endproperties

  properties (GetAccess = public , SetAccess = protected)
    ParameterValues
    Truncation
    IsTruncated
    InputData = [];
  endproperties

  methods (Hidden)

    function this = UniformDistribution (Lower, Upper)
      if (nargin == 0)
        Lower = 0;
        Upper = 1;
      endif
      checkparams (Lower, Upper)
      this.IsTruncated = false;
      this.ParameterValues = [Lower, Upper];
    endfunction

    function display (this)
      fprintf ("%s =\n", inputname(1));
      __disp__ (this, "Uniform distribution (continuous)");
    endfunction

    function disp (this)
      __disp__ (this, "Uniform distribution (continuous)");
    endfunction

    function this = set.Lower (this, Lower)
      checkparams (Lower, this.Upper)
      this.ParameterValues(1) = Lower;
    endfunction

    function Lower = get.Lower (this)
      Lower = this.ParameterValues(1);
    endfunction

    function this = set.Upper (this, Upper)
      checkparams (this.Lower, Upper)
      this.ParameterValues(2) = Upper;
    endfunction

    function Upper = get.Upper (this)
      Upper = this.ParameterValues(2);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {UniformDistribution} {@var{p} =} cdf (@var{pd}, @var{x})
    ## @deftypefnx {UniformDistribution} {@var{p} =} cdf (@var{pd}, @var{x}, @qcode{"upper"})
    ##
    ## Compute the cumulative distribution function (CDF).
    ##
    ## @code{@var{p} = cdf (@var{pd}, @var{x})} computes the CDF of the
    ## probability distribution object, @var{pd}, evaluated at the values in
    ## @var{x}.
    ##
    ## @code{@var{p} = cdf (@dots{}, @qcode{"upper"})} returns the complement of
    ## the CDF of the probability distribution object, @var{pd}, evaluated at
    ## the values in @var{x}.
    ##
    ## @end deftypefn
    function p = cdf (this, x, uflag)
      if (! isscalar (this))
        error ("cdf: requires a scalar probability distribution.");
      endif
      ## Check for "upper" flag
      if (nargin > 2 && strcmpi (uflag, "upper"))
        utail = true;
      elseif (nargin > 2 && ! strcmpi (uflag, "upper"))
        error ("cdf: invalid argument for upper tail.");
      else
        utail = false;
      endif
      ## Do the computations
      p = unifcdf (x, this.Lower, this.Upper);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        p(lb) = 0;
        p(ub) = 1;
        p(! (lb | ub)) -= unifcdf (lx, this.Lower, this.Upper);
        p(! (lb | ub)) /= diff (unifcdf ([lx, ux], this.Lower, this.Upper));
      endif
      ## Apply uflag
      if (utail)
        p = 1 - p;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {UniformDistribution} {@var{p} =} icdf (@var{pd}, @var{p})
    ##
    ## Compute the cumulative distribution function (CDF).
    ##
    ## @code{@var{p} = icdf (@var{pd}, @var{x})} computes the quantile (the
    ## inverse of the CDF) of the probability distribution object, @var{pd},
    ## evaluated at the values in @var{x}.
    ##
    ## @end deftypefn
    function x = icdf (this, p)
      if (! isscalar (this))
        error ("icdf: requires a scalar probability distribution.");
      endif
      if (this.IsTruncated)
        lp = unifcdf (this.Truncation(1), this.Lower, this.Upper);
        up = unifcdf (this.Truncation(2), this.Lower, this.Upper);
        ## Adjust p values within range of p @ lower limit and p @ upper limit
        np = lp + (up - lp) .* p;
        x = unifinv (np, this.Lower, this.Upper);
      else
        x = unifinv (p, this.Lower, this.Upper);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {UniformDistribution} {@var{r} =} iqr (@var{pd})
    ##
    ## Compute the interquartile range of a probability distribution.
    ##
    ## @code{@var{r} = iqr (@var{pd})} computes the interquartile range of the
    ## probability distribution object, @var{pd}.
    ##
    ## @end deftypefn
    function r = iqr (this)
      if (! isscalar (this))
        error ("iqr: requires a scalar probability distribution.");
      endif
        r = diff (icdf (this, [0.25, 0.75]));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {UniformDistribution} {@var{m} =} mean (@var{pd})
    ##
    ## Compute the mean of a probability distribution.
    ##
    ## @code{@var{m} = mean (@var{pd})} computes the mean of the probability
    ## distribution object, @var{pd}.
    ##
    ## @end deftypefn
    function m = mean (this)
      if (! isscalar (this))
        error ("mean: requires a scalar probability distribution.");
      endif
      if (this.IsTruncated)
        fm = @(x) x .* pdf (this, x);
        m = integral (fm, this.Truncation(1), this.Truncation(2));
      else
        m = unifstat (this.Lower, this.Upper);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {UniformDistribution} {@var{m} =} median (@var{pd})
    ##
    ## Compute the median of a probability distribution.
    ##
    ## @code{@var{m} = median (@var{pd})} computes the median of the probability
    ## distribution object, @var{pd}.
    ##
    ## @end deftypefn
    function m = median (this)
      if (! isscalar (this))
        error ("median: requires a scalar probability distribution.");
      endif
      if (this.IsTruncated)
        lx = this.Truncation(1);
        ux = this.Truncation(2);
        Fa_b = unifcdf ([lx, ux], this.Lower, this.Upper);
        m = unifinv (sum (Fa_b) / 2, this.Lower, this.Upper);
      else
        m = unifstat (this.Lower, this.Upper);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {UniformDistribution} {@var{y} =} pdf (@var{pd}, @var{x})
    ##
    ## Compute the probability distribution function (PDF).
    ##
    ## @code{@var{y} = pdf (@var{pd}, @var{x})} computes the PDF of the
    ## probability distribution object, @var{pd}, evaluated at the values in
    ## @var{x}.
    ##
    ## @end deftypefn
    function y = pdf (this, x)
      if (! isscalar (this))
        error ("pdf: requires a scalar probability distribution.");
      endif
      y = unifpdf (x, this.Lower, this.Upper);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        y(lb | ub) = 0;
        y(! (lb | ub)) /= diff (unifcdf ([lx, ux], this.Lower, this.Upper));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {UniformDistribution} {} plot (@var{pd})
    ## @deftypefnx {UniformDistribution} {} plot (@var{pd}, @var{Name}, @var{Value})
    ## @deftypefnx {UniformDistribution} {@var{h} =} plot (@dots{})
    ##
    ## Plot a probability distribution object.
    ##
    ## @code{plot (@var{pd}} plots a probability density function (PDF) of the
    ## probability distribution object @var{pd}.  If @var{pd} contains data,
    ## which have been fitted by @code{fitdist}, the PDF is superimposed over a
    ## histogram of the data.
    ##
    ## @code{plot (@var{pd}, @var{Name}, @var{Value})} specifies additional
    ## options with the @qcode{Name-Value} pair arguments listed below.
    ##
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @tab @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{"PlotType"} @tab @tab A character vector specifying the plot
    ## type.  @qcode{"pdf"} plots the probability density function (PDF).  When
    ## @var{pd} is fit to data, the PDF is superimposed on a histogram of the
    ## data.  @qcode{"cdf"} plots the cumulative density function (CDF).  When
    ## @var{pd} is fit to data, the CDF is superimposed over an empirical CDF.
    ## @qcode{"probability"} plots a probability plot using a CDF of the data
    ## and a CDF of the fitted probability distribution.  This option is
    ## available only when @var{pd} is fitted to data.
    ##
    ## @item @qcode{"Discrete"} @tab @tab A logical scalar to specify whether to
    ## plot the PDF or CDF of a discrete distribution object as a line plot or a
    ## stem plot, by specifying @qcode{false} or @qcode{true}, respectively.  By
    ## default, it is @qcode{true} for discrete distributions and @qcode{false}
    ## for continuous distributions.  When @var{pd} is a continuous distribution
    ## object, option is ignored.
    ##
    ## @item @qcode{"Parent"} @tab @tab An axes graphics object for plot.  If
    ## not specified, the @code{plot} function plots into the current axes or
    ## creates a new axes object if one does not exist.
    ## @end multitable
    ##
    ## @code{@var{h} = plot (@dots{})} returns a graphics handle to the plotted
    ## objects.
    ##
    ## @end deftypefn
    function [varargout] = plot (this, varargin)
      if (! isscalar (this))
        error ("plot: requires a scalar probability distribution.");
      endif
      h = __plot__ (this, false, varargin{:});
      if (nargout > 0)
        varargout{1} = h;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {UniformDistribution} {@var{y} =} random (@var{pd})
    ## @deftypefnx {UniformDistribution} {@var{y} =} random (@var{pd}, @var{rows})
    ## @deftypefnx {UniformDistribution} {@var{y} =} random (@var{pd}, @var{rows}, @var{cols}, @dots{})
    ## @deftypefnx {UniformDistribution} {@var{y} =} random (@var{pd}, [@var{sz}])
    ##
    ## Generate random arrays from the probability distribution object.
    ##
    ## @code{@var{r} = random (@var{pd})} returns a random number from the
    ## distribution object @var{pd}.
    ##
    ## When called with a single size argument, @code{betarnd} returns a square
    ## matrix with the dimension specified.  When called with more than one
    ## scalar argument, the first two arguments are taken as the number of rows
    ## and columns and any further arguments specify additional matrix
    ## dimensions.  The size may also be specified with a row vector of
    ## dimensions, @var{sz}.
    ##
    ## @end deftypefn
    function r = random (this, varargin)
      if (! isscalar (this))
        error ("random: requires a scalar probability distribution.");
      endif
      if (this.IsTruncated)
        r = unifrnd (this.Truncation(1), this.Truncation(2), varargin{:});
      else
        r = unifrnd (this.Lower, this.Upper, varargin{:});
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {UniformDistribution} {@var{s} =} std (@var{pd})
    ##
    ## Compute the standard deviation of a probability distribution.
    ##
    ## @code{@var{s} = std (@var{pd})} computes the standard deviation of the
    ## probability distribution object, @var{pd}.
    ##
    ## @end deftypefn
    function s = std (this)
      if (! isscalar (this))
        error ("std: requires a scalar probability distribution.");
      endif
      v = var (this.Lower, this.Upper);
      s = sqrt (v);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {UniformDistribution} {@var{t} =} truncate (@var{pd}, @var{lower}, @var{upper})
    ##
    ## Truncate a probability distribution.
    ##
    ## @code{@var{t} = truncate (@var{pd})} returns a probability distribution
    ## @var{t}, which is the probability distribution @var{pd} truncated to the
    ## specified interval with lower limit, @var{lower}, and upper limit,
    ## @var{upper}.  If @var{pd} is fitted to data with @code{fitdist}, the
    ## returned probability distribution @var{t} is not fitted, does not contain
    ## any data or estimated values, and it is as it has been created with the
    ## @var{makedist} function, but it includes the truncation interval.
    ##
    ## @end deftypefn
    function this = truncate (this, lower, upper)
      if (! isscalar (this))
        error ("truncate: requires a scalar probability distribution.");
      endif
      if (nargin < 3)
        error ("truncate: missing input argument.");
      elseif (lower >= upper)
        error ("truncate: invalid lower upper limits.");
      endif
      ## Check boundaries and constrain within support [Lower, Upper]
      lower(lower < this.Lower) = this.Lower;
      upper(upper > this.Upper) = this.Upper;
      this.Truncation = [lower, upper];
      this.IsTruncated = true;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {UniformDistribution} {@var{v} =} var (@var{pd})
    ##
    ## Compute the variance of a probability distribution.
    ##
    ## @code{@var{v} = var (@var{pd})} computes the standard deviation of the
    ## probability distribution object, @var{pd}.
    ##
    ## @end deftypefn
    function v = var (this)
      if (! isscalar (this))
        error ("var: requires a scalar probability distribution.");
      endif
      if (this.IsTruncated)
        fm = @(x) x .* pdf (this, x);
        m = integral (fm, this.Truncation(1), this.Truncation(2));
        fv =  @(x) ((x - m) .^ 2) .* pdf (pd, x);
        v = integral (fv, this.Truncation(1), this.Truncation(2));
      else
        [~, v] = unifstat (this.Lower, this.Upper);
      endif
    endfunction

  endmethods

endclassdef

function checkparams (Lower, Upper)
  if (! (isscalar (Lower) && isnumeric (Lower) && isreal (Lower)
                          && isfinite (Lower)))
    error ("UniformDistribution: LOWER must be a real scalar.")
  endif
  if (! (isscalar (Upper) && isnumeric (Upper) && isreal (Upper)
                          && isfinite (Upper)))
    error ("UniformDistribution: UPPER must be a real scalar.")
  endif
  if (! (Lower < Upper))
    error ("UniformDistribution: LOWER must be less than UPPER.")
  endif
endfunction

## Test input validation
## 'UniformDistribution' constructor
%!error <UniformDistribution: LOWER must be a real scalar.> ...
%! UniformDistribution (i, 1)
%!error <UniformDistribution: LOWER must be a real scalar.> ...
%! UniformDistribution (Inf, 1)
%!error <UniformDistribution: LOWER must be a real scalar.> ...
%! UniformDistribution ([1, 2], 1)
%!error <UniformDistribution: LOWER must be a real scalar.> ...
%! UniformDistribution ("a", 1)
%!error <UniformDistribution: LOWER must be a real scalar.> ...
%! UniformDistribution (NaN, 1)
%!error <UniformDistribution: UPPER must be a real scalar.> ...
%! UniformDistribution (1, i)
%!error <UniformDistribution: UPPER must be a real scalar.> ...
%! UniformDistribution (1, Inf)
%!error <UniformDistribution: UPPER must be a real scalar.> ...
%! UniformDistribution (1, [1, 2])
%!error <UniformDistribution: UPPER must be a real scalar.> ...
%! UniformDistribution (1, "a")
%!error <UniformDistribution: UPPER must be a real scalar.> ...
%! UniformDistribution (1, NaN)
%!error <UniformDistribution: LOWER must be less than UPPER.> ...
%! UniformDistribution (2, 1)

## 'cdf' method
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (UniformDistribution, 2, "uper")
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (UniformDistribution, 2, 3)

## 'plot' method
%!error <plot: optional arguments must be in NAME-VALUE pairs.> ...
%! plot (UniformDistribution, "Parent")
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (UniformDistribution, "PlotType", 12)
%!error <plot: invalid VALUE size for 'Parameter' argument.> ...
%! plot (UniformDistribution, "PlotType", {"pdf", "cdf"})
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (UniformDistribution, "PlotType", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (UniformDistribution, "Discrete", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (UniformDistribution, "Discrete", [1, 0])
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (UniformDistribution, "Discrete", {true})
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (UniformDistribution, "Parent", 12)
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (UniformDistribution, "Parent", "hax")

## 'truncate' method
%!error <truncate: missing input argument.> ...
%! truncate (UniformDistribution)
%!error <truncate: missing input argument.> ...
%! truncate (UniformDistribution, 2)
%!error <truncate: invalid lower upper limits.> ...
%! truncate (UniformDistribution, 4, 2)

## Catch errors when using array of probability objects with available methods
%!shared pd
%! pd = UniformDistribution (0, 1);
%! pd(2) = UniformDistribution (0, 2);
%!error <cdf: requires a scalar probability distribution.> cdf (pd, 1)
%!error <icdf: requires a scalar probability distribution.> icdf (pd, 0.5)
%!error <iqr: requires a scalar probability distribution.> iqr (pd)
%!error <mean: requires a scalar probability distribution.> mean (pd)
%!error <median: requires a scalar probability distribution.> median (pd)
%!error <pdf: requires a scalar probability distribution.> pdf (pd, 1)
%!error <plot: requires a scalar probability distribution.> plot (pd)
%!error <random: requires a scalar probability distribution.> random (pd)
%!error <std: requires a scalar probability distribution.> std (pd)
%!error <truncate: requires a scalar probability distribution.> ...
%! truncate (pd, 2, 4)
%!error <var: requires a scalar probability distribution.> var (pd)
