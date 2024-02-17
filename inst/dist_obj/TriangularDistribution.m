classdef TriangularDistribution
  ## -*- texinfo -*-
  ## @deftypefn {statistics} TriangularDistribution
  ##
  ## Triangular probability distribution object.
  ##
  ## A @code{TriangularDistribution} object consists of parameters, a model
  ## description, and sample data for a triangular probability distribution.
  ##
  ## The triangular distribution uses the following parameters.
  ##
  ## @multitable @columnfractions 0.25 0.48 0.27
  ## @headitem @var{Parameter} @tab @var{Description} @tab @var{Support}
  ##
  ## @item @qcode{A} @tab Lower limit @tab @math{-Inf < A < Inf}
  ## @item @qcode{B} @tab Peak location @tab @math{A <= B <= C}
  ## @item @qcode{C} @tab Upper limit @tab @math{C > A}
  ## @end multitable
  ##
  ## There are several ways to create a @code{TriangularDistribution} object.
  ##
  ## @itemize
  ## @item Create a distribution with specified parameter values using the
  ## @code{makedist} function.
  ## @item Use the constructor @qcode{TriangularDistribution (@var{A}, @var{B}
  ## @var{C})} to create a triangular distribution with specified parameter
  ## values.
  ## @end itemize
  ##
  ## It is highly recommended to use @code{makedist} function to create
  ## probability distribution objects, instead of the constructor.
  ##
  ## A @code{TriangularDistribution} object contains the following
  ## properties, which can be accessed using dot notation.
  ##
  ## @multitable @columnfractions 0.25 0.25 0.25 0.25
  ## @item @qcode{DistributionName} @tab @qcode{DistributionCode} @tab
  ## @qcode{NumParameters} @tab @qcode{ParameterNames}
  ## @item @qcode{ParameterDescription} @tab @qcode{ParameterValues} @tab
  ## @qcode{Truncation} @tab @qcode{IsTruncated}
  ## @end multitable
  ##
  ## A @code{TriangularDistribution} object contains the following methods:
  ## @code{cdf}, @code{icdf}, @code{iqr}, @code{mean}, @code{median},
  ## @code{pdf}, @code{plot}, @code{random}, @code{std}, @code{truncate},
  ## @code{var}.
  ##
  ## Further information about the continuous triangular distribution can be found
  ## at @url{https://en.wikipedia.org/wiki/Continuous_uniform_distribution}
  ##
  ## @seealso{makedist, tricdf, triinv, tripdf, trirnd, tristat}
  ## @end deftypefn

  properties (Dependent = true)
    A
    B
    C
  endproperties

  properties (GetAccess = public, Constant = true)
    DistributionName = "TriangularDistribution";
    DistributionCode = "tri";
    NumParameters = 3;
    ParameterNames = {"A", "B", "C"};
    ParameterDescription = {"Lower limit", "Peak location", "Upper limit"};
  endproperties

  properties (GetAccess = public , SetAccess = protected)
    ParameterValues
    Truncation
    IsTruncated
    InputData = [];
  endproperties

  methods (Hidden)

    function this = TriangularDistribution (A, B, C)
      if (nargin == 0)
        A = 0;
        B = 0.5;
        C = 1;
      endif
      checkparams (A, B, C)
      this.IsTruncated = false;
      this.ParameterValues = [A, B, C];
    endfunction

    function display (this)
      fprintf ("%s =\n", inputname(1));
      __disp__ (this, "Triangular distribution");
    endfunction

    function disp (this)
      __disp__ (this, "Triangular distribution");
    endfunction

    function this = set.A (this, A)
      checkparams (A, this.B, this.C)
      this.ParameterValues(1) = A;
    endfunction

    function A = get.A (this)
      A = this.ParameterValues(1);
    endfunction

    function this = set.B (this, B)
      checkparams (this.A, B, this.C)
      this.ParameterValues(2) = B;
    endfunction

    function B = get.B (this)
      B = this.ParameterValues(2);
    endfunction

    function this = set.C (this, C)
      checkparams (this.A, this.B, C)
      this.ParameterValues(3) = C;
    endfunction

    function C = get.C (this)
      C = this.ParameterValues(3);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {TriangularDistribution} {@var{p} =} cdf (@var{pd}, @var{x})
    ## @deftypefnx {TriangularDistribution} {@var{p} =} cdf (@var{pd}, @var{x}, @qcode{"upper"})
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
      p = tricdf (x, this.A, this.B, this.C);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        p(lb) = 0;
        p(ub) = 1;
        p(! (lb | ub)) -= tricdf (lx, this.A, this.B, this.C);
        p(! (lb | ub)) /= diff (tricdf ([lx, ux], this.A, this.B, this.C));
      endif
      ## Apply uflag
      if (utail)
        p = 1 - p;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TriangularDistribution} {@var{p} =} icdf (@var{pd}, @var{p})
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
        lp = tricdf (this.Truncation(1), this.A, this.B, this.C);
        up = tricdf (this.Truncation(2), this.A, this.B, this.C);
        ## Adjust p values within range of p @ lower limit and p @ upper limit
        np = lp + (up - lp) .* p;
        x = triinv (np, this.A, this.B, this.C);
      else
        x = triinv (p, this.A, this.B, this.C);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TriangularDistribution} {@var{r} =} iqr (@var{pd})
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
    ## @deftypefn  {TriangularDistribution} {@var{m} =} mean (@var{pd})
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
        m = tristat (this.A, this.B, this.C);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TriangularDistribution} {@var{m} =} median (@var{pd})
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
        Fa_b = tricdf ([lx, ux], this.A, this.B, this.C);
        m = triinv (sum (Fa_b) / 2, this.A, this.B, this.C);
      else
        m = tristat (this.A, this.B, this.C);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TriangularDistribution} {@var{y} =} pdf (@var{pd}, @var{x})
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
      y = tripdf (x, this.A, this.B, this.C);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        y(lb | ub) = 0;
        y(! (lb | ub)) /= diff (tricdf ([lx, ux], this.A, this.B, this.C));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TriangularDistribution} {} plot (@var{pd})
    ## @deftypefnx {TriangularDistribution} {} plot (@var{pd}, @var{Name}, @var{Value})
    ## @deftypefnx {TriangularDistribution} {@var{h} =} plot (@dots{})
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
    ## @deftypefn  {TriangularDistribution} {@var{y} =} random (@var{pd})
    ## @deftypefnx {TriangularDistribution} {@var{y} =} random (@var{pd}, @var{rows})
    ## @deftypefnx {TriangularDistribution} {@var{y} =} random (@var{pd}, @var{rows}, @var{cols}, @dots{})
    ## @deftypefnx {TriangularDistribution} {@var{y} =} random (@var{pd}, [@var{sz}])
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
        sz = [varargin{:}];
        ps = prod (sz);
        ## Get an estimate of how many more random numbers we need to randomly
        ## pick the appropriate size from
        lx = this.Truncation(1);
        ux = this.Truncation(2);
        ratio = 1 / diff (tricdf ([ux, lx], this.A, this.B, this.C));
        nsize = 2 * ratio * ps;     # times 2 to be on the safe side
        ## Generate the numbers and remove out-of-bound random samples
        r = trirnd (this.A, this.B, this.C, nsize, 1);
        r(r < lx | r > ux) = [];
        ## Randomly select the required size and reshape to requested dimensions
        r = randperm (r, ps);
        r = reshape (r, sz);
      else
        r = trirnd (this.A, this.B, this.C, varargin{:});
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TriangularDistribution} {@var{s} =} std (@var{pd})
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
      v = var (this.A, this.B, this.C);
      s = sqrt (v);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TriangularDistribution} {@var{t} =} truncate (@var{pd}, @var{lower}, @var{upper})
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
      ## Check boundaries and constrain within support [A, C]
      lower(lower < this.A) = this.A;
      upper(upper > this.C) = this.C;
      this.Truncation = [lower, upper];
      this.IsTruncated = true;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TriangularDistribution} {@var{v} =} var (@var{pd})
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
        [~, v] = tristat (this.A, this.B, this.C);
      endif
    endfunction

  endmethods

endclassdef

function checkparams (A, B, C)
  if (! (isscalar (A) && isnumeric (A) && isreal (A) && isfinite (A)))
    error ("TriangularDistribution: lower limit A must be a real scalar.")
  endif
  if (! (isscalar (B) && isnumeric (B) && isreal (B) && isfinite (B)))
    error ("TriangularDistribution: mode B must be a real scalar.")
  endif
  if (! (isscalar (C) && isnumeric (C) && isreal (C) && isfinite (C)))
    error ("TriangularDistribution: upper limit C must be a real scalar.")
  endif
  if (! (A < C))
    error (strcat (["TriangularDistribution: lower limit A must"], ...
                   [" be less than upper limit C."]))
  endif
  if (! (A <= B && B <= C))
    error (strcat (["TriangularDistribution: mode B must be within"], ...
                   [" lower limit A and upper limit C."]))
  endif
endfunction

## Test input validation
## 'TriangularDistribution' constructor
%!error <TriangularDistribution: lower limit A must be a real scalar.> ...
%! TriangularDistribution (i, 1, 2)
%!error <TriangularDistribution: lower limit A must be a real scalar.> ...
%! TriangularDistribution (Inf, 1, 2)
%!error <TriangularDistribution: lower limit A must be a real scalar.> ...
%! TriangularDistribution ([1, 2], 1, 2)
%!error <TriangularDistribution: lower limit A must be a real scalar.> ...
%! TriangularDistribution ("a", 1, 2)
%!error <TriangularDistribution: lower limit A must be a real scalar.> ...
%! TriangularDistribution (NaN, 1, 2)
%!error <TriangularDistribution: mode B must be a real scalar.> ...
%! TriangularDistribution (1, i, 2)
%!error <TriangularDistribution: mode B must be a real scalar.> ...
%! TriangularDistribution (1, Inf, 2)
%!error <TriangularDistribution: mode B must be a real scalar.> ...
%! TriangularDistribution (1, [1, 2], 2)
%!error <TriangularDistribution: mode B must be a real scalar.> ...
%! TriangularDistribution (1, "a", 2)
%!error <TriangularDistribution: mode B must be a real scalar.> ...
%! TriangularDistribution (1, NaN, 2)
%!error <TriangularDistribution: upper limit C must be a real scalar.> ...
%! TriangularDistribution (1, 2, i)
%!error <TriangularDistribution: upper limit C must be a real scalar.> ...
%! TriangularDistribution (1, 2, Inf)
%!error <TriangularDistribution: upper limit C must be a real scalar.> ...
%! TriangularDistribution (1, 2, [1, 2])
%!error <TriangularDistribution: upper limit C must be a real scalar.> ...
%! TriangularDistribution (1, 2, "a")
%!error <TriangularDistribution: upper limit C must be a real scalar.> ...
%! TriangularDistribution (1, 2, NaN)
%!error <TriangularDistribution: lower limit A must be less than upper limit C.> ...
%! TriangularDistribution (1, 1, 1)
%!error <TriangularDistribution: mode B must be within lower limit A and upper limit C.> ...
%! TriangularDistribution (1, 0.5, 2)

## 'cdf' method
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (TriangularDistribution, 2, "uper")
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (TriangularDistribution, 2, 3)

## 'plot' method
%!error <plot: optional arguments must be in NAME-VALUE pairs.> ...
%! plot (TriangularDistribution, "Parent")
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (TriangularDistribution, "PlotType", 12)
%!error <plot: invalid VALUE size for 'Parameter' argument.> ...
%! plot (TriangularDistribution, "PlotType", {"pdf", "cdf"})
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (TriangularDistribution, "PlotType", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (TriangularDistribution, "Discrete", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (TriangularDistribution, "Discrete", [1, 0])
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (TriangularDistribution, "Discrete", {true})
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (TriangularDistribution, "Parent", 12)
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (TriangularDistribution, "Parent", "hax")

## 'truncate' method
%!error <truncate: missing input argument.> ...
%! truncate (TriangularDistribution)
%!error <truncate: missing input argument.> ...
%! truncate (TriangularDistribution, 2)
%!error <truncate: invalid lower upper limits.> ...
%! truncate (TriangularDistribution, 4, 2)

## Catch errors when using array of probability objects with available methods
%!shared pd
%! pd = TriangularDistribution (0, 1, 2);
%! pd(2) = TriangularDistribution (0, 1, 2);
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
