## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

classdef LoguniformDistribution
  ## -*- texinfo -*-
  ## @deftypefn {statistics} LoguniformDistribution
  ##
  ## Log-uniform probability distribution object.
  ##
  ## A @code{LoguniformDistribution} object consists of parameters, a model
  ## description, and sample data for a log-uniform probability distribution.
  ##
  ## The log-uniform distribution uses the following parameters.
  ##
  ## @multitable @columnfractions 0.25 0.48 0.27
  ## @headitem @var{Parameter} @tab @var{Description} @tab @var{Support}
  ##
  ## @item @qcode{Lower} @tab Lower limit @tab @math{0 < Lower < Upper}
  ## @item @qcode{Upper} @tab Upper limit @tab @math{Lower < Upper < Inf}
  ## @end multitable
  ##
  ## There are several ways to create a @code{LoguniformDistribution} object.
  ##
  ## @itemize
  ## @item Create a distribution with specified parameter values using the
  ## @code{makedist} function.
  ## @item Use the constructor @qcode{LoguniformDistribution (@var{Lower})}
  ## to create a log-uniform distribution with specified parameter values.
  ## @end itemize
  ##
  ## It is highly recommended to use the @code{makedist} function to create
  ## probability distribution objects, instead of the constructor.
  ##
  ## A @code{LoguniformDistribution} object contains the following properties,
  ## which can be accessed using dot notation.
  ##
  ## @multitable @columnfractions 0.25 0.25 0.25 0.25
  ## @item @qcode{DistributionName} @tab @qcode{DistributionCode} @tab
  ## @qcode{NumParameters} @tab @qcode{ParameterNames}
  ## @item @qcode{ParameterDescription} @tab @qcode{ParameterValues} @tab
  ## @qcode{Truncation} @tab @qcode{IsTruncated}
  ## @end multitable
  ##
  ## A @code{LoguniformDistribution} object contains the following methods:
  ## @code{cdf}, @code{icdf}, @code{iqr}, @code{mean}, @code{median},
  ## @code{pdf}, @code{plot}, @code{random}, @code{std}, @code{truncate},
  ## @code{var}.
  ##
  ## Further information about the log-uniform distribution can be found at
  ## @url{https://en.wikipedia.org/wiki/Reciprocal_distribution}
  ##
  ## @seealso{fitdist, makedist}
  ## @end deftypefn

  properties (Dependent = true)
    Lower
    Upper
  endproperties

  properties (GetAccess = public, Constant = true)
    CensoringAllowed = false;
    DistributionName = "LoguniformDistribution";
    DistributionCode = "logu";
    NumParameters = 2;
    ParameterNames = {"Lower", "Upper"};
    ParameterDescription = {"Outcome probabilities"};
  endproperties

  properties (GetAccess = public , SetAccess = protected)
    ParameterValues
    Truncation
    IsTruncated
  endproperties

  methods (Hidden)

    function this = LoguniformDistribution (Lower, Upper)
      if (nargin == 0)
        Lower = 1;
        Upper = 4;
      endif
      checkparams (Lower, Upper);
      this.IsTruncated = false;
      this.ParameterValues = [Lower, Upper];
    endfunction

    function display (this)
      fprintf ("%s =\n", inputname(1));
      __disp__ (this, "loguniform distribution");
    endfunction

    function disp (this)
      __disp__ (this, "loguniform distribution");
    endfunction

    function this = set.Lower (this, Lower)
      checkparams (Lower, this.Upper);
      this.ParameterValues(1) = Lower;
    endfunction

    function Lower = get.Lower (this)
      Lower = this.ParameterValues(1);
    endfunction

    function this = set.Upper (this, Upper)
      checkparams (this.Lower, Upper);
      this.ParameterValues(2) = Upper;
    endfunction

    function Upper = get.Upper (this)
      Upper = this.ParameterValues(2);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {LoguniformDistribution} {@var{p} =} cdf (@var{pd}, @var{x})
    ## @deftypefnx {LoguniformDistribution} {@var{p} =} cdf (@var{pd}, @var{x}, @qcode{"upper"})
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
      d = log (this.Upper / this.Lower);
      p = log (x / this.Lower) / d;
      p(x<this.Lower) = 0;
      p(x>this.Upper) = 1;
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        p(lb) = 0;
        p(ub) = 1;
        p(! (lb | ub)) -= log (lx / this.Lower) / d;
        p(! (lb | ub)) /= diff (log ([lx, ux] ./ this.Lower) ./ d);
      endif
      ## Apply uflag
      if (utail)
        p = 1 - p;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LoguniformDistribution} {@var{p} =} icdf (@var{pd}, @var{p})
    ##
    ## Compute the inverse cumulative distribution function (iCDF).
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
      ## Do the computations
      is_nan = p < 0 | p > 1 | isnan (p);
      is_val = p >= 0 & p <= 1;
      if (this.IsTruncated)
        lx = this.Truncation(1);
        ux = this.Truncation(2);
        d = log (this.Upper / this.Lower);
        lp = log (lx / this.Lower) / d;
        up = log (ux / this.Lower) / d;
        ## Adjust p values within range of p @ lower limit and p @ upper limit
        p = lp + (up - lp) .* p;
        is_nan = p < lp | p > up | isnan (p);
        is_val = p >= lp & p <= up;
      endif
      x = p;
      x(is_nan) = NaN;
      x(is_val) = (this.Upper .^ p(is_val)) ./ (this.Lower .^ (p(is_val) - 1));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LoguniformDistribution} {@var{r} =} iqr (@var{pd})
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
    ## @deftypefn  {LoguniformDistribution} {@var{m} =} mean (@var{pd})
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
        m = (this.Upper - this.Lower) / log (this.Upper / this.Lower);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LoguniformDistribution} {@var{m} =} median (@var{pd})
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
      m = icdf (this, 0.5);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LoguniformDistribution} {@var{y} =} pdf (@var{pd}, @var{x})
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
      d = log (this.Upper / this.Lower);
      y = 1 ./ (x .* d);
      y(x < this.Lower | x > this.Upper) = 0;
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        y(lb | ub) = 0;
        y(! (lb | ub)) /= diff (log ([lx, ux] ./ this.Lower) ./ d);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LoguniformDistribution} {} plot (@var{pd})
    ## @deftypefnx {LoguniformDistribution} {} plot (@var{pd}, @var{Name}, @var{Value})
    ## @deftypefnx {LoguniformDistribution} {@var{h} =} plot (@dots{})
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
    ## @deftypefn  {LoguniformDistribution} {@var{y} =} random (@var{pd})
    ## @deftypefnx {LoguniformDistribution} {@var{y} =} random (@var{pd}, @var{rows})
    ## @deftypefnx {LoguniformDistribution} {@var{y} =} random (@var{pd}, @var{rows}, @var{cols}, @dots{})
    ## @deftypefnx {LoguniformDistribution} {@var{y} =} random (@var{pd}, [@var{sz}])
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
      sz = [varargin{:}];
      r = icdf (this, rand (sz));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LoguniformDistribution} {@var{s} =} std (@var{pd})
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
      v = var (this);
      s = sqrt (v);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LoguniformDistribution} {@var{t} =} truncate (@var{pd}, @var{lower}, @var{upper})
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
      endif
      if (lower >= upper)
        error ("truncate: invalid lower upper limits.");
      endif
      this.Truncation = [lower, upper];
      this.IsTruncated = true;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LoguniformDistribution} {@var{v} =} var (@var{pd})
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
        fv =  @(x) ((x - m) .^ 2) .* pdf (this, x);
        v = integral (fv, this.Truncation(1), this.Truncation(2));
      else
        a = this.Lower;
        b = this.Upper;
        l = log (b / a);
        v = (b ^ 2 - a ^2) / (2 * l) - ((b - a) / l) ^ 2;
      endif
    endfunction

  endmethods

endclassdef

function checkparams (Lower, Upper)
  if (! (isscalar (Lower) && isnumeric (Lower) && isreal (Lower)
                          && isfinite (Lower) && Lower > 0))
    error ("LoguniformDistribution: LOWER must be a positive real scalar.")
  endif
  if (! (isscalar (Upper) && isnumeric (Upper) && isreal (Upper)
                          && isfinite (Upper)))
    error ("LoguniformDistribution: UPPER must be a real scalar.")
  endif
  if (! (Lower < Upper))
    error ("LoguniformDistribution: LOWER must be less than UPPER.")
  endif
endfunction

## Test output
%!shared pd, t
%! pd = LoguniformDistribution (1, 4);
%! t = truncate (pd, 2, 4);
%!assert (cdf (pd, [0, 1, 2, 3, 4, 5]), [0, 0, 0.5, 0.7925, 1, 1], 1e-4);
%!assert (cdf (t, [0, 1, 2, 3, 4, 5]), [0, 0, 0, 0.5850, 1, 1], 1e-4);
%!assert (cdf (pd, [1.5, 2, 3, 4]), [0.2925, 0.5, 0.7925, 1], 1e-4);
%!assert (cdf (t, [1.5, 2, 3, 4]), [0, 0, 0.5850, 1], 1e-4);
%!assert (icdf (pd, [0:0.2:1]), [1, 1.3195, 1.7411, 2.2974, 3.0314, 4], 1e-4);
%!assert (icdf (t, [0:0.2:1]), [2, 2.2974, 2.6390, 3.0314, 3.4822, 4], 1e-4);
%!assert (icdf (pd, [-1, 0.4:0.2:1, NaN]), [NaN, 1.7411, 2.2974, 3.0314, 4, NaN], 1e-4);
%!assert (icdf (t, [-1, 0.4:0.2:1, NaN]), [NaN, 2.6390, 3.0314, 3.4822, 4, NaN], 1e-4);
%!assert (iqr (pd), 1.4142, 1e-4);
%!assert (iqr (t), 0.9852, 1e-4);
%!assert (mean (pd), 2.1640, 1e-4);
%!assert (mean (t), 2.8854, 1e-4);
%!assert (median (pd), 2);
%!assert (median (t), 2.8284, 1e-4);
%!assert (pdf (pd, [0, 1, 2, 3, 4, 5]), [0, 0.7213, 0.3607, 0.2404, 0.1803, 0], 1e-4);
%!assert (pdf (t, [0, 1, 2, 3, 4, 5]), [0, 0, 0.7213, 0.4809, 0.3607, 0], 1e-4);
%!assert (pdf (pd, [-1, 1, 2, 3, 4, NaN]), [0, 0.7213, 0.3607, 0.2404, 0.1803, NaN], 1e-4);
%!assert (pdf (t, [-1, 1, 2, 3, 4, NaN]), [0, 0, 0.7213, 0.4809, 0.3607, NaN], 1e-4);
%!assert (isequal (size (random (pd, 100, 50)), [100, 50]))
%!assert (any (random (pd, 1000, 1) < 1), false);
%!assert (any (random (pd, 1000, 1) > 4), false);
%!assert (any (random (t, 1000, 1) < 2), false);
%!assert (any (random (t, 1000, 1) > 4), false);
%!assert (std (pd), 0.8527, 1e-4);
%!assert (std (t), 0.5751, 1e-4);
%!assert (var (pd), 0.7270, 1e-4);
%!assert (var (t), 0.3307, 1e-4);

## Test input validation
## 'LoguniformDistribution' constructor
%!error <LoguniformDistribution: LOWER must be a positive real scalar.> ...
%! LoguniformDistribution (i, 1)
%!error <LoguniformDistribution: LOWER must be a positive real scalar.> ...
%! LoguniformDistribution (Inf, 1)
%!error <LoguniformDistribution: LOWER must be a positive real scalar.> ...
%! LoguniformDistribution ([1, 2], 1)
%!error <LoguniformDistribution: LOWER must be a positive real scalar.> ...
%! LoguniformDistribution ("a", 1)
%!error <LoguniformDistribution: LOWER must be a positive real scalar.> ...
%! LoguniformDistribution (NaN, 1)
%!error <LoguniformDistribution: UPPER must be a real scalar.> ...
%! LoguniformDistribution (1, i)
%!error <LoguniformDistribution: UPPER must be a real scalar.> ...
%! LoguniformDistribution (1, Inf)
%!error <LoguniformDistribution: UPPER must be a real scalar.> ...
%! LoguniformDistribution (1, [1, 2])
%!error <LoguniformDistribution: UPPER must be a real scalar.> ...
%! LoguniformDistribution (1, "a")
%!error <LoguniformDistribution: UPPER must be a real scalar.> ...
%! LoguniformDistribution (1, NaN)
%!error <LoguniformDistribution: LOWER must be less than UPPER.> ...
%! LoguniformDistribution (2, 1)

## 'cdf' method
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (LoguniformDistribution, 2, "uper")
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (LoguniformDistribution, 2, 3)

## 'plot' method
%!error <plot: optional arguments must be in NAME-VALUE pairs.> ...
%! plot (LoguniformDistribution, "Parent")
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (LoguniformDistribution, "PlotType", 12)
%!error <plot: invalid VALUE size for 'Parameter' argument.> ...
%! plot (LoguniformDistribution, "PlotType", {"pdf", "cdf"})
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (LoguniformDistribution, "PlotType", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (LoguniformDistribution, "Discrete", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (LoguniformDistribution, "Discrete", [1, 0])
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (LoguniformDistribution, "Discrete", {true})
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (LoguniformDistribution, "Parent", 12)
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (LoguniformDistribution, "Parent", "hax")
%!error <plot: invalid NAME for optional argument.> ...
%! plot (LoguniformDistribution, "invalidNAME", "pdf")
%!error <plot: 'probability' PlotType is not supported for 'LoguniformDistribution'.> ...
%! plot (LoguniformDistribution, "PlotType", "probability")

## 'truncate' method
%!error <truncate: missing input argument.> ...
%! truncate (LoguniformDistribution)
%!error <truncate: missing input argument.> ...
%! truncate (LoguniformDistribution, 2)
%!error <truncate: invalid lower upper limits.> ...
%! truncate (LoguniformDistribution, 4, 2)

## Catch errors when using array of probability objects with available methods
%!shared pd
%! pd = LoguniformDistribution(1, 4);
%! pd(2) = LoguniformDistribution(2, 5);
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
