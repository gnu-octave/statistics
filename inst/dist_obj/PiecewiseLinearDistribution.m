## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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

classdef PiecewiseLinearDistribution
  ## -*- texinfo -*-
  ## @deftp {statistics} PiecewiseLinearDistribution
  ##
  ## Piecewise linear probability distribution object.
  ##
  ## A @code{PiecewiseLinearDistribution} object consists of parameters, a model
  ## description, and sample data for a piecewise linear probability
  ## distribution.
  ##
  ## The piecewise linear distribution is a continuous probability distribution
  ## that is defined by a set of points where the cumulative distribution
  ## function (CDF) changes slope.  It is defined by a vector of @math{x} values
  ## and a corresponding vector of CDF values @var{Fx}.
  ##
  ## There are several ways to create a @code{PiecewiseLinearDistribution}
  ## object.
  ##
  ## @itemize
  ## @item Create a distribution with specified parameter values using the
  ## @code{makedist} function.
  ## @item Use the constructor @qcode{PiecewiseLinearDistribution (@var{x},
  ## @var{Fx})} to create a piecewise linear distribution with specified
  ## parameter values @var{x} and @var{Fx}.
  ## @end itemize
  ##
  ## It is highly recommended to use @code{fitdist} and @code{makedist}
  ## functions to create probability distribution objects, instead of the class
  ## constructor or the aforementioned static method.
  ##
  ## Further information about the piecewise linear distribution can be found at
  ## @url{https://en.wikipedia.org/wiki/Piecewise_linear_function}
  ##
  ## @seealso{makedist, plcdf, plinv, plpdf, plrnd, plstat}
  ## @end deftp

  properties (Dependent = true)
    ## -*- texinfo -*-
    ## @deftp {PiecewiseLinearDistribution} {property} x
    ##
    ## Vector of x values
    ##
    ## A numeric vector of @math{x} values at which the CDF changes slope.  You
    ## can access the @qcode{x} property using dot name assignment.
    ##
    ## @end deftp
    x
    ## -*- texinfo -*-
    ## @deftp {PiecewiseLinearDistribution} {property} Fx
    ##
    ## Vector of CDF values
    ##
    ## A numeric vector of CDF values that correspond to each value in @math{x}.
    ## You can access the @qcode{Fx} property using dot name assignment.
    ##
    ## @end deftp
    Fx
  endproperties

  properties (GetAccess = public, Constant = true)
    ## -*- texinfo -*-
    ## @deftp {PiecewiseLinearDistribution} {property} DistributionName
    ##
    ## Probability distribution name
    ##
    ## A character vector specifying the name of the probability distribution
    ## object.  This property is read-only.
    ##
    ## @end deftp
    DistributionName = "PiecewiseLinearDistribution";

    ## -*- texinfo -*-
    ## @deftp {PiecewiseLinearDistribution} {property} NumParameters
    ##
    ## Number of parameters
    ##
    ## A scalar integer value specifying the number of parameters characterizing
    ## the probability distribution.  This property is read-only.
    ##
    ## @end deftp
    NumParameters = 2;
    ## -*- texinfo -*-
    ## @deftp {PiecewiseLinearDistribution} {property} ParameterNames
    ##
    ## Names of parameters
    ##
    ## A @math{2x1} cell array of character vectors with each element containing
    ## the name of a distribution parameter.  This property is read-only.
    ##
    ## @end deftp
    ParameterNames = {"x", "Fx"};
    ## -*- texinfo -*-
    ## @deftp {PiecewiseLinearDistribution} {property} ParameterDescription
    ##
    ## Description of parameters
    ##
    ## A @math{2x1} cell array of character vectors with each element containing
    ## a short description of a distribution parameter.  This property is
    ## read-only.
    ##
    ## @end deftp
    ParameterDescription = {"x", "cdf = F(x)"};
  endproperties

  properties (GetAccess = public, Constant = true, Hidden)
    CensoringAllowed = false;
    DistributionCode = "pl";
    ParameterRange = [-Inf, Inf; -Inf, Inf; 0, 1; 0, 1];
    ParameterLogCI = [false, false, false, false];
  endproperties

  properties (GetAccess = public , SetAccess = protected)
    ## -*- texinfo -*-
    ## @deftp {PiecewiseLinearDistribution} {property} ParameterValues
    ##
    ## Distribution parameter values
    ##
    ## A @math{2x1} numeric vector containing the values of the distribution
    ## parameters.  This property is read-only. You can change the distribution
    ## parameters by assigning new values to the @qcode{x} and @qcode{Fx}
    ## properties.
    ##
    ## @end deftp
    ParameterValues
    ## -*- texinfo -*-
    ## @deftp {PiecewiseLinearDistribution} {property} Truncation
    ##
    ## Truncation interval
    ##
    ## A @math{1x2} numeric vector specifying the truncation interval for the
    ## probability distribution.  First element contains the lower boundary,
    ## second element contains the upper boundary.  This property is read-only.
    ## You can only truncate a probability distribution with the
    ## @qcode{truncate} method.
    ##
    ## @end deftp
    Truncation
    ## -*- texinfo -*-
    ## @deftp {PiecewiseLinearDistribution} {property} IsTruncated
    ##
    ## Flag for truncated probability distribution
    ##
    ## A logical scalar value specifying whether a probability distribution is
    ## truncated or not.  This property is read-only.
    ##
    ## @end deftp
    IsTruncated
  endproperties

  methods (Hidden)

    function this = PiecewiseLinearDistribution (x, Fx)
      if (nargin == 0)
        x = [0; 1];
        Fx = [0; 1];
      else
        x = x(:);
        Fx = Fx(:);
      endif
      checkparams (x, Fx);
      this.IsTruncated = false;
      this.ParameterValues = [x, Fx];
    endfunction

    function display (this)
      fprintf ("%s =\n", inputname(1));
      __disp__ (this, "Piecewise Linear distribution");
    endfunction

    function disp (this)
      __disp__ (this, "Piecewise Linear distribution");
    endfunction

    function this = set.x (this, x)
      checkparams (x, this.Fx);
      this.ParameterValues(:,1) = x;
    endfunction

    function x = get.x (this)
      x = this.ParameterValues(:,1);
    endfunction

    function this = set.Fx (this, Fx)
      checkparams (this.x, Fx);
      this.ParameterValues(:,2) = Fx;
    endfunction

    function Fx = get.Fx (this)
      Fx = this.ParameterValues(:,2);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {PiecewiseLinearDistribution} {@var{p} =} cdf (@var{pd}, @var{x})
    ## @deftypefnx {PiecewiseLinearDistribution} {@var{p} =} cdf (@var{pd}, @var{x}, @qcode{"upper"})
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
      p = plcdf (x, this.x, this.Fx);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        p(lb) = 0;
        p(ub) = 1;
        p(! (lb | ub)) -= plcdf (lx, this.x, this.Fx);
        p(! (lb | ub)) /= diff (plcdf ([lx, ux], this.x, this.Fx));
      endif
      ## Apply uflag
      if (utail)
        p = 1 - p;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {PiecewiseLinearDistribution} {@var{x} =} icdf (@var{pd}, @var{p})
    ##
    ## Compute the inverse cumulative distribution function (iCDF).
    ##
    ## @code{@var{x} = icdf (@var{pd}, @var{p})} computes the quantile (the
    ## inverse of the CDF) of the probability distribution object, @var{pd},
    ## evaluated at the values in @var{p}.
    ##
    ## @end deftypefn
    function x = icdf (this, p)
      if (! isscalar (this))
        error ("icdf: requires a scalar probability distribution.");
      endif
      if (this.IsTruncated)
        lp = plcdf (this.Truncation(1), this.x, this.Fx);
        up = plcdf (this.Truncation(2), this.x, this.Fx);
        ## Adjust p values within range of p @ lower limit and p @ upper limit
        is_nan = p < 0 | p > 1;
        p(is_nan) = NaN;
        np = lp + (up - lp) .* p;
        x = plinv (np, this.x, this.Fx);
        x(x < this.Truncation(1)) = this.Truncation(1);
        x(x > this.Truncation(2)) = this.Truncation(2);
      else
        x = plinv (p, this.x, this.Fx);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {PiecewiseLinearDistribution} {@var{r} =} iqr (@var{pd})
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
    ## @deftypefn  {PiecewiseLinearDistribution} {@var{m} =} mean (@var{pd})
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
        m = integral (fm, this.Truncation(1), this.Truncation(2), ...
                      "ArrayValued", 1);
      else
        m = plstat (this.x, this.Fx);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {PiecewiseLinearDistribution} {@var{m} =} median (@var{pd})
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
        Fa_b = plcdf ([lx, ux], this.x, this.Fx);
        m = plinv (sum (Fa_b) / 2, this.x, this.Fx);
      else
        m = plinv (0.5, this.x, this.Fx);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {PiecewiseLinearDistribution} {@var{y} =} pdf (@var{pd}, @var{x})
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
      y = plpdf (x, this.x, this.Fx);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        y(lb | ub) = 0;
        y(! (lb | ub)) /= diff (plcdf ([lx, ux], this.x, this.Fx));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {PiecewiseLinearDistribution} {} plot (@var{pd})
    ## @deftypefnx {PiecewiseLinearDistribution} {} plot (@var{pd}, @var{Name}, @var{Value})
    ## @deftypefnx {PiecewiseLinearDistribution} {@var{h} =} plot (@dots{})
    ##
    ## Plot a probability distribution object.
    ##
    ## @code{plot (@var{pd})} plots a probability density function (PDF) of the
    ## probability distribution object @var{pd}.  If @var{pd} contains data,
    ## which have been fitted by @code{fitdist}, the PDF is superimposed over a
    ## histogram of the data.
    ##
    ## @code{plot (@var{pd}, @var{Name}, @var{Value})} specifies additional
    ## options with the @qcode{Name-Value} pair arguments listed below.
    ##
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Name} @tab @tab @var{Value}
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
    ## @deftypefn  {PiecewiseLinearDistribution} {@var{r} =} random (@var{pd})
    ## @deftypefnx {PiecewiseLinearDistribution} {@var{r} =} random (@var{pd}, @var{rows})
    ## @deftypefnx {PiecewiseLinearDistribution} {@var{r} =} random (@var{pd}, @var{rows}, @var{cols}, @dots{})
    ## @deftypefnx {PiecewiseLinearDistribution} {@var{r} =} random (@var{pd}, [@var{sz}])
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
        x = this.x(:)';
        Fx = this.Fx(:)';
        lp = plcdf (this.Truncation(1), x, Fx);
        up = plcdf (this.Truncation(2), x, Fx);
        u = unifrnd (lp, up, varargin{:});
        r = zeros (size (u));
        [~, bin] = histc (u(:)', Fx);
        r0 = x(bin);
        dx = diff (x);
        dF = diff (Fx);
        dr = (u(:)' - Fx(bin)) .* dx(bin) ./ dF(bin);
        r(:) = r0 + dr;
      else
        r = plrnd (this.x, this.Fx, varargin{:});
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {PiecewiseLinearDistribution} {@var{s} =} std (@var{pd})
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
    ## @deftypefn  {PiecewiseLinearDistribution} {@var{t} =} truncate (@var{pd}, @var{lower}, @var{upper})
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
      this.Truncation = [lower, upper];
      this.IsTruncated = true;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {PiecewiseLinearDistribution} {@var{v} =} var (@var{pd})
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
        m = integral (fm, this.Truncation(1), this.Truncation(2), ...
                      "ArrayValued", 1);
        fv =  @(x) ((x - m) .^ 2) .* pdf (this, x);
        v = integral (fv, this.Truncation(1), this.Truncation(2), ...
                      "ArrayValued", 1);
      else
        [~, v] = plstat (this.x, this.Fx);
      endif
    endfunction

  endmethods

endclassdef

function checkparams (x, Fx)
  if (! (isvector (x) && isnumeric (x) && isreal (x) && isfinite (x)))
    error ("PiecewiseLinearDistribution: X must be a real vector.")
  endif
  if (! (isvector (Fx) && isnumeric (Fx) && isreal (Fx) && isfinite (Fx)))
    error ("PiecewiseLinearDistribution: Fx must be a real vector.")
  endif
  if (! isvector (x) || ! isvector (Fx) || ! isequal (size (x), size (Fx)))
    error (strcat (["PiecewiseLinearDistribution: X and FX must"], ...
                   [" be vectors of equal size."]));
  endif
  if (length (x) < 2 || length (Fx) < 2)
    error (strcat (["PiecewiseLinearDistribution: X and FX must"], ...
                   [" be at least two-elements long."]));
  endif
  if (any (Fx < 0) || any (Fx > 1))
    error (strcat (["PiecewiseLinearDistribution: FX must be"], ...
                   [" bounded in the range [0, 1]."]));
  endif
endfunction

## Test output
%!shared pd, t
%! load patients
%! [f, x] = ecdf (Weight);
%! f = f(1:5:end);
%! x = x(1:5:end);
%! pd = PiecewiseLinearDistribution (x, f);
%! t = truncate (pd, 130, 180);
%!assert (cdf (pd, [120, 130, 140, 150, 200]), [0.0767, 0.25, 0.4629, 0.5190, 0.9908], 1e-4);
%!assert (cdf (t, [120, 130, 140, 150, 200]), [0, 0, 0.4274, 0.5403, 1], 1e-4);
%!assert (cdf (pd, [100, 250, NaN]), [0, 1, NaN], 1e-4);
%!assert (cdf (t, [115, 290, NaN]), [0, 1, NaN], 1e-4);
%!assert (icdf (pd, [0:0.2:1]), [111, 127.5, 136.62, 169.67, 182.17, 202], 1e-2);
%!assert (icdf (t, [0:0.2:1]), [130, 134.15, 139.26, 162.5, 173.99, 180], 1e-2);
%!assert (icdf (pd, [-1, 0.4:0.2:1, NaN]), [NA, 136.62, 169.67, 182.17, 202, NA], 1e-2);
%!assert (icdf (t, [-1, 0.4:0.2:1, NaN]), [NA, 139.26, 162.5, 173.99, 180, NA], 1e-2);
%!assert (iqr (pd), 50.0833, 1e-4);
%!assert (iqr (t), 36.8077, 1e-4);
%!assert (mean (pd), 153.61, 1e-10);
%!assert (mean (t), 152.311, 1e-3);
%!assert (median (pd), 142, 1e-10);
%!assert (median (t), 141.9462, 1e-4);
%!assert (pdf (pd, [120, 130, 140, 150, 200]), [0.0133, 0.0240, 0.0186, 0.0024, 0.0004], 1e-4);
%!assert (pdf (t, [120, 130, 140, 150, 200]), [0, 0.0482, 0.0373, 0.0048, 0], 1e-4);
%!assert (pdf (pd, [100, 250, NaN]), [0, 0, NaN], 1e-4);
%!assert (pdf (t, [100, 250, NaN]), [0, 0, NaN], 1e-4);
%!assert (isequal (size (random (pd, 100, 50)), [100, 50]))
%!assert (any (random (t, 1000, 1) < 130), false);
%!assert (any (random (t, 1000, 1) > 180), false);
%!assert (std (pd), 26.5196, 1e-4);
%!assert (std (t), 18.2941, 1e-4);
%!assert (var (pd), 703.2879, 1e-4);
%!assert (var (t), 334.6757, 1e-4);

## Test input validation
## 'PiecewiseLinearDistribution' constructor
%!error <PiecewiseLinearDistribution: X must be a real vector.> ...
%! PiecewiseLinearDistribution ([0, i], [0, 1])
%!error <PiecewiseLinearDistribution: X must be a real vector.> ...
%! PiecewiseLinearDistribution ([0, Inf], [0, 1])
%!error <PiecewiseLinearDistribution: X must be a real vector.> ...
%! PiecewiseLinearDistribution (["a", "c"], [0, 1])
%!error <PiecewiseLinearDistribution: X must be a real vector.> ...
%! PiecewiseLinearDistribution ([NaN, 1], [0, 1])
%!error <PiecewiseLinearDistribution: Fx must be a real vector.> ...
%! PiecewiseLinearDistribution ([0, 1], [0, i])
%!error <PiecewiseLinearDistribution: Fx must be a real vector.> ...
%! PiecewiseLinearDistribution ([0, 1], [0, Inf])
%!error <PiecewiseLinearDistribution: Fx must be a real vector.> ...
%! PiecewiseLinearDistribution ([0, 1], ["a", "c"])
%!error <PiecewiseLinearDistribution: Fx must be a real vector.> ...
%! PiecewiseLinearDistribution ([0, 1], [NaN, 1])
%!error <PiecewiseLinearDistribution: X and FX must be vectors of equal size.> ...
%! PiecewiseLinearDistribution ([0, 1], [0, 0.5, 1])
%!error <PiecewiseLinearDistribution: X and FX must be at least two-elements long.> ...
%! PiecewiseLinearDistribution ([0], [1])
%!error <PiecewiseLinearDistribution: FX must be bounded in the range> ...
%! PiecewiseLinearDistribution ([0, 0.5, 1], [0, 1, 1.5])

## 'cdf' method
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (PiecewiseLinearDistribution, 2, "uper")
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (PiecewiseLinearDistribution, 2, 3)

## 'plot' method
%!error <plot: optional arguments must be in NAME-VALUE pairs.> ...
%! plot (PiecewiseLinearDistribution, "Parent")
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (PiecewiseLinearDistribution, "PlotType", 12)
%!error <plot: invalid VALUE size for 'Parameter' argument.> ...
%! plot (PiecewiseLinearDistribution, "PlotType", {"pdf", "cdf"})
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (PiecewiseLinearDistribution, "PlotType", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (PiecewiseLinearDistribution, "Discrete", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (PiecewiseLinearDistribution, "Discrete", [1, 0])
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (PiecewiseLinearDistribution, "Discrete", {true})
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (PiecewiseLinearDistribution, "Parent", 12)
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (PiecewiseLinearDistribution, "Parent", "hax")
%!error <plot: invalid NAME for optional argument.> ...
%! plot (PiecewiseLinearDistribution, "invalidNAME", "pdf")
%!error <plot: 'probability' PlotType is not supported for 'PiecewiseLinearDistribution'.> ...
%! plot (PiecewiseLinearDistribution, "PlotType", "probability")

## 'truncate' method
%!error <truncate: missing input argument.> ...
%! truncate (PiecewiseLinearDistribution)
%!error <truncate: missing input argument.> ...
%! truncate (PiecewiseLinearDistribution, 2)
%!error <truncate: invalid lower upper limits.> ...
%! truncate (PiecewiseLinearDistribution, 4, 2)

## Catch errors when using array of probability objects with available methods
%!shared pd
%! pd = PiecewiseLinearDistribution ();
%! pd(2) = PiecewiseLinearDistribution ();
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
