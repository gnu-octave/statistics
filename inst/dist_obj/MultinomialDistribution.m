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

classdef MultinomialDistribution
  ## -*- texinfo -*-
  ## @deftp {statistics} MultinomialDistribution
  ##
  ## Multinomial probability distribution object.
  ##
  ## A @code{MultinomialDistribution} object consists of parameters, a model
  ## description, and sample data for a multinomial probability distribution.
  ##
  ## The multinomial distribution is a discrete probability distribution that
  ## models the outcomes of n independent trials of a k-category system, where
  ## each trial has a probability of falling into each category.  It is
  ## defined by the vector of probabilities for each outcome.
  ##
  ## There are several ways to create a @code{MultinomialDistribution} object.
  ##
  ## @itemize
  ## @item Create a distribution with specified parameter values using the
  ## @code{makedist} function.
  ## @item Use the constructor @qcode{MultinomialDistribution (@var{Probabilities})}
  ## to create a multinomial distribution with specified parameter values.
  ## @end itemize
  ##
  ## It is highly recommended to use the @code{makedist} function to create
  ## probability distribution objects, instead of the constructor.
  ##
  ## Further information about the multinomial distribution can be found at
  ## @url{https://en.wikipedia.org/wiki/Multinomial_distribution}
  ##
  ## @seealso{makedist, mnpdf, mnrnd}
  ## @end deftp

  properties (Dependent = true)
    ## -*- texinfo -*-
    ## @deftp {MultinomialDistribution} {property} Probabilities
    ##
    ## Outcome probabilities
    ##
    ## A row vector of probabilities for each outcome.  You can access the
    ## @qcode{Probabilities} property using dot name assignment.
    ##
    ## @end deftp
    Probabilities
  endproperties

  properties (GetAccess = public, Constant = true)
    ## -*- texinfo -*-
    ## @deftp {MultinomialDistribution} {property} DistributionName
    ##
    ## Probability distribution name
    ##
    ## A character vector specifying the name of the probability distribution
    ## object.  This property is read-only.
    ##
    ## @end deftp
    DistributionName = "MultinomialDistribution";

    ## -*- texinfo -*-
    ## @deftp {MultinomialDistribution} {property} NumParameters
    ##
    ## Number of parameters
    ##
    ## A scalar integer value specifying the number of parameters characterizing
    ## the probability distribution.  This property is read-only.
    ##
    ## @end deftp
    NumParameters = 1;

    ## -*- texinfo -*-
    ## @deftp {MultinomialDistribution} {property} ParameterNames
    ##
    ## Names of parameters
    ##
    ## A @math{1x1} cell array of character vectors with each element containing
    ## the name of a distribution parameter.  This property is read-only.
    ##
    ## @end deftp
    ParameterNames = {"Probabilities"};

    ## -*- texinfo -*-
    ## @deftp {MultinomialDistribution} {property} ParameterDescription
    ##
    ## Description of parameters
    ##
    ## A @math{1x1} cell array of character vectors with each element containing
    ## a short description of a distribution parameter.  This property is
    ## read-only.
    ##
    ## @end deftp
    ParameterDescription = {"Outcome probabilities"};
  endproperties

  properties (GetAccess = public, Constant = true, Hidden)
    CensoringAllowed = false;
    DistributionCode = "mn";
  endproperties

  properties (GetAccess = public , SetAccess = protected)
    ## -*- texinfo -*-
    ## @deftp {MultinomialDistribution} {property} ParameterValues
    ##
    ## Distribution parameter values
    ##
    ## A numeric vector containing the values of the distribution
    ## parameters.  This property is read-only. You can change the distribution
    ## parameters by assigning new values to the @qcode{Probabilities}
    ## property.
    ##
    ## @end deftp
    ParameterValues

    ## -*- texinfo -*-
    ## @deftp {MultinomialDistribution} {property} Truncation
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
    ## @deftp {MultinomialDistribution} {property} IsTruncated
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

    function this = MultinomialDistribution (Probabilities)
      if (nargin == 0)
        Probabilities = [0.5, 0.5];
      endif
      checkparams (Probabilities);
      this.IsTruncated = false;
      this.ParameterValues = Probabilities;
    endfunction

    function display (this)
      fprintf ("%s =\n", inputname(1));
      __disp__ (this, "multinomial distribution");
    endfunction

    function disp (this)
      __disp__ (this, "multinomial distribution");
    endfunction

    function this = set.Probabilities (this, Probabilities)
      checkparams (Probabilities);
      this.ParameterValues = Probabilities(:)';
    endfunction

    function Probabilities = get.Probabilities (this)
      Probabilities = this.ParameterValues(:)';
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialDistribution} {@var{p} =} cdf (@var{pd}, @var{x})
    ## @deftypefnx {MultinomialDistribution} {@var{p} =} cdf (@var{pd}, @var{x}, @qcode{"upper"})
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
      ## Check input data
      if (! isreal (x))
        error ("cdf: X must be real.");
      endif
      probs = this.Probabilities;
      ## Check for truncation and normalize truncated probabilities vector
      if (this.IsTruncated)
        lx = this.Truncation(1);
        ux = this.Truncation(2);
        probs = this.Probabilities([lx:ux]);
        probs = probs .* (1 / sum (probs));
        x = x - lx + 1;
      endif
      ## Do the computations
      is_nan = isnan (x);
      sz = size (x);
      xf = floor (x);
      pk = length (probs);
      ## Create cumulative probability vector
      pc = cumsum (probs);
      pc(end) = 1;  # Force last element to 1
      xf(xf > pk) = pk;
      xf(xf < 1) = 1;
      xf(is_nan) = 1;
      p = pc(xf);
      p(x < 1) = 0;
      p(is_nan) = NaN;
      p = reshape (p, sz);
      ## Apply uflag
      if (utail)
        p = 1 - p;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialDistribution} {@var{p} =} icdf (@var{pd}, @var{p})
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
      probs = this.Probabilities;
      ## Do the computations
      sz = size(p);
      p = p(:);
      x = zeros (numel (p), 1);
      pc = cumsum(this.Probabilities);
      pc = [0 pc(1:(end-1))];
      is_one = p == 0;
      is_nan = isnan (p) | p > 1 | p < 0;
      for i = 1:length (pc)
        x(p > pc(i)) = i;
      endfor
      x(is_one) = 1;
      x(is_nan) = NaN;
      ## Check for truncation and clip edges
      if (this.IsTruncated)
        lx = this.Truncation(1);
        ux = this.Truncation(2);
        lp = pc(lx);
        up = pc(ux);
        lb = p >= 0 & p <= lp;
        ub = p <= 1 & p >= up;
        x(lb) = lx;
        x(ub) = ux;
      endif
      x = reshape (x, sz);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialDistribution} {@var{r} =} iqr (@var{pd})
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
    ## @deftypefn  {MultinomialDistribution} {@var{m} =} mean (@var{pd})
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
      probs = this.Probabilities;
      if (this.IsTruncated)
        x = 1:numel (probs);
        w = x >= this.Truncation(1) & x <= this.Truncation(2);
        x = x(w);
        y = pdf (this, x);
        m = sum (y .* x);
      else
        m = sum (this.Probabilities .* (1:numel (probs)));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialDistribution} {@var{m} =} median (@var{pd})
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
    ## @deftypefn  {MultinomialDistribution} {@var{y} =} pdf (@var{pd}, @var{x})
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
      probs = this.Probabilities;
      size_x = size (x);
      is_nan = isnan (x);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        ux = this.Truncation(2);
        is_out = x < lx | x > ux | (x-floor (x)) > 0 | is_nan;
        tprobs = probs([lx:ux]);
        tprobs = tprobs .* (1 / sum (tprobs));
        probs([lx:ux]) = tprobs;
        copy_x = x;
        copy_x(is_out) = 1;
        y = probs(copy_x);
        y(is_out) = 0;
        y(is_nan) = NaN;
        y = reshape (y, size_x);
        return
      else
        is_out = x < 1 | x > length (probs) | (x-floor (x)) > 0 | is_nan;
        copy_x = x;
        copy_x(is_out) = 1;
        y = probs(copy_x);
        y(is_out) = 0;
        y(is_nan) = NaN;
        y = reshape (y, size_x);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialDistribution} {} plot (@var{pd})
    ## @deftypefnx {MultinomialDistribution} {} plot (@var{pd}, @var{Name}, @var{Value})
    ## @deftypefnx {MultinomialDistribution} {@var{h} =} plot (@dots{})
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
      h = __plot__ (this, true, varargin{:});
      if (nargout > 0)
        varargout{1} = h;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialDistribution} {@var{y} =} random (@var{pd})
    ## @deftypefnx {MultinomialDistribution} {@var{y} =} random (@var{pd}, @var{rows})
    ## @deftypefnx {MultinomialDistribution} {@var{y} =} random (@var{pd}, @var{rows}, @var{cols}, @dots{})
    ## @deftypefnx {MultinomialDistribution} {@var{y} =} random (@var{pd}, [@var{sz}])
    ##
    ## Generate random arrays from the probability distribution object.
    ##
    ## @code{@var{r} = random (@var{pd})} returns a random number from the
    ## distribution object @var{pd}.
    ##
    ## When called with a single size argument, @code{mnrnd} returns a square
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
      ps = prod (sz);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        ux = this.Truncation(2);
        tprobs = pdf (this, [1:numel(this.Probabilities)]);
        tprobs = tprobs([lx:ux]);
        cp = cumsum (tprobs);
      else
        cp = cumsum (this.Probabilities);
      endif
      bins = min ([0, cp], 1);
      bins(end) = 1;
      [~, r] = histc (rand (ps, 1), bins);
      r = reshape (r, sz);
      if (this.IsTruncated)
        r = r + lx - 1;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialDistribution} {@var{s} =} std (@var{pd})
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
    ## @deftypefn  {MultinomialDistribution} {@var{t} =} truncate (@var{pd}, @var{lower}, @var{upper})
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
        error ("truncate: is_nan input argument.");
      endif
      ## Constrain within the length of Probabilities vector
      lower = round (lower);
      upper = round (upper);
      k = numel (this.Probabilities);
      lower(lower < 1) = 1;
      upper(upper > k) = k;
      if (lower >= upper)
        error ("truncate: invalid lower upper limits.");
      endif
      this.Truncation = [lower, upper];
      this.IsTruncated = true;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialDistribution} {@var{v} =} var (@var{pd})
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
      probs = this.Probabilities;
      if (this.IsTruncated)
        x = 1:numel (probs);
        w = x >= this.Truncation(1) & x <= this.Truncation(2);
        x = x(w);
        y = pdf (this, x);
        m = sum (y .* x);
        v = sum (y .* (x - m) .^ 2);
      else
        v = sum (probs .* (1:numel (probs)) .^ 2) - mean (this) ^ 2;
      endif
    endfunction

  endmethods

endclassdef

function checkparams (Probabilities)
  if (! (isvector (Probabilities) && isnumeric (Probabilities) &&
         isreal (Probabilities) && isfinite (Probabilities) &&
         abs (sum (Probabilities) - 1) < eps * 100))
    error (["MultinomialDistribution: PROBABILITIES must be a vector", ...
            " of positive real scalars that sum up to 1."])
  endif
endfunction

%!demo
%! ## Generate a data set of 5000 random samples from a Multinomial distribution
%! ## with parameters Probabilities = [0.1, 0.2, 0.3, 0.2, 0.1, 0.1].  Create
%! ## the distribution and plot the PDF superimposed on a histogram of the data.
%!
%! probs = [0.1, 0.2, 0.3, 0.2, 0.1, 0.1];
%! pd = makedist ("Multinomial", "Probabilities", probs);
%! rand ("seed", 2);
%! data = random (pd, 5000, 1);
%! hist (data, length (probs));
%! hold on
%! x = 1:length (probs);
%! y = pdf (pd, x) * 5000;
%! stem (x, y, "r", "LineWidth", 2);
%! hold off
%! msg = "Multinomial distribution with Probabilities = [%s]";
%! probs_str = num2str (probs, "%0.1f ");
%! title (sprintf (msg, probs_str))

## Test output
%!shared pd, t
%! pd = MultinomialDistribution ([0.1, 0.2, 0.3, 0.2, 0.1, 0.1]);
%! t = truncate (pd, 2, 4);
%!assert (cdf (pd, [2, 3, 4]), [0.3, 0.6, 0.8], eps);
%!assert (cdf (t, [2, 3, 4]), [0.2857, 0.7143, 1], 1e-4);
%!assert (cdf (pd, [1.5, 2, 3, 4]), [0.1, 0.3, 0.6, 0.8], eps);
%!assert (cdf (pd, [1.5, 2-eps, 3, 4]), [0.1, 0.1, 0.6, 0.8], eps);
%!assert (cdf (t, [1.5, 2, 3, 4]), [0, 0.2857, 0.7143, 1], 1e-4);
%!assert (cdf (t, [1.5, 2-eps, 3, 4]), [0, 0, 0.7143, 1], 1e-4);
%!assert (cdf (pd, [1, 2.5, 4, 6]), [0.1, 0.3, 0.8, 1], eps);
%!assert (icdf (pd, [0, 0.2857, 0.7143, 1]), [1, 2, 4, 6]);
%!assert (icdf (t, [0, 0.2857, 0.7143, 1]), [2, 2, 4, 4]);
%!assert (icdf (t, [0, 0.35, 0.7143, 1]), [2, 3, 4, 4]);
%!assert (icdf (t, [0, 0.35, 0.7143, 1, NaN]), [2, 3, 4, 4, NaN]);
%!assert (icdf (t, [-0.5, 0, 0.35, 0.7143, 1, NaN]), [NaN, 2, 3, 4, 4, NaN]);
%!assert (icdf (pd, [-0.5, 0, 0.35, 0.7143, 1, NaN]), [NaN, 1, 3, 4, 6, NaN]);
%!assert (iqr (pd), 2);
%!assert (iqr (t), 2);
%!assert (mean (pd), 3.3, 1e-14);
%!assert (mean (t), 3, eps);
%!assert (median (pd), 3);
%!assert (median (t), 3);
%!assert (pdf (pd, [-5, 1, 2.5, 4, 6, NaN, 9]), [0, 0.1, 0, 0.2, 0.1, NaN, 0]);
%!assert (pdf (pd, [-5, 1, 2, 3, 4, 6, NaN, 9]), ...
%! [0, 0.1, 0.2, 0.3, 0.2, 0.1, NaN, 0]);
%!assert (pdf (t, [-5, 1, 2, 3, 4, 6, NaN, 0]), ...
%! [0, 0, 0.2857, 0.4286, 0.2857, 0, NaN, 0], 1e-4);
%!assert (pdf (t, [-5, 1, 2, 4, 6, NaN, 0]), ...
%! [0, 0, 0.2857, 0.2857, 0, NaN, 0], 1e-4);
%!assert (unique (random (pd, 1000, 5)), [1, 2, 3, 4, 5, 6]');
%!assert (unique (random (t, 1000, 5)), [2, 3, 4]');
%!assert (std (pd), 1.4177, 1e-4);
%!assert (std (t), 0.7559, 1e-4);
%!assert (var (pd), 2.0100, 1e-4);
%!assert (var (t), 0.5714, 1e-4);

## Test input validation
## 'MultinomialDistribution' constructor
%!error <MultinomialDistribution: PROBABILITIES must be a vector of positive real scalars that sum up to 1.> ...
%! MultinomialDistribution(0)
%!error <MultinomialDistribution: PROBABILITIES must be a vector of positive real scalars that sum up to 1.> ...
%! MultinomialDistribution(-1)
%!error <MultinomialDistribution: PROBABILITIES must be a vector of positive real scalars that sum up to 1.> ...
%! MultinomialDistribution(Inf)
%!error <MultinomialDistribution: PROBABILITIES must be a vector of positive real scalars that sum up to 1.> ...
%! MultinomialDistribution(i)
%!error <MultinomialDistribution: PROBABILITIES must be a vector of positive real scalars that sum up to 1.> ...
%! MultinomialDistribution("a")
%!error <MultinomialDistribution: PROBABILITIES must be a vector of positive real scalars that sum up to 1.> ...
%! MultinomialDistribution([1, 2])
%!error <MultinomialDistribution: PROBABILITIES must be a vector of positive real scalars that sum up to 1.> ...
%! MultinomialDistribution(NaN)

## 'cdf' method
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (MultinomialDistribution, 2, "uper")
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (MultinomialDistribution, 2, 3)
%!error <cdf: X must be real.> ...
%! cdf (MultinomialDistribution, i)

## 'plot' method
%!error <plot: optional arguments must be in NAME-VALUE pairs.> ...
%! plot (MultinomialDistribution, "Parent")
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (MultinomialDistribution, "PlotType", 12)
%!error <plot: invalid VALUE size for 'Parameter' argument.> ...
%! plot (MultinomialDistribution, "PlotType", {"pdf", "cdf"})
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (MultinomialDistribution, "PlotType", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (MultinomialDistribution, "Discrete", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (MultinomialDistribution, "Discrete", [1, 0])
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (MultinomialDistribution, "Discrete", {true})
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (MultinomialDistribution, "Parent", 12)
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (MultinomialDistribution, "Parent", "hax")
%!error <plot: invalid NAME for optional argument.> ...
%! plot (MultinomialDistribution, "invalidNAME", "pdf")
%!error <plot: 'probability' PlotType is not supported for 'MultinomialDistribution'.> ...
%! plot (MultinomialDistribution, "PlotType", "probability")

## 'truncate' method
%!error <truncate: is_nan input argument.> ...
%! truncate (MultinomialDistribution)
%!error <truncate: is_nan input argument.> ...
%! truncate (MultinomialDistribution, 2)
%!error <truncate: invalid lower upper limits.> ...
%! truncate (MultinomialDistribution, 4, 2)

## Catch errors when using array of probability objects with available methods
%!shared pd
%! pd = MultinomialDistribution([0.1, 0.2, 0.3, 0.4]);
%! pd(2) = MultinomialDistribution([0.1, 0.2, 0.3, 0.4]);
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
