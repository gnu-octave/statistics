## Copyright (C) 2024-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

classdef ExponentialDistribution
  ## -*- texinfo -*-
  ## @deftp {statistics} ExponentialDistribution
  ##
  ## Exponential probability distribution object.
  ##
  ## A @code{ExponentialDistribution} object consists of parameters, a model
  ## description, and sample data for a exponential probability distribution.
  ##
  ## The exponential distribution is a continuous probability distribution with
  ## mean parameter @var{mu} that models the time between events in a Poisson
  ## process.
  ##
  ## There are several ways to create a @code{ExponentialDistribution} object.
  ##
  ## @itemize
  ## @item Fit a distribution to data using the @code{fitdist} function.
  ## @item Create a distribution with fixed parameter values using the
  ## @code{makedist} function.
  ## @item Use the constructor @qcode{ExponentialDistribution (@var{mu})}
  ## to create a exponential distribution with fixed parameter value @var{mu}.
  ## @item Use the static method @qcode{ExponentialDistribution.fit (@var{x},
  ## @var{alpha}, @var{censor}, @var{freq}, @var{options})} to fit a
  ## distribution to the data in @var{x} using the same input arguments as the
  ## @code{expfit} function.
  ## @end itemize
  ##
  ## It is highly recommended to use @code{fitdist} and @code{makedist}
  ## functions to create probability distribution objects, instead of the class
  ## constructor or the aforementioned static method.
  ##
  ## Further information about the exponential distribution can be found at
  ## @url{https://en.wikipedia.org/wiki/Exponential_distribution}
  ##
  ## @seealso{fitdist, makedist, expcdf, expinv, exppdf, exprnd, expfit,
  ## explike, expstat}
  ## @end deftp

  properties (Dependent = true)
    ## -*- texinfo -*-
    ## @deftp {ExponentialDistribution} {property} mu
    ##
    ## Mean parameter
    ##
    ## A positive scalar value characterizing the mean of the
    ## exponential distribution. You can access the @qcode{mu}
    ## property using dot name assignment.
    ##
    ## @end deftp
    mu
  endproperties

  properties (GetAccess = public, Constant = true)
    ## -*- texinfo -*-
    ## @deftp {ExponentialDistribution} {property} DistributionName
    ##
    ## Probability distribution name
    ##
    ## A character vector specifying the name of the probability distribution
    ## object. This property is read-only.
    ##
    ## @end deftp
    DistributionName = "ExponentialDistribution";

    ## -*- texinfo -*-
    ## @deftp {ExponentialDistribution} {property} NumParameters
    ##
    ## Number of parameters
    ##
    ## A scalar integer value specifying the number of parameters characterizing
    ## the probability distribution. This property is read-only.
    ##
    ## @end deftp
    NumParameters = 1;

    ## -*- texinfo -*-
    ## @deftp {ExponentialDistribution} {property} ParameterNames
    ##
    ## Names of parameters
    ##
    ## A @math{1x1} cell array of character vectors with each element containing
    ## the name of a distribution parameter. This property is read-only.
    ##
    ## @end deftp
    ParameterNames = {"mu"};

    ## -*- texinfo -*-
    ## @deftp {ExponentialDistribution} {property} ParameterDescription
    ##
    ## Description of parameters
    ##
    ## A @math{1x1} cell array of character vectors with each element containing
    ## a short description of a distribution parameter. This property is
    ## read-only.
    ##
    ## @end deftp
    ParameterDescription = {"Mean"};
  endproperties

  properties (GetAccess = public, Constant = true, Hidden)
    CensoringAllowed = true;
    DistributionCode = "exp";
    ParameterRange = [realmin; Inf];
    ParameterLogCI = true;
  endproperties

  properties (GetAccess = public , SetAccess = protected)
    ## -*- texinfo -*-
    ## @deftp {ExponentialDistribution} {property} ParameterValues
    ##
    ## Distribution parameter values
    ##
    ## A @math{1x1} numeric vector containing the value of the distribution
    ## parameter. This property is read-only. You can change the distribution
    ## parameter by assigning a new value to the @qcode{mu}
    ## property.
    ##
    ## @end deftp
    ParameterValues

    ## -*- texinfo -*-
    ## @deftp {ExponentialDistribution} {property} ParameterCovariance
    ##
    ## Covariance matrix of the parameter estimates
    ##
    ## A scalar numeric value containing the variance-covariance of the
    ## parameter estimate. The covariance matrix is only meaningful
    ## when the distribution was fitted to data. If the distribution object was
    ## created with fixed parameters, or a parameter of a fitted distribution is
    ## modified, then the variance-covariance is zero. This
    ## property is read-only.
    ##
    ## @end deftp
    ParameterCovariance

    ## -*- texinfo -*-
    ## @deftp {ExponentialDistribution} {property} ParameterIsFixed
    ##
    ## Flag for fixed parameters
    ##
    ## A @math{1x1} logical vector specifying whether the parameter is fixed or
    ## estimated. @qcode{true} value corresponds to fixed parameter,
    ## @qcode{false} value corresponds to parameter estimate. This property is
    ## read-only.
    ##
    ## @end deftp
    ParameterIsFixed

    ## -*- texinfo -*-
    ## @deftp {ExponentialDistribution} {property} Truncation
    ##
    ## Truncation interval
    ##
    ## A @math{1x2} numeric vector specifying the truncation interval for the
    ## probability distribution. First element contains the lower boundary,
    ## second element contains the upper boundary. This property is read-only.
    ## You can only truncate a probability distribution with the
    ## @qcode{truncate} method.
    ##
    ## @end deftp
    Truncation

    ## -*- texinfo -*-
    ## @deftp {ExponentialDistribution} {property} IsTruncated
    ##
    ## Flag for truncated probability distribution
    ##
    ## A logical scalar value specifying whether a probability distribution is
    ## truncated or not. This property is read-only.
    ##
    ## @end deftp
    IsTruncated

    ## -*- texinfo -*-
    ## @deftp {ExponentialDistribution} {property} InputData
    ##
    ## Data used for fitting a probability distribution
    ##
    ## A scalar structure containing the following fields:
    ## @itemize
    ## @item @qcode{data}: a numeric vector containing the data used for
    ## distribution fitting.
    ## @item @qcode{cens}: a numeric vector of logical values indicating
    ## censoring information corresponding to the elements of the data used for
    ## distribution fitting. If no censoring vector was used for distribution
    ## fitting, then this field defaults to an empty array.
    ## @item @qcode{freq}: a numeric vector of non-negative integer values
    ## containing the frequency information corresponding to the elements of the
    ## data used for distribution fitting. If no frequency vector was used for
    ## distribution fitting, then this field defaults to an empty array.
    ## @end itemize
    ##
    ## @end deftp
    InputData
  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    ParameterCI
  endproperties

  methods (Hidden)

    function this = ExponentialDistribution (mu)
      if (nargin == 0)
        mu = 1;
      endif
      checkparams (mu);
      this.InputData = [];
      this.IsTruncated = false;
      this.ParameterValues = mu;
      this.ParameterIsFixed = true;
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function display (this)
      fprintf ("%s =\n", inputname(1));
      __disp__ (this, "exponential distribution");
    endfunction

    function disp (this)
      __disp__ (this, "exponential distribution");
    endfunction

    function this = set.mu (this, mu)
      checkparams (mu);
      this.InputData = [];
      this.ParameterValues(1) = mu;
      this.ParameterIsFixed = true;
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function mu = get.mu (this)
      mu = this.ParameterValues(1);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {ExponentialDistribution} {@var{p} =} cdf (@var{pd}, @var{x})
    ## @deftypefnx {ExponentialDistribution} {@var{p} =} cdf (@var{pd}, @var{x}, @qcode{"upper"})
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
      p = expcdf (x, this.mu);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        p(lb) = 0;
        p(ub) = 1;
        p(! (lb | ub)) -= expcdf (lx, this.mu);
        p(! (lb | ub)) /= diff (expcdf ([lx, ux], this.mu));
      endif
      ## Apply uflag
      if (utail)
        p = 1 - p;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ExponentialDistribution} {@var{x} =} icdf (@var{pd}, @var{p})
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
        lp = expcdf (this.Truncation(1), this.mu);
        up = expcdf (this.Truncation(2), this.mu);
        ## Adjust p values within range of p @ lower limit and p @ upper limit
        is_nan = p < 0 | p > 1;
        p(is_nan) = NaN;
        np = lp + (up - lp) .* p;
        x = expinv (np, this.mu);
        x(x < this.Truncation(1)) = this.Truncation(1);
        x(x > this.Truncation(2)) = this.Truncation(2);
      else
        x = expinv (p, this.mu);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ExponentialDistribution} {@var{r} =} iqr (@var{pd})
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
    ## @deftypefn  {ExponentialDistribution} {@var{m} =} mean (@var{pd})
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
        m = expstat (this.mu);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ExponentialDistribution} {@var{m} =} median (@var{pd})
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
        Fa_b = expcdf ([lx, ux], this.mu);
        m = expinv (sum (Fa_b) / 2, this.mu);
      else
        m = this.mu .* log (2);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ExponentialDistribution} {@var{nlogL} =} negloglik (@var{pd})
    ##
    ## Compute the negative loglikelihood of a probability distribution.
    ##
    ## @code{@var{nlogL} = negloglik (@var{pd})} computes the negative
    ## loglikelihood of the probability distribution object, @var{pd}.
    ##
    ## @end deftypefn
    function nlogL = negloglik (this)
      if (! isscalar (this))
        error ("negloglik: requires a scalar probability distribution.");
      endif
      if (isempty (this.InputData))
        nlogL = [];
        return
      endif
      nlogL = - explike (this.mu, this.InputData.data, ...
                         this.InputData.cens, this.InputData.freq);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ExponentialDistribution} {@var{ci} =} paramci (@var{pd})
    ## @deftypefnx {ExponentialDistribution} {@var{ci} =} paramci (@var{pd}, @var{Name}, @var{Value})
    ##
    ## Compute the confidence intervals for probability distribution parameters.
    ##
    ## @code{@var{ci} = paramci (@var{pd})} computes the lower and upper
    ## boundaries of the 95% confidence interval for each parameter of the
    ## probability distribution object, @var{pd}.
    ##
    ## @code{@var{ci} = paramci (@var{pd}, @var{Name}, @var{Value})} computes the
    ## confidence intervals with additional options specified by
    ## @qcode{Name-Value} pair arguments listed below.
    ##
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"Alpha"} @tab @tab A scalar value in the range @math{(0,1)}
    ## specifying the significance level for the confidence interval.  The
    ## default value 0.05 corresponds to a 95% confidence interval.
    ##
    ## @item @qcode{"Parameter"} @tab @tab A character vector or a cell array of
    ## character vectors specifying the parameter names for which to compute
    ## confidence intervals.  By default, @code{paramci} computes confidence
    ## intervals for all distribution parameters.
    ## @end multitable
    ##
    ## @code{paramci} is meaningful only when @var{pd} is fitted to data,
    ## otherwise an empty array, @qcode{[]}, is returned.
    ##
    ## @end deftypefn
    function ci = paramci (this, varargin)
      if (! isscalar (this))
        error ("paramci: requires a scalar probability distribution.");
      endif
      if (isempty (this.InputData))
        ci = [this.ParameterValues; this.ParameterValues];
      else
        ci = __paramci__ (this, varargin{:});
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ExponentialDistribution} {@var{y} =} pdf (@var{pd}, @var{x})
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
      y = exppdf (x, this.mu);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        y(lb | ub) = 0;
        y(! (lb | ub)) /= diff (expcdf ([lx, ux], this.mu));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ExponentialDistribution} {} plot (@var{pd})
    ## @deftypefnx {ExponentialDistribution} {} plot (@var{pd}, @var{Name}, @var{Value})
    ## @deftypefnx {ExponentialDistribution} {@var{h} =} plot (@dots{})
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
    ## @deftypefn  {ExponentialDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum})
    ## @deftypefnx {ExponentialDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @qcode{"Display"}, @var{display})
    ## @deftypefnx {ExponentialDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @var{setparam})
    ## @deftypefnx {ExponentialDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @var{setparam}, @qcode{"Display"}, @var{display})
    ##
    ## Profile likelihood function for a probability distribution object.
    ##
    ## @code{[@var{nlogL}, @var{param}] = proflik (@var{pd}, @var{pnum})}
    ## returns a vector @var{nlogL} of negative loglikelihood values and a
    ## vector @var{param} of corresponding parameter values for the parameter in
    ## the position indicated by @var{pnum}.  By default, @code{proflik} uses
    ## the lower and upper bounds of the 95% confidence interval and computes
    ## 100 equispaced values for the selected parameter.  @var{pd} must be
    ## fitted to data.
    ##
    ## @code{[@var{nlogL}, @var{param}] = proflik (@var{pd}, @var{pnum},
    ## @qcode{"Display"}, @qcode{"on"})} also plots the profile likelihood
    ## against the default range of the selected parameter.
    ##
    ## @code{[@var{nlogL}, @var{param}] = proflik (@var{pd}, @var{pnum},
    ## @var{setparam})} defines a user-defined range of the selected parameter.
    ##
    ## @code{[@var{nlogL}, @var{param}] = proflik (@var{pd}, @var{pnum},
    ## @var{setparam}, @qcode{"Display"}, @qcode{"on"})} also plots the profile
    ## likelihood against the user-defined range of the selected parameter.
    ##
    ## For the exponential distribution, @qcode{@var{pnum} = 1} selects the
    ## parameter @qcode{mu}.
    ##
    ## When opted to display the profile likelihood plot, @code{proflik} also
    ## plots the baseline loglikelihood computed at the lower bound of the 95%
    ## confidence interval and estimated maximum likelihood.  The latter might
    ## not be observable if it is outside of the used-defined range of parameter
    ## values.
    ##
    ## @end deftypefn
    function [varargout] = proflik (this, pnum, varargin)
      if (! isscalar (this))
        error ("proflik: requires a scalar probability distribution.");
      endif
      if (isempty (this.InputData))
        error ("proflik: no fitted data available.");
      endif
      [varargout{1:nargout}] = __proflik__ (this, pnum, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ExponentialDistribution} {@var{r} =} random (@var{pd})
    ## @deftypefnx {ExponentialDistribution} {@var{r} =} random (@var{pd}, @var{rows})
    ## @deftypefnx {ExponentialDistribution} {@var{r} =} random (@var{pd}, @var{rows}, @var{cols}, @dots{})
    ## @deftypefnx {ExponentialDistribution} {@var{r} =} random (@var{pd}, [@var{sz}])
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
        ratio = 1 / diff (expcdf ([lx, ux], this.mu));
        nsize = fix (2 * ratio * ps);       # times 2 to be on the safe side
        ## Generate the numbers and remove out-of-bound random samples
        r = exprnd (this.mu, nsize, 1);
        r(r < lx | r > ux) = [];
        ## Randomly select the required size and reshape to requested dimensions
        idx = randperm (numel (r), ps);
        r = reshape (r(idx), sz);
      else
        r = exprnd (this.mu, varargin{:});
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ExponentialDistribution} {@var{s} =} std (@var{pd})
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
    ## @deftypefn  {ExponentialDistribution} {@var{t} =} truncate (@var{pd}, @var{lower}, @var{upper})
    ##
    ## Truncate a probability distribution.
    ##
    ## @code{@var{t} = truncate (@var{pd}, @var{lower}, @var{upper})} returns a
    ## probability distribution @var{t}, which is the probability distribution
    ## @var{pd} truncated to the specified interval with lower limit, @var{lower},
    ## and upper limit, @var{upper}.  If @var{pd} is fitted to data with
    ## @code{fitdist}, the returned probability distribution @var{t} is not
    ## fitted, does not contain any data or estimated values, and it is as it
    ## has been created with the @var{makedist} function, but it includes the
    ## truncation interval.
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
      ## Check boundaries and constrain within support [0, Inf)
      lower(lower < 0) = 0;
      this.Truncation = [lower, upper];
      this.IsTruncated = true;
      this.InputData = [];
      this.ParameterIsFixed = true;
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ExponentialDistribution} {@var{v} =} var (@var{pd})
    ##
    ## Compute the variance of a probability distribution.
    ##
    ## @code{@var{v} = var (@var{pd})} computes the variance of the
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
        [~, v] = expstat (this.mu);
      endif
    endfunction

  endmethods

  methods (Static, Hidden)

    function pd = fit (x, varargin)
      ## Check input arguments
      if (nargin < 2)
        alpha = 0.05;
      else
        alpha = varargin{1};
      endif
      if (nargin < 3)
        censor = [];
      else
        censor = varargin{2};
      endif
      if (nargin < 4)
        freq = [];
      else
        freq = varargin{3};
      endif
      ## Fit data
      [phat, pci] = expfit (x, alpha, censor, freq);
      [~, acov] = explike (phat, x, censor, freq);
      ## Create fitted distribution object
      pd = ExponentialDistribution.makeFitted ...
           (phat, pci, acov, x, censor, freq);
    endfunction

    function pd = makeFitted (phat, pci, acov, x, censor, freq)
      mu = phat(1);
      pd = ExponentialDistribution (mu);
      pd.ParameterCI = pci;
      pd.ParameterIsFixed = false;
      pd.ParameterCovariance = acov;
      pd.InputData = struct ("data", x, "cens", censor, "freq", freq);
    endfunction

  endmethods

endclassdef

function checkparams (mu)
  if (! (isscalar (mu) && isnumeric (mu) && isreal (mu)
                          && isfinite (mu) && mu > 0))
    error ("ExponentialDistribution: MU must be a positive real scalar.")
  endif
endfunction

## Test output
%!shared pd, t
%! pd = ExponentialDistribution (1);
%! t = truncate (pd, 2, 4);
%!assert (cdf (pd, [0:5]), [0, 0.6321, 0.8647, 0.9502, 0.9817, 0.9933], 1e-4);
%!assert (cdf (t, [0:5]), [0, 0, 0, 0.7311, 1, 1], 1e-4);
%!assert (cdf (pd, [1.5, 2, 3, 4]), [0.7769, 0.8647, 0.9502, 0.9817], 1e-4);
%!assert (cdf (t, [1.5, 2, 3, 4]), [0, 0, 0.7311, 1], 1e-4);
%!assert (icdf (pd, [0:0.2:1]), [0, 0.2231, 0.5108, 0.9163, 1.6094, Inf], 1e-4);
%!assert (icdf (t, [0:0.2:1]), [2, 2.1899, 2.4244, 2.7315, 3.1768, 4], 1e-4);
%!assert (icdf (pd, [-1, 0.4:0.2:1, NaN]), [NaN, 0.5108, 0.9163, 1.6094, Inf, NaN], 1e-4);
%!assert (icdf (t, [-1, 0.4:0.2:1, NaN]), [NaN, 2.4244, 2.7315, 3.1768, 4, NaN], 1e-4);
%!assert (iqr (pd), 1.0986, 1e-4);
%!assert (iqr (t), 0.8020, 1e-4);
%!assert (mean (pd), 1);
%!assert (mean (t), 2.6870, 1e-4);
%!assert (median (pd), 0.6931, 1e-4);
%!assert (median (t), 2.5662, 1e-4);
%!assert (pdf (pd, [0:5]), [1, 0.3679, 0.1353, 0.0498, 0.0183, 0.0067], 1e-4);
%!assert (pdf (t, [0:5]), [0, 0, 1.1565, 0.4255, 0.1565, 0], 1e-4);
%!assert (pdf (pd, [-1, 1:4, NaN]), [0, 0.3679, 0.1353, 0.0498, 0.0183, NaN], 1e-4);
%!assert (pdf (t, [-1, 1:4, NaN]), [0, 0, 1.1565, 0.4255, 0.1565, NaN], 1e-4);
%!assert (isequal (size (random (pd, 100, 50)), [100, 50]))
%!assert (any (random (t, 1000, 1) < 2), false);
%!assert (any (random (t, 1000, 1) > 4), false);
%!assert (std (pd), 1);
%!assert (std (t), 0.5253, 1e-4);
%!assert (var (pd), 1);
%!assert (var (t), 0.2759, 1e-4);

## Test input validation
## 'ExponentialDistribution' constructor
%!error <ExponentialDistribution: MU must be a positive real scalar.> ...
%! ExponentialDistribution(0)
%!error <ExponentialDistribution: MU must be a positive real scalar.> ...
%! ExponentialDistribution(-1)
%!error <ExponentialDistribution: MU must be a positive real scalar.> ...
%! ExponentialDistribution(Inf)
%!error <ExponentialDistribution: MU must be a positive real scalar.> ...
%! ExponentialDistribution(i)
%!error <ExponentialDistribution: MU must be a positive real scalar.> ...
%! ExponentialDistribution("a")
%!error <ExponentialDistribution: MU must be a positive real scalar.> ...
%! ExponentialDistribution([1, 2])
%!error <ExponentialDistribution: MU must be a positive real scalar.> ...
%! ExponentialDistribution(NaN)

## 'cdf' method
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (ExponentialDistribution, 2, "uper")
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (ExponentialDistribution, 2, 3)

## 'paramci' method
%!shared x
%! x = exprnd (1, [100, 1]);
%!error <paramci: optional arguments must be in NAME-VALUE pairs.> ...
%! paramci (ExponentialDistribution.fit (x), "alpha")
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (ExponentialDistribution.fit (x), "alpha", 0)
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (ExponentialDistribution.fit (x), "alpha", 1)
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (ExponentialDistribution.fit (x), "alpha", [0.5 2])
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (ExponentialDistribution.fit (x), "alpha", "")
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (ExponentialDistribution.fit (x), "alpha", {0.05})
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (ExponentialDistribution.fit (x), "parameter", "mu", ...
%!          "alpha", {0.05})
%!error <paramci: invalid VALUE size for 'Parameter' argument.> ...
%! paramci (ExponentialDistribution.fit (x), "parameter", {"mu", "param"})
%!error <paramci: invalid VALUE size for 'Parameter' argument.> ...
%! paramci (ExponentialDistribution.fit (x), "alpha", 0.01, ...
%!          "parameter", {"mu", "param"})
%!error <paramci: unknown distribution parameter.> ...
%! paramci (ExponentialDistribution.fit (x), "parameter", "param")
%!error <paramci: unknown distribution parameter.> ...
%! paramci (ExponentialDistribution.fit (x), "alpha", 0.01, "parameter", "parm")
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (ExponentialDistribution.fit (x), "NAME", "value")
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (ExponentialDistribution.fit (x), "alpha", 0.01, "NAME", "value")
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (ExponentialDistribution.fit (x), "alpha", 0.01, ...
%!          "parameter", "mu", "NAME", "value")

## 'plot' method
%!error <plot: optional arguments must be in NAME-VALUE pairs.> ...
%! plot (ExponentialDistribution, "Parent")
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (ExponentialDistribution, "PlotType", 12)
%!error <plot: invalid VALUE size for 'Parameter' argument.> ...
%! plot (ExponentialDistribution, "PlotType", {"pdf", "cdf"})
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (ExponentialDistribution, "PlotType", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (ExponentialDistribution, "Discrete", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (ExponentialDistribution, "Discrete", [1, 0])
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (ExponentialDistribution, "Discrete", {true})
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (ExponentialDistribution, "Parent", 12)
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (ExponentialDistribution, "Parent", "hax")
%!error <plot: invalid NAME for optional argument.> ...
%! plot (ExponentialDistribution, "invalidNAME", "pdf")
%!error <plot: no fitted DATA to plot a probability plot.> ...
%! plot (ExponentialDistribution, "PlotType", "probability")

## 'proflik' method
%!error <proflik: no fitted data available.> ...
%! proflik (ExponentialDistribution, 2)
%!error <proflik: PNUM must be a scalar number indexing a non-fixed parameter.> ...
%! proflik (ExponentialDistribution.fit (x), 3)
%!error <proflik: PNUM must be a scalar number indexing a non-fixed parameter.> ...
%! proflik (ExponentialDistribution.fit (x), [1, 2])
%!error <proflik: PNUM must be a scalar number indexing a non-fixed parameter.> ...
%! proflik (ExponentialDistribution.fit (x), {1})
%!error <proflik: SETPARAM must be a numeric vector.> ...
%! proflik (ExponentialDistribution.fit (x), 1, ones (2))
%!error <proflik: missing VALUE for 'Display' argument.> ...
%! proflik (ExponentialDistribution.fit (x), 1, "Display")
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (ExponentialDistribution.fit (x), 1, "Display", 1)
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (ExponentialDistribution.fit (x), 1, "Display", {1})
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (ExponentialDistribution.fit (x), 1, "Display", {"on"})
%!error <proflik: invalid VALUE size for 'Display' argument.> ...
%! proflik (ExponentialDistribution.fit (x), 1, "Display", ["on"; "on"])
%!error <proflik: invalid VALUE for 'Display' argument.> ...
%! proflik (ExponentialDistribution.fit (x), 1, "Display", "onnn")
%!error <proflik: invalid NAME for optional arguments.> ...
%! proflik (ExponentialDistribution.fit (x), 1, "NAME", "on")
%!error <proflik: invalid optional argument.> ...
%! proflik (ExponentialDistribution.fit (x), 1, {"NAME"}, "on")
%!error <proflik: invalid optional argument.> ...
%! proflik (ExponentialDistribution.fit (x), 1, {[1 2 3 4]}, "Display", "on")

## 'truncate' method
%!error <truncate: missing input argument.> ...
%! truncate (ExponentialDistribution)
%!error <truncate: missing input argument.> ...
%! truncate (ExponentialDistribution, 2)
%!error <truncate: invalid lower upper limits.> ...
%! truncate (ExponentialDistribution, 4, 2)

## Catch errors when using array of probability objects with available methods
%!shared pd
%! pd = ExponentialDistribution(1);
%! pd(2) = ExponentialDistribution(3);
%!error <cdf: requires a scalar probability distribution.> cdf (pd, 1)
%!error <icdf: requires a scalar probability distribution.> icdf (pd, 0.5)
%!error <iqr: requires a scalar probability distribution.> iqr (pd)
%!error <mean: requires a scalar probability distribution.> mean (pd)
%!error <median: requires a scalar probability distribution.> median (pd)
%!error <negloglik: requires a scalar probability distribution.> negloglik (pd)
%!error <paramci: requires a scalar probability distribution.> paramci (pd)
%!error <pdf: requires a scalar probability distribution.> pdf (pd, 1)
%!error <plot: requires a scalar probability distribution.> plot (pd)
%!error <proflik: requires a scalar probability distribution.> proflik (pd, 2)
%!error <random: requires a scalar probability distribution.> random (pd)
%!error <std: requires a scalar probability distribution.> std (pd)
%!error <truncate: requires a scalar probability distribution.> ...
%! truncate (pd, 2, 4)
%!error <var: requires a scalar probability distribution.> var (pd)
