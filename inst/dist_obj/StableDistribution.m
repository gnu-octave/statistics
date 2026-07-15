## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

classdef StableDistribution
  ## -*- texinfo -*-
  ## @deftp {statistics} StableDistribution
  ##
  ## Stable probability distribution object.
  ##
  ## A @code{StableDistribution} object consists of parameters, a model
  ## description, and sample data for a stable probability distribution.
  ##
  ## The stable distribution is a continuous probability distribution family
  ## closed under linear combinations, generalizing the normal, Cauchy, and Levy
  ## distributions.  It is parameterized, in the Nolan @qcode{S0}
  ## parameterization, by a tail index (first shape parameter) @var{alpha} in
  ## @math{(0, 2]}, a skewness (second shape parameter) @var{beta} in
  ## @math{[-1, 1]}, a scale parameter @var{gam} greater than zero, and a
  ## location parameter @var{delta}.
  ##
  ## There are several ways to create a @code{StableDistribution} object.
  ##
  ## @itemize
  ## @item Fit a distribution to data using the @code{fitdist} function.
  ## @item Create a distribution with fixed parameter values using the
  ## @code{makedist} function.
  ## @item Use the constructor @qcode{StableDistribution (@var{alpha},
  ## @var{beta}, @var{gam}, @var{delta})} to create a stable distribution with
  ## fixed parameter values @var{alpha}, @var{beta}, @var{gam}, and @var{delta}.
  ## @item Use the static method @qcode{StableDistribution.fit (@var{x},
  ## @var{alpha}, @var{freq}, @var{options})} to fit a distribution to the data
  ## in @var{x} using the same input arguments as the @code{stblfit} function.
  ## @end itemize
  ##
  ## It is highly recommended to use @code{fitdist} and @code{makedist}
  ## functions to create probability distribution objects, instead of the class
  ## constructor or the aforementioned static method.
  ##
  ## Fitting is by maximum likelihood.  Because the stable density has no closed
  ## form, it is evaluated by numerical inversion of the characteristic function,
  ## which makes fitting considerably slower than for the closed-form
  ## distributions.
  ##
  ## Further information about the stable distribution can be found at
  ## @url{https://en.wikipedia.org/wiki/Stable_distribution}
  ##
  ## @seealso{fitdist, makedist, stblpdf, stblcdf, stblinv, stblrnd, stblfit,
  ## stbllike}
  ## @end deftp

  properties (Dependent = true)
    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} alpha
    ##
    ## Tail index (first shape parameter)
    ##
    ## A scalar value in the range @math{(0, 2]} characterizing the tail
    ## behaviour of the stable distribution.  You can access the @qcode{alpha}
    ## property using dot name assignment.
    ##
    ## @end deftp
    alpha

    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} beta
    ##
    ## Skewness (second shape parameter)
    ##
    ## A scalar value in the range @math{[-1, 1]} characterizing the skewness of
    ## the stable distribution.  You can access the @qcode{beta} property using
    ## dot name assignment.
    ##
    ## @end deftp
    beta

    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} gam
    ##
    ## Scale parameter
    ##
    ## A positive scalar value characterizing the scale of the stable
    ## distribution.  You can access the @qcode{gam} property using dot name
    ## assignment.
    ##
    ## @end deftp
    gam

    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} delta
    ##
    ## Location parameter
    ##
    ## A scalar value characterizing the location of the stable distribution.
    ## You can access the @qcode{delta} property using dot name assignment.
    ##
    ## @end deftp
    delta
  endproperties

  properties (GetAccess = public, Constant = true)
    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} DistributionName
    ##
    ## Probability distribution name
    ##
    ## A character vector specifying the name of the probability distribution
    ## object.  This property is read-only.
    ##
    ## @end deftp
    DistributionName = "StableDistribution";

    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} NumParameters
    ##
    ## Number of parameters
    ##
    ## A scalar integer value specifying the number of parameters characterizing
    ## the probability distribution.  This property is read-only.
    ##
    ## @end deftp
    NumParameters = 4;

    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} ParameterNames
    ##
    ## Names of parameters
    ##
    ## A @math{4x1} cell array of character vectors with each element containing
    ## the name of a distribution parameter.  This property is read-only.
    ##
    ## @end deftp
    ParameterNames = {"alpha", "beta", "gam", "delta"};

    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} ParameterDescription
    ##
    ## Description of parameters
    ##
    ## A @math{4x1} cell array of character vectors with each element containing
    ## a short description of a distribution parameter.  This property is
    ## read-only.
    ##
    ## @end deftp
    ParameterDescription = {"First shape parameter", ...
                            "Second shape parameter", "Scale", "Location"};
  endproperties

  properties (GetAccess = public, Constant = true, Hidden)
    CensoringAllowed = false;
    DistributionCode = "stbl";
    ParameterRange = [realmin, -1, realmin, -Inf; 2, 1, Inf, Inf];
    ParameterLogCI = [false, false, false, false];
  endproperties

  properties (GetAccess = public, SetAccess = protected)
    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} ParameterValues
    ##
    ## Distribution parameter values
    ##
    ## A @math{4x1} numeric vector containing the values of the distribution
    ## parameters, matching the order in @qcode{ParameterNames}.  This property
    ## is read-only; use dot name assignment on the @qcode{alpha}, @qcode{beta},
    ## @qcode{gam}, and @qcode{delta} properties.
    ##
    ## @end deftp
    ParameterValues

    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} ParameterCovariance
    ##
    ## Covariance matrix of the parameter estimates
    ##
    ## A @math{4x4} numeric matrix containing the variance-covariance of the
    ## distribution parameters.  This property is read-only.
    ##
    ## @end deftp
    ParameterCovariance

    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} ParameterCI
    ##
    ## Confidence intervals for the parameter estimates
    ##
    ## A @math{2x4} numeric matrix containing the lower and upper boundaries of
    ## the confidence interval for each parameter, when the distribution has been
    ## fitted to data.  This property is read-only.
    ##
    ## @end deftp
    ParameterCI

    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} ParameterIsFixed
    ##
    ## Flags for fixed parameters
    ##
    ## A @math{4x1} logical vector specifying which parameters are held fixed
    ## rather than estimated.  This property is read-only.
    ##
    ## @end deftp
    ParameterIsFixed

    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} Truncation
    ##
    ## Truncation interval
    ##
    ## A two-element numeric vector with the truncation interval, if the
    ## distribution is truncated.  This property is read-only.
    ##
    ## @end deftp
    Truncation

    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} IsTruncated
    ##
    ## Flag for truncated distribution
    ##
    ## A logical scalar that is true when the distribution is truncated.  This
    ## property is read-only.
    ##
    ## @end deftp
    IsTruncated

    ## -*- texinfo -*-
    ## @deftp {StableDistribution} {property} InputData
    ##
    ## Data used for fitting the distribution
    ##
    ## A structure containing the data used to fit the distribution.  It is empty
    ## unless the distribution was fitted with @code{fitdist} or the static
    ## @code{fit} method.  This property is read-only.
    ##
    ## @end deftp
    InputData
  endproperties

  methods (Hidden)

    function this = StableDistribution (alpha, beta, gam, delta)
      if (nargin == 0)
        alpha = 2;
        beta = 0;
        gam = 1;
        delta = 0;
      endif
      checkparams (alpha, beta, gam, delta);
      this.InputData = [];
      this.IsTruncated = false;
      this.ParameterValues = [alpha, beta, gam, delta];
      this.ParameterIsFixed = [true, true, true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function display (this)
      __disp__ (this, "Stable distribution");
    endfunction

    function disp (this)
      __disp__ (this, "Stable distribution");
    endfunction

    function this = set.alpha (this, alpha)
      checkparams (alpha, this.beta, this.gam, this.delta);
      this.InputData = [];
      this.ParameterValues(1) = alpha;
      this.ParameterIsFixed = [true, true, true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function alpha = get.alpha (this)
      alpha = this.ParameterValues(1);
    endfunction

    function this = set.beta (this, beta)
      checkparams (this.alpha, beta, this.gam, this.delta);
      this.InputData = [];
      this.ParameterValues(2) = beta;
      this.ParameterIsFixed = [true, true, true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function beta = get.beta (this)
      beta = this.ParameterValues(2);
    endfunction

    function this = set.gam (this, gam)
      checkparams (this.alpha, this.beta, gam, this.delta);
      this.InputData = [];
      this.ParameterValues(3) = gam;
      this.ParameterIsFixed = [true, true, true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function gam = get.gam (this)
      gam = this.ParameterValues(3);
    endfunction

    function this = set.delta (this, delta)
      checkparams (this.alpha, this.beta, this.gam, delta);
      this.InputData = [];
      this.ParameterValues(4) = delta;
      this.ParameterIsFixed = [true, true, true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function delta = get.delta (this)
      delta = this.ParameterValues(4);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {StableDistribution} {@var{p} =} cdf (@var{pd}, @var{x})
    ## @deftypefnx {StableDistribution} {@var{p} =} cdf (@var{pd}, @var{x}, @qcode{"upper"})
    ##
    ## Compute the cumulative distribution function (CDF).
    ##
    ## @code{@var{p} = cdf (@var{pd}, @var{x})} computes the CDF of the
    ## probability distribution object, @var{pd}, evaluated at the values in
    ## @var{x}.  The optional @qcode{"upper"} flag computes the upper tail
    ## probability.
    ##
    ## @end deftypefn
    function p = cdf (this, x, uflag)
      if (! isscalar (this))
        error ("cdf: requires a scalar probability distribution.");
      endif
      if (nargin > 2 && strcmpi (uflag, "upper"))
        utail = true;
      elseif (nargin > 2 && ! strcmpi (uflag, "upper"))
        error ("cdf: invalid argument for upper tail.");
      else
        utail = false;
      endif
      p = stblcdf (x, this.alpha, this.beta, this.gam, this.delta);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        p(lb) = 0;
        p(ub) = 1;
        p(! (lb | ub)) -= stblcdf (lx, this.alpha, this.beta, this.gam, ...
                                   this.delta);
        p(! (lb | ub)) /= diff (stblcdf ([lx, ux], this.alpha, this.beta, ...
                                         this.gam, this.delta));
      endif
      if (utail)
        p = 1 - p;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {StableDistribution} {@var{x} =} icdf (@var{pd}, @var{p})
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
        lp = stblcdf (this.Truncation(1), this.alpha, this.beta, this.gam, ...
                      this.delta);
        up = stblcdf (this.Truncation(2), this.alpha, this.beta, this.gam, ...
                      this.delta);
        p(p < 0 | p > 1) = NaN;
        np = lp + (up - lp) .* p;
        x = stblinv (np, this.alpha, this.beta, this.gam, this.delta);
        x(x < this.Truncation(1)) = this.Truncation(1);
        x(x > this.Truncation(2)) = this.Truncation(2);
      else
        x = stblinv (p, this.alpha, this.beta, this.gam, this.delta);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {StableDistribution} {@var{r} =} iqr (@var{pd})
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
    ## @deftypefn  {StableDistribution} {@var{m} =} mean (@var{pd})
    ##
    ## Compute the mean of a probability distribution.
    ##
    ## @code{@var{m} = mean (@var{pd})} computes the mean of the probability
    ## distribution object, @var{pd}.  The mean is @qcode{NaN} for
    ## @code{@var{alpha} <= 1}, where it is undefined.
    ##
    ## @end deftypefn
    function m = mean (this)
      if (! isscalar (this))
        error ("mean: requires a scalar probability distribution.");
      endif
      if (this.IsTruncated)
        fm = @(x) x .* pdf (this, x);
        m = integral (fm, this.Truncation(1), this.Truncation(2));
      elseif (this.alpha == 2)
        m = this.delta;
      elseif (this.alpha > 1)
        m = this.delta - this.beta .* this.gam .* tan (pi .* this.alpha ./ 2);
      else
        m = NaN;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {StableDistribution} {@var{m} =} median (@var{pd})
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
        Fa_b = stblcdf ([lx, ux], this.alpha, this.beta, this.gam, this.delta);
        m = stblinv (sum (Fa_b) / 2, this.alpha, this.beta, this.gam, ...
                     this.delta);
      else
        m = stblinv (0.5, this.alpha, this.beta, this.gam, this.delta);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {StableDistribution} {@var{nlogL} =} negloglik (@var{pd})
    ##
    ## Compute the negative loglikelihood of a probability distribution.
    ##
    ## @code{@var{nlogL} = negloglik (@var{pd})} computes the negative
    ## loglikelihood of the probability distribution object, @var{pd}.  It
    ## returns an empty value when @var{pd} is not fitted to data.
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
      nlogL = stbllike ([this.alpha, this.beta, this.gam, this.delta], ...
                        this.InputData.data, this.InputData.freq);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {StableDistribution} {@var{ci} =} paramci (@var{pd})
    ## @deftypefnx {StableDistribution} {@var{ci} =} paramci (@var{pd}, @var{Name}, @var{Value})
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
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'Alpha'} @tab A scalar value in the range @math{(0,1)}
    ## specifying the significance level for the confidence interval.  The
    ## default value 0.05 corresponds to a 95% confidence interval.
    ##
    ## @item @qcode{'Parameter'} @tab A character vector or a cell array of
    ## character vectors specifying the parameter names for which to compute
    ## confidence intervals.  By default, @code{paramci} computes confidence
    ## intervals for all distribution parameters.
    ## @end multitable
    ##
    ## @code{paramci} is meaningful only when @var{pd} is fitted to data,
    ## otherwise the parameter values are returned in both rows.
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
    ## @deftypefn  {StableDistribution} {@var{y} =} pdf (@var{pd}, @var{x})
    ##
    ## Compute the probability density function (PDF).
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
      y = stblpdf (x, this.alpha, this.beta, this.gam, this.delta);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        ux = this.Truncation(2);
        y(x < lx | x > ux) = 0;
        y /= diff (stblcdf ([lx, ux], this.alpha, this.beta, this.gam, ...
                            this.delta));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {StableDistribution} {} plot (@var{pd})
    ## @deftypefnx {StableDistribution} {} plot (@var{pd}, @var{Name}, @var{Value})
    ## @deftypefnx {StableDistribution} {@var{h} =} plot (@dots{})
    ##
    ## Plot a probability distribution object.
    ##
    ## @code{plot (@var{pd})} plots the probability density function (PDF) of
    ## the probability distribution object @var{pd}.  Name-value pair arguments
    ## select the plotted function and its appearance, as documented in
    ## @code{__plot__}.
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
    ## @deftypefn  {StableDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum})
    ## @deftypefnx {StableDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @qcode{'Display'}, @var{display})
    ## @deftypefnx {StableDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @var{setparam})
    ## @deftypefnx {StableDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @var{setparam}, @qcode{'Display'}, @var{display})
    ##
    ## Profile likelihood function for a probability distribution object.
    ##
    ## @code{[@var{nlogL}, @var{param}] = proflik (@var{pd}, @var{pnum})} returns
    ## a vector @var{nlogL} of negative loglikelihood values and a vector
    ## @var{param} of corresponding parameter values for the parameter in the
    ## position indicated by @var{pnum}.  By default, @code{proflik} uses the
    ## lower and upper bounds of the 95% confidence interval and computes 100
    ## equispaced values for the selected parameter.  @var{pd} must be fitted to
    ## data.
    ##
    ## @code{[@var{nlogL}, @var{param}] = proflik (@var{pd}, @var{pnum},
    ## @qcode{'Display'}, @qcode{'on'})} also plots the profile likelihood
    ## against the default range of the selected parameter.
    ##
    ## @code{[@var{nlogL}, @var{param}] = proflik (@var{pd}, @var{pnum},
    ## @var{setparam})} defines a user-defined range of the selected parameter.
    ##
    ## @code{[@var{nlogL}, @var{param}] = proflik (@var{pd}, @var{pnum},
    ## @var{setparam}, @qcode{'Display'}, @qcode{'on'})} also plots the profile
    ## likelihood against the user-defined range of the selected parameter.
    ##
    ## For the stable distribution, @qcode{@var{pnum} = 1} selects the tail index
    ## @qcode{alpha}, @qcode{@var{pnum} = 2} selects the skewness @qcode{beta},
    ## @qcode{@var{pnum} = 3} selects the scale @qcode{gam}, and
    ## @qcode{@var{pnum} = 4} selects the location @qcode{delta}.
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
    ## @deftypefn  {StableDistribution} {@var{r} =} random (@var{pd})
    ## @deftypefnx {StableDistribution} {@var{r} =} random (@var{pd}, @var{rows})
    ## @deftypefnx {StableDistribution} {@var{r} =} random (@var{pd}, @var{rows}, @var{cols}, @dots{})
    ## @deftypefnx {StableDistribution} {@var{r} =} random (@var{pd}, [@var{sz}])
    ##
    ## Generate random arrays from the probability distribution object.
    ##
    ## @code{@var{r} = random (@var{pd})} returns a random number from the
    ## distribution object @var{pd}, following the size conventions of
    ## @code{stblrnd}.
    ##
    ## @end deftypefn
    function r = random (this, varargin)
      if (! isscalar (this))
        error ("random: requires a scalar probability distribution.");
      endif
      if (this.IsTruncated)
        sz = [varargin{:}];
        ps = prod (sz);
        lx = this.Truncation(1);
        ux = this.Truncation(2);
        ratio = 1 / diff (stblcdf ([lx, ux], this.alpha, this.beta, ...
                                   this.gam, this.delta));
        nsize = fix (2 * ratio * ps);
        r = stblrnd (this.alpha, this.beta, this.gam, this.delta, nsize, 1);
        r(r < lx | r > ux) = [];
        while (numel (r) < ps)
          r = [r; stblrnd(this.alpha, this.beta, this.gam, this.delta, ...
                          nsize, 1)];
          r(r < lx | r > ux) = [];
        endwhile
        idx = randperm (numel (r), ps);
        r = reshape (r(idx), sz);
      else
        r = stblrnd (this.alpha, this.beta, this.gam, this.delta, varargin{:});
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {StableDistribution} {@var{s} =} std (@var{pd})
    ##
    ## Compute the standard deviation of a probability distribution.
    ##
    ## @code{@var{s} = std (@var{pd})} computes the standard deviation of the
    ## probability distribution object, @var{pd}.  It is @qcode{NaN} for
    ## @code{@var{alpha} < 2}, where the variance is infinite.
    ##
    ## @end deftypefn
    function s = std (this)
      if (! isscalar (this))
        error ("std: requires a scalar probability distribution.");
      endif
      s = sqrt (var (this));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {StableDistribution} {@var{t} =} truncate (@var{pd}, @var{lower}, @var{upper})
    ##
    ## Truncate a probability distribution.
    ##
    ## @code{@var{t} = truncate (@var{pd}, @var{lower}, @var{upper})} returns
    ## the probability distribution @var{pd} truncated to the interval with
    ## lower limit @var{lower} and upper limit @var{upper}.
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
      this.InputData = [];
      this.ParameterIsFixed = [true, true, true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {StableDistribution} {@var{v} =} var (@var{pd})
    ##
    ## Compute the variance of a probability distribution.
    ##
    ## @code{@var{v} = var (@var{pd})} computes the variance of the probability
    ## distribution object, @var{pd}.  It is @qcode{NaN} for @code{@var{alpha} <
    ## 2}, where the variance is infinite.
    ##
    ## @end deftypefn
    function v = var (this)
      if (! isscalar (this))
        error ("var: requires a scalar probability distribution.");
      endif
      if (this.IsTruncated)
        fm = @(x) x .* pdf (this, x);
        m = integral (fm, this.Truncation(1), this.Truncation(2));
        fv = @(x) ((x - m) .^ 2) .* pdf (this, x);
        v = integral (fv, this.Truncation(1), this.Truncation(2));
      elseif (this.alpha == 2)
        v = 2 .* this.gam .^ 2;
      else
        v = NaN;
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
        freq = [];
      else
        freq = varargin{2};
      endif
      if (nargin < 4)
        options.Display = "off";
        options.MaxFunEvals = 400;
        options.MaxIter = 200;
        options.TolX = 1e-6;
      else
        options = varargin{3};
      endif
      ## Fit data
      [phat, pci] = stblfit (x, alpha, freq, options);
      [~, acov] = stbllike (phat, x, freq);
      ## Create fitted distribution object
      pd = StableDistribution.makeFitted (phat, pci, acov, x, freq);
    endfunction

    function pd = makeFitted (phat, pci, acov, x, freq)
      alpha = phat(1);
      beta = phat(2);
      gam = phat(3);
      delta = phat(4);
      pd = StableDistribution (alpha, beta, gam, delta);
      pd.ParameterCI = pci;
      pd.ParameterIsFixed = [false, false, false, false];
      pd.ParameterCovariance = acov;
      pd.InputData = struct ("data", x, "cens", [], "freq", freq);
    endfunction

  endmethods

endclassdef

function checkparams (alpha, beta, gam, delta)
  if (! (isscalar (alpha) && isnumeric (alpha) && isreal (alpha) ...
         && alpha > 0 && alpha <= 2))
    error ("StableDistribution: ALPHA must be a real scalar in (0, 2].");
  endif
  if (! (isscalar (beta) && isnumeric (beta) && isreal (beta) ...
         && beta >= -1 && beta <= 1))
    error ("StableDistribution: BETA must be a real scalar in [-1, 1].");
  endif
  if (! (isscalar (gam) && isnumeric (gam) && isreal (gam) ...
         && isfinite (gam) && gam > 0))
    error ("StableDistribution: GAM must be a positive real scalar.");
  endif
  if (! (isscalar (delta) && isnumeric (delta) && isreal (delta) ...
         && isfinite (delta)))
    error ("StableDistribution: DELTA must be a real scalar.");
  endif
endfunction

%!demo
%! ## Create a stable distribution and plot its pdf
%! pd = makedist ("Stable", "alpha", 1.5, "beta", 0.5, "gam", 1, "delta", 0);
%! plot (pd);
%! title ("Stable distribution, alpha = 1.5, beta = 0.5");

## Test output against MATLAB (created via makedist)
%!shared pd
%! pd = makedist ("Stable", "alpha", 1.5, "beta", 0.5, "gam", 1, "delta", 0);
%!test
%! assert (pd.alpha, 1.5);
%! assert (pd.beta, 0.5);
%! assert (pd.gam, 1);
%! assert (pd.delta, 0);
%! assert (pd.DistributionName, "StableDistribution");
%! assert (pd.NumParameters, 4);
%!test
%! x = -5:5;
%! exp_p = [0.00961772128347771, 0.0143422747723476, 0.0257902242195547, ...
%!          0.0657154294128386, 0.201576145758624, 0.462186560100778, ...
%!          0.712063555515659, 0.855535196378772, 0.921201224725992, ...
%!          0.951409668616683, 0.966845678836178];
%! assert (cdf (pd, x), exp_p, 1e-8);
%!test
%! assert (icdf (pd, [0.1, 0.5, 0.9]), ...
%!         [-1.63127009138493, 0.133853042315326, 2.58231785139714], 1e-6);
%!test  # mean, variance, median of a skewed stable (alpha = 1.5)
%! assert (mean (pd), 0.5, 1e-12);
%! assert (isnan (var (pd)));
%! assert (median (pd), 0.133853042315326, 1e-6);

## alpha = 2 is the normal distribution with variance 2*gam^2
%!test
%! pn = makedist ("Stable", "alpha", 2, "beta", 0, "gam", 1, "delta", 0);
%! assert (mean (pn), 0);
%! assert (var (pn), 2, 1e-12);
%! assert (std (pn), sqrt (2), 1e-12);
%! assert (pdf (pn, 0:2), normpdf (0:2, 0, sqrt (2)), 1e-12);

## alpha <= 1 has no finite mean
%!test
%! pc = makedist ("Stable", "alpha", 0.8, "beta", 0.5, "gam", 1, "delta", 0);
%! assert (isnan (mean (pc)));
%! assert (isnan (var (pc)));
%! assert (median (pc), 0.250487323305453, 1e-5);

## Scale and location, and a truncated distribution
%!test
%! ps = makedist ("Stable", "alpha", 1.5, "beta", 0.5, "gam", 2, "delta", 3);
%! assert (mean (ps), 4, 1e-12);
%! assert (median (ps), 3.26770608463065, 1e-6);
%!test  # truncation renormalizes and bounds the support
%! pt = truncate (pd, -1, 3);
%! assert (pt.IsTruncated);
%! assert (cdf (pt, [-2, 3]), [0, 1]);
%! assert (isfinite (mean (pt)));
%! assert (all (random (pt, 100, 1) >= -1 & random (pt, 100, 1) <= 3));

## Test input validation
%!error <StableDistribution: ALPHA must be a real scalar in \(0, 2\].> ...
%! StableDistribution (2.5, 0, 1, 0)
%!error <StableDistribution: BETA must be a real scalar in \[-1, 1\].> ...
%! StableDistribution (1.5, 2, 1, 0)
%!error <StableDistribution: GAM must be a positive real scalar.> ...
%! StableDistribution (1.5, 0, 0, 0)
%!error <StableDistribution: DELTA must be a real scalar.> ...
%! StableDistribution (1.5, 0, 1, Inf)
%!error <cdf: requires a scalar probability distribution.> ...
%! cdf ([pd, pd], 1)
