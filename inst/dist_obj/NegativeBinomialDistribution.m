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

classdef NegativeBinomialDistribution
  ## -*- texinfo -*-
  ## @deftypefn {statistics} NegativeBinomialDistribution
  ##
  ## Negative binomial probability distribution object.
  ##
  ## A @code{NegativeBinomialDistribution} object consists of parameters, a
  ## model description, and sample data for a negative binomial probability
  ## distribution.
  ##
  ## The negative binomial distribution uses the following parameters.
  ##
  ## @multitable @columnfractions 0.25 0.48 0.27
  ## @headitem @var{Parameter} @tab @var{Description} @tab @var{Support}
  ##
  ## @item @qcode{R} @tab Number of successes @tab @math{R > 0}
  ## @item @qcode{P} @tab Probability of success @tab @math{0 < P <= 1}
  ## @end multitable
  ##
  ## There are several ways to create a @code{NegativeBinomialDistribution}
  ## object.
  ##
  ## @itemize
  ## @item Fit a distribution to data using the @code{fitdist} function.
  ## @item Create a distribution with specified parameter values using the
  ## @code{makedist} function.
  ## @item Use the constructor @qcode{NegativeBinomialDistribution (@var{R},
  ## @var{P})} to create a negative binomial distribution with specified
  ## parameter values.
  ## @item Use the static method @qcode{NegativeBinomialDistribution.fit
  ## (@var{x}, @var{censor}, @var{freq}, @var{options})} to a distribution to
  ## data @var{x}.
  ## @end itemize
  ##
  ## It is highly recommended to use @code{fitdist} and @code{makedist}
  ## functions to create probability distribution objects, instead of the
  ## constructor and the aforementioned static method.
  ##
  ## A @code{NegativeBinomialDistribution} object contains the following
  ## properties, which can be accessed using dot notation.
  ##
  ## @multitable @columnfractions 0.25 0.25 0.25 0.25
  ## @item @qcode{DistributionName} @tab @qcode{DistributionCode} @tab
  ## @qcode{NumParameters} @tab @qcode{ParameterNames}
  ## @item @qcode{ParameterDescription} @tab @qcode{ParameterValues} @tab
  ## @qcode{ParameterValues} @tab @qcode{ParameterCI}
  ## @item @qcode{ParameterIsFixed} @tab @qcode{Truncation} @tab
  ## @qcode{IsTruncated} @tab @qcode{InputData}
  ## @end multitable
  ##
  ## A @code{NegativeBinomialDistribution} object contains the following methods:
  ## @code{cdf}, @code{icdf}, @code{iqr}, @code{mean}, @code{median},
  ## @code{negloglik}, @code{paramci}, @code{pdf}, @code{plot}, @code{proflik},
  ## @code{random}, @code{std}, @code{truncate}, @code{var}.
  ##
  ## Further information about the negative binomial distribution can be found
  ## at @url{https://en.wikipedia.org/wiki/Negative_binomial_distribution}
  ##
  ## @seealso{fitdist, makedist, nbincdf, nbininv, nbinpdf, nbinrnd, nbinfit,
  ## nbinlike, nbinstat}
  ## @end deftypefn

  properties (Dependent = true)
    R
    P
  endproperties

  properties (GetAccess = public, Constant = true)
    CensoringAllowed = false;
    DistributionName = "NegativeBinomialDistribution";
    DistributionCode = "nbin";
    NumParameters = 2;
    ParameterNames = {"R", "P"};
    ParameterDescription = {"Number of successes", "Probability of success"};
  endproperties

  properties (GetAccess = public, Constant = true)
    ParameterRange = [realmin, realmin; Inf, 1];
    ParameterLogCI = [true, true];
  endproperties

  properties (GetAccess = public , SetAccess = protected)
    ParameterValues
    ParameterCI
    ParameterCovariance
    ParameterIsFixed
    Truncation
    IsTruncated
    InputData
  endproperties

  methods (Hidden)

    function this = NegativeBinomialDistribution (R, P)
      if (nargin == 0)
        R = 1;
        P = 0.5;
      endif
      checkparams (R, P);
      this.InputData = [];
      this.IsTruncated = false;
      this.ParameterValues = [R, P];
      this.ParameterIsFixed = [true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function display (this)
      fprintf ("%s =\n", inputname(1));
      __disp__ (this, "negative binomial distribution");
    endfunction

    function disp (this)
      __disp__ (this, "negative binomial distribution");
    endfunction

    function this = set.R (this, R)
      checkparams (R, this.P);
      this.InputData = [];
      this.ParameterValues(1) = R;
      this.ParameterIsFixed = [true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function R = get.R (this)
      R = this.ParameterValues(1);
    endfunction

    function this = set.P (this, P)
      checkparams (this.R, P);
      this.InputData = [];
      this.ParameterValues(2) = P;
      this.ParameterIsFixed = [true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function P = get.P (this)
      P = this.ParameterValues(2);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {NegativeBinomialDistribution} {@var{p} =} cdf (@var{pd}, @var{x})
    ## @deftypefnx {NegativeBinomialDistribution} {@var{p} =} cdf (@var{pd}, @var{x}, @qcode{"upper"})
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
      p = nbincdf (x, this.R, this.P);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        p(lb) = 0;
        p(ub) = 1;
        p(! (lb | ub)) -= nbincdf (lx, this.R, this.P);
        p(! (lb | ub)) /= diff (nbincdf ([lx, ux], this.R, this.P));
      endif
      ## Apply uflag
      if (utail)
        p = 1 - p;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {NegativeBinomialDistribution} {@var{p} =} icdf (@var{pd}, @var{p})
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
        lp = nbincdf (this.Truncation(1), this.R, this.P);
        up = nbincdf (this.Truncation(2), this.R, this.P);
        ## Adjust p values within range of p @ lower limit and p @ upper limit
        is_nan = p < 0 | p > 1;
        p(is_nan) = NaN;
        np = lp + (up - lp) .* p;
        x = nbininv (np, this.R, this.P);
        x(x < this.Truncation(1)) = this.Truncation(1);
        x(x > this.Truncation(2)) = this.Truncation(2);
      else
        x = nbininv (p, this.R, this.P);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {NegativeBinomialDistribution} {@var{r} =} iqr (@var{pd})
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
    ## @deftypefn  {NegativeBinomialDistribution} {@var{m} =} mean (@var{pd})
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
        m = nbinstat (this.R, this.P);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {NegativeBinomialDistribution} {@var{m} =} median (@var{pd})
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
        Fa_b = nbincdf ([lx, ux], this.R, this.P);
        m = nbininv (sum (Fa_b) / 2, this.R, this.P);
      else
        m = nbininv (0.5, this.R, this.P);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {NegativeBinomialDistribution} {@var{nlogL} =} negloglik (@var{pd})
    ##
    ## Compute the negative loglikelihood of a probability distribution.
    ##
    ## @code{@var{m} = negloglik (@var{pd})} computes the negative loglikelihood
    ## of the probability distribution object, @var{pd}.
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
      nlogL = - nbinlike ([this.R, this.P], this.InputData.data);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {NegativeBinomialDistribution} {@var{ci} =} paramci (@var{pd})
    ## @deftypefnx {NegativeBinomialDistribution} {@var{ci} =} paramci (@var{pd}, @var{Name}, @var{Value})
    ##
    ## Compute the confidence intervals for probability distribution parameters.
    ##
    ## @code{@var{ci} = paramci (@var{pd})} computes the lower and upper
    ## boundaries of the 95% confidence interval for each parameter of the
    ## probability distribution object, @var{pd}.
    ##
    ## @code{@var{ci} = paramci (@var{pd}, @var{Name}, @var{Value})} computes the
    ## confidence intervals with additional options specified specified by
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
    ## @deftypefn  {NegativeBinomialDistribution} {@var{y} =} pdf (@var{pd}, @var{x})
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
      y = nbinpdf (x, this.R, this.P);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        y(lb | ub) = 0;
        y(! (lb | ub)) /= diff (nbincdf ([lx, ux], this.R, this.P));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {NegativeBinomialDistribution} {} plot (@var{pd})
    ## @deftypefnx {NegativeBinomialDistribution} {} plot (@var{pd}, @var{Name}, @var{Value})
    ## @deftypefnx {NegativeBinomialDistribution} {@var{h} =} plot (@dots{})
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
    ## @deftypefn  {NegativeBinomialDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum})
    ## @deftypefnx {NegativeBinomialDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @qcode{"Display"}, @var{display})
    ## @deftypefnx {NegativeBinomialDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @var{setparam})
    ## @deftypefnx {NegativeBinomialDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @var{setparam}, @qcode{"Display"}, @var{display})
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
    ## For the negative binomial distribution, @qcode{@var{pnum} = 1} selects
    ## the parameter @qcode{R} and @qcode{@var{pnum} = 2} selects the parameter
    ## @var{P}.
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
    ## @deftypefn  {NegativeBinomialDistribution} {@var{y} =} random (@var{pd})
    ## @deftypefnx {NegativeBinomialDistribution} {@var{y} =} random (@var{pd}, @var{rows})
    ## @deftypefnx {NegativeBinomialDistribution} {@var{y} =} random (@var{pd}, @var{rows}, @var{cols}, @dots{})
    ## @deftypefnx {NegativeBinomialDistribution} {@var{y} =} random (@var{pd}, [@var{sz}])
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
        ratio = 1 / diff (poisscdf ([lx, ux], this.R, this.P));
        nsize = fix (2 * ratio * ps);       # times 2 to be on the safe side
        ## Generate the numbers and remove out-of-bound random samples
        r = nbinrnd (this.R, this.P, nsize, 1);
        r(r < lx | r > ux) = [];
        ## Randomly select the required size and reshape to requested dimensions
        idx = randperm (numel (r), ps);
        r = reshape (r(idx), sz);
      else
        r = nbinrnd (this.R, this.P, varargin{:});
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {NegativeBinomialDistribution} {@var{s} =} std (@var{pd})
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
    ## @deftypefn  {NegativeBinomialDistribution} {@var{t} =} truncate (@var{pd}, @var{lower}, @var{upper})
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
      ## Check boundaries and constrain within support: Natural numbers
      lower = round (lower);
      upper = round (upper);
      lower(lower < 0) = 0;
      if (lower >= upper)
        error ("truncate: invalid lower upper limits.");
      endif
      this.Truncation = [lower, upper];
      this.IsTruncated = true;
      this.InputData = [];
      this.ParameterIsFixed = [true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {NegativeBinomialDistribution} {@var{v} =} var (@var{pd})
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
        [~, v] = nbinstat (this.R, this.P);
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
      [phat, pci] = nbinfit (x, alpha, freq, options);
      [~, acov] = nbinlike (phat, x, freq);
      ## Create fitted distribution object
      pd = NegativeBinomialDistribution.makeFitted (phat, pci, acov, x, freq);
    endfunction

    function pd = makeFitted (phat, pci, acov, x, freq)
      R = phat(1);
      P = phat(2);
      pd = NegativeBinomialDistribution (R, P);
      pd.ParameterCI = pci;
      pd.ParameterIsFixed = [false, false];
      pd.ParameterCovariance = acov;
      pd.InputData = struct ("data", x, "cens", [], "freq", freq);
    endfunction

  endmethods

endclassdef

function checkparams (R, P)
  if (! (isscalar (R) && isnumeric (R) && isreal (R) && isfinite (R) && R > 0))
    error ("NegativeBinomialDistribution: R must be a positive scalar.")
  endif
  if (! (isscalar (P) && isnumeric (P) && isreal (P) && isfinite (P)
                      && P > 0 && P <= 1))
    error (strcat (["NegativeBinomialDistribution: P must be a real"], ...
                   [" scalar bounded in the range (0, 1]."]))
  endif
endfunction

## Test output
%!shared pd, t
%! pd = NegativeBinomialDistribution;
%! t = truncate (pd, 2, 4);
%!assert (cdf (pd, [0:5]), [0.5, 0.75, 0.875, 0.9375, 0.9688, 0.9844], 1e-4);
#%!assert (cdf (t, [0:5]), [0, 0, 0.5714, 0.8571, 1, 1], 1e-4);
%!assert (cdf (pd, [1.5, 2, 3, 4]), [0.75, 0.875, 0.9375, 0.9688], 1e-4);
#%!assert (cdf (t, [1.5, 2, 3, 4]), [0, 0.5714, 0.8571, 1], 1e-4);
%!assert (icdf (pd, [0:0.2:1]), [0, 0, 0, 1, 2, Inf], 1e-4);
#%!assert (icdf (t, [0:0.2:1]), [2, 2, 2, 3, 3, 4], 1e-4);
%!assert (icdf (pd, [-1, 0.4:0.2:1, NaN]), [NaN, 0, 1, 2, Inf, NaN], 1e-4);
#%!assert (icdf (t, [-1, 0.4:0.2:1, NaN]), [NaN, 2, 3, 3, 4, NaN], 1e-4);
%!assert (iqr (pd), 1);
#%!assert (iqr (t), 1);
%!assert (mean (pd), 1);
#%!assert (mean (t), 2.5714, 1e-4);
%!assert (median (pd), 0);
#%!assert (median (t), 2);
%!assert (pdf (pd, [0:5]), [0.5, 0.25, 0.125, 0.0625, 0.0312, 0.0156], 1e-4);
#%!assert (pdf (t, [0:5]), [0, 0, 0.5714, 0.2857, 0.1429, 0], 1e-4);
%!assert (pdf (pd, [-1, 1:4, NaN]), [0, 0.25, 0.125, 0.0625, 0.0312, NaN], 1e-4);
#%!assert (pdf (t, [-1, 1:4, NaN]), [0, 0, 0.5714, 0.2857, 0.1429, NaN], 1e-4);
%!assert (isequal (size (random (pd, 100, 50)), [100, 50]))
#%!assert (any (random (t, 1000, 1) < 2), false);
#%!assert (any (random (t, 1000, 1) > 4), false);
%!assert (std (pd), 1.4142, 1e-4);
#%!assert (std (t), 0.7284, 1e-4);
%!assert (var (pd), 2);
#%!assert (var (t), 0.5306, 1e-4);

## Test input validation
## 'NegativeBinomialDistribution' constructor
%!error <NegativeBinomialDistribution: R must be a positive scalar.> ...
%! NegativeBinomialDistribution(Inf, 1)
%!error <NegativeBinomialDistribution: R must be a positive scalar.> ...
%! NegativeBinomialDistribution(i, 1)
%!error <NegativeBinomialDistribution: R must be a positive scalar.> ...
%! NegativeBinomialDistribution("a", 1)
%!error <NegativeBinomialDistribution: R must be a positive scalar.> ...
%! NegativeBinomialDistribution([1, 2], 1)
%!error <NegativeBinomialDistribution: R must be a positive scalar.> ...
%! NegativeBinomialDistribution(NaN, 1)
%!error <NegativeBinomialDistribution: P must be a real scalar bounded in the range> ...
%! NegativeBinomialDistribution(1, 0)
%!error <NegativeBinomialDistribution: P must be a real scalar bounded in the range> ...
%! NegativeBinomialDistribution(1, -1)
%!error <NegativeBinomialDistribution: P must be a real scalar bounded in the range> ...
%! NegativeBinomialDistribution(1, Inf)
%!error <NegativeBinomialDistribution: P must be a real scalar bounded in the range> ...
%! NegativeBinomialDistribution(1, i)
%!error <NegativeBinomialDistribution: P must be a real scalar bounded in the range> ...
%! NegativeBinomialDistribution(1, "a")
%!error <NegativeBinomialDistribution: P must be a real scalar bounded in the range> ...
%! NegativeBinomialDistribution(1, [1, 2])
%!error <NegativeBinomialDistribution: P must be a real scalar bounded in the range> ...
%! NegativeBinomialDistribution(1, NaN)
%!error <NegativeBinomialDistribution: P must be a real scalar bounded in the range> ...
%! NegativeBinomialDistribution(1, 1.2)

## 'cdf' method
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (NegativeBinomialDistribution, 2, "uper")
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (NegativeBinomialDistribution, 2, 3)

## 'paramci' method
%!shared x
%! x = nbinrnd (1, 0.5, [1, 100]);
%!error <paramci: optional arguments must be in NAME-VALUE pairs.> ...
%! paramci (NegativeBinomialDistribution.fit (x), "alpha")
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (NegativeBinomialDistribution.fit (x), "alpha", 0)
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (NegativeBinomialDistribution.fit (x), "alpha", 1)
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (NegativeBinomialDistribution.fit (x), "alpha", [0.5 2])
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (NegativeBinomialDistribution.fit (x), "alpha", "")
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (NegativeBinomialDistribution.fit (x), "alpha", {0.05})
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (NegativeBinomialDistribution.fit (x), "parameter", "R", ...
%!          "alpha", {0.05})
%!error <paramci: invalid VALUE size for 'Parameter' argument.> ...
%! paramci (NegativeBinomialDistribution.fit (x), ...
%!          "parameter", {"R", "P", "param"})
%!error <paramci: invalid VALUE size for 'Parameter' argument.> ...
%! paramci (NegativeBinomialDistribution.fit (x), "alpha", 0.01, ...
%!          "parameter", {"R", "P", "param"})
%!error <paramci: unknown distribution parameter.> ...
%! paramci (NegativeBinomialDistribution.fit (x), "parameter", "param")
%!error <paramci: unknown distribution parameter.> ...
%! paramci (NegativeBinomialDistribution.fit (x), "alpha", 0.01, ...
%!          "parameter", "param")
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (NegativeBinomialDistribution.fit (x), "NAME", "value")
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (NegativeBinomialDistribution.fit (x), "alpha", 0.01, ...
%!          "NAME", "value")
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (NegativeBinomialDistribution.fit (x), "alpha", 0.01, ...
%!          "parameter", "R", "NAME", "value")

## 'plot' method
%!error <plot: optional arguments must be in NAME-VALUE pairs.> ...
%! plot (NegativeBinomialDistribution, "Parent")
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (NegativeBinomialDistribution, "PlotType", 12)
%!error <plot: invalid VALUE size for 'Parameter' argument.> ...
%! plot (NegativeBinomialDistribution, "PlotType", {"pdf", "cdf"})
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (NegativeBinomialDistribution, "PlotType", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (NegativeBinomialDistribution, "Discrete", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (NegativeBinomialDistribution, "Discrete", [1, 0])
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (NegativeBinomialDistribution, "Discrete", {true})
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (NegativeBinomialDistribution, "Parent", 12)
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (NegativeBinomialDistribution, "Parent", "hax")

## 'proflik' method
%!error <proflik: no fitted data available.> ...
%! proflik (NegativeBinomialDistribution, 2)
%!error <proflik: PNUM must be a scalar number indexing a non-fixed parameter.> ...
%! proflik (NegativeBinomialDistribution.fit (x), 3)
%!error <proflik: PNUM must be a scalar number indexing a non-fixed parameter.> ...
%! proflik (NegativeBinomialDistribution.fit (x), [1, 2])
%!error <proflik: PNUM must be a scalar number indexing a non-fixed parameter.> ...
%! proflik (NegativeBinomialDistribution.fit (x), {1})
%!error <proflik: SETPARAM must be a numeric vector.> ...
%! proflik (NegativeBinomialDistribution.fit (x), 1, ones (2))
%!error <proflik: missing VALUE for 'Display' argument.> ...
%! proflik (NegativeBinomialDistribution.fit (x), 1, "Display")
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (NegativeBinomialDistribution.fit (x), 1, "Display", 1)
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (NegativeBinomialDistribution.fit (x), 1, "Display", {1})
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (NegativeBinomialDistribution.fit (x), 1, "Display", {"on"})
%!error <proflik: invalid VALUE size for 'Display' argument.> ...
%! proflik (NegativeBinomialDistribution.fit (x), 1, "Display", ["on"; "on"])
%!error <proflik: invalid VALUE for 'Display' argument.> ...
%! proflik (NegativeBinomialDistribution.fit (x), 1, "Display", "onnn")
%!error <proflik: invalid NAME for optional arguments.> ...
%! proflik (NegativeBinomialDistribution.fit (x), 1, "NAME", "on")
%!error <proflik: invalid optional argument.> ...
%! proflik (NegativeBinomialDistribution.fit (x), 1, {"NAME"}, "on")
%!error <proflik: invalid optional argument.> ...
%! proflik (NegativeBinomialDistribution.fit (x), 1, {[1 2 3]}, "Display", "on")

## 'truncate' method
%!error <truncate: missing input argument.> ...
%! truncate (NegativeBinomialDistribution)
%!error <truncate: missing input argument.> ...
%! truncate (NegativeBinomialDistribution, 2)
%!error <truncate: invalid lower upper limits.> ...
%! truncate (NegativeBinomialDistribution, 4, 2)

## Catch errors when using array of probability objects with available methods
%!shared pd
%! pd = NegativeBinomialDistribution(1, 0.5);
%! pd(2) = NegativeBinomialDistribution(1, 0.6);
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
