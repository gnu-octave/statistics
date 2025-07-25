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

classdef BinomialDistribution
  ## -*- texinfo -*-
  ## @deftypefn {statistics} BinomialDistribution
  ##
  ## Binomial probability distribution object.
  ##
  ## A @code{BinomialDistribution} object consists of parameters, a model
  ## description, and sample data for a binomial probability distribution.
  ##
  ## The binomial distribution is a discrete probability distribution that
  ## models the number of successes in a sequence of @math{n} independent
  ## trials, eachwith a probability of success @math{p}.
  ##
  ## There are several ways to create a @code{BinomialDistribution} object.
  ##
  ## @itemize
  ## @item Fit a distribution to data using the @code{fitdist} function.
  ## @item Create a distribution with specified parameter values using the
  ## @code{makedist} function.
  ## @item Use the constructor @qcode{BinomialDistribution (@var{N}, @var{p})}
  ## to create a binomial distribution with specified parameter values.
  ## @item Use the static method @qcode{BinomialDistribution.fit (@var{x},
  ## @var{censor}, @var{freq}, @var{options})} to fit a distribution to data
  ## @var{x}.
  ## @end itemize
  ##
  ## It is highly recommended to use @code{fitdist} and @code{makedist}
  ## functions to create probability distribution objects, instead of the
  ## constructor and the aforementioned static method.
  ##
  ## Further information about the binomial distribution can be found at
  ## @url{https://en.wikipedia.org/wiki/Binomial_distribution}
  ##
  ## @seealso{fitdist, makedist, binocdf, binoinv, binopdf, binornd, binofit,
  ## binolike, binostat}
  ## @end deftypefn

  properties (Dependent = true)
    ## -*- texinfo -*-
    ## @deftp {BinomialDistribution} {property} N
    ##
    ## Number of trials
    ##
    ## A positive integer value characterizing the number of trials in the
    ## binomial distribution. You can access the @qcode{N} property using dot
    ## name assignment.
    ##
    ## @end deftp
    N

    ## -*- texinfo -*-
    ## @deftp {BinomialDistribution} {property} p
    ##
    ## Probability of success
    ##
    ## A scalar value in the range @math{[0, 1]} characterizing the probability
    ## of success in each trial of the binomial distribution. You can access
    ## the @qcode{p} property using dot name assignment.
    ##
    ## @end deftp
    p
  endproperties

  properties (GetAccess = public, Constant = true)
    ## -*- texinfo -*-
    ## @deftp {BinomialDistribution} {property} DistributionName
    ##
    ## Probability distribution name
    ##
    ## A character vector specifying the name of the probability distribution
    ## object. This property is read-only.
    ##
    ## @end deftp
    DistributionName = "BinomialDistribution";

    ## -*- texinfo -*-
    ## @deftp {BinomialDistribution} {property} NumParameters
    ##
    ## Number of parameters
    ##
    ## A scalar integer value specifying the number of parameters characterizing
    ## the probability distribution. This property is read-only.
    ##
    ## @end deftp
    NumParameters = 2;

    ## -*- texinfo -*-
    ## @deftp {BinomialDistribution} {property} ParameterNames
    ##
    ## Names of parameters
    ##
    ## A @math{2x1} cell array of character vectors with each element containing
    ## the name of a distribution parameter. This property is read-only.
    ##
    ## @end deftp
    ParameterNames = {"N", "p"};

    ## -*- texinfo -*-
    ## @deftp {BinomialDistribution} {property} ParameterDescription
    ##
    ## Description of parameters
    ##
    ## A @math{2x1} cell array of character vectors with each element containing
    ## a short description of a distribution parameter. This property is
    ## read-only.
    ##
    ## @end deftp
    ParameterDescription = {"Number of trials", "Probability of success"};
  endproperties

  properties (GetAccess = public, Constant = true, Hidden)
    CensoringAllowed = false;
    DistributionCode = "bino";
    ParameterRange = [realmin, realmin; Inf, 1];
    ParameterLogCI = [true, true];
  endproperties

  properties (GetAccess = public, SetAccess = protected)
    ## -*- texinfo -*-
    ## @deftp {BinomialDistribution} {property} ParameterValues
    ##
    ## Distribution parameter values
    ##
    ## A @math{2x1} numeric vector containing the values of the distribution
    ## parameters. This property is read-only. You can change the distribution
    ## parameters by assigning new values to the @qcode{N} and @qcode{p}
    ## properties.
    ##
    ## @end deftp
    ParameterValues

    ## -*- texinfo -*-
    ## @deftp {BinomialDistribution} {property} ParameterCovariance
    ##
    ## Covariance matrix of the parameter estimates
    ##
    ## A @math{2x2} numeric matrix containing the variance-covariance of the
    ## parameter estimates. Diagonal elements contain the variance of each
    ## estimated parameter and non-diagonal elements contain the covariance
    ## between the parameter estimates. The covariance matrix is only
    ## meaningful when the distribution was fitted to data. If the distribution
    ## object was created with fixed parameters, or a parameter of a fitted
    ## distribution is modified, then all elements of the variance-covariance
    ## are zero. This property is read-only.
    ##
    ## @end deftp
    ParameterCovariance

    ## -*- texinfo -*-
    ## @deftp {BinomialDistribution} {property} ParameterIsFixed
    ##
    ## Flag for fixed parameters
    ##
    ## A @math{1x2} logical vector specifying which parameters are fixed and
    ## which are estimated. @qcode{true} values correspond to fixed parameters,
    ## @qcode{false} values correspond to parameter estimates. This property is
    ## read-only.
    ##
    ## @end deftp
    ParameterIsFixed

    ## -*- texinfo -*-
    ## @deftp {BinomialDistribution} {property} Truncation
    ##
    ## Truncation interval
    ##
    ## A @math{1x2} numeric vector specifying the truncation interval for the
    ## probability distribution. First element contains the lower boundary,
    ## second element contains the upper boundary. This property is
    ## read-only. You can only truncate a probability distribution with the
    ## @qcode{truncate} method.
    ##
    ## @end deftp
    Truncation

    ## -*- texinfo -*-
    ## @deftp {BinomialDistribution} {property} IsTruncated
    ##
    ## Flag for truncated probability distribution
    ##
    ## A logical scalar value specifying whether a probability distribution is
    ## truncated or not. This property is read-only.
    ##
    ## @end deftp
    IsTruncated

    ## -*- texinfo -*-
    ## @deftp {BinomialDistribution} {property} InputData
    ##
    ## Data used for fitting a probability distribution
    ##
    ## A scalar structure containing the following fields:
    ## @itemize
    ## @item @qcode{data} : a numeric vector containing the data used for
    ## distribution fitting.
    ## @item @qcode{cens} : an empty array, since @qcode{BinomialDistribution}
    ## does not allow censoring.
    ## @item @qcode{frequency} : a numeric vector of non-negative integer values
    ## containing the frequency information corresponding to the elements of the
    ## data used for distribution fitting. If no frequency vector was used for
    ## distribution fitting, then this field defaults to an empty array.
    ##
    ## @end deftp
    InputData
  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    ParameterCI
  endproperties

  methods (Hidden)

    function this = BinomialDistribution (N, p)
      if (nargin == 0)
        N = 1;
        p = 0.5;
      endif
      checkparams (N, p);
      this.InputData = [];
      this.IsTruncated = false;
      this.ParameterValues = [N, p];
      this.ParameterIsFixed = [true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function display (this)
      fprintf ("%s =\n", inputname(1));
      __disp__ (this, "binomial distribution");
    endfunction

    function disp (this)
      __disp__ (this, "binomial distribution");
    endfunction

    function this = set.N (this, N)
      checkparams (N, this.p);
      this.InputData = [];
      this.ParameterValues(1) = N;
      this.ParameterIsFixed = [true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function N = get.N (this)
      N = this.ParameterValues(1);
    endfunction

    function this = set.p (this, p)
      checkparams (this.N, p);
      this.InputData = [];
      this.ParameterValues(2) = p;
      this.ParameterIsFixed = [true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function p = get.p (this)
      p = this.ParameterValues(2);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {BinomialDistribution} {@var{p} =} cdf (@var{pd}, @var{x})
    ## @deftypefnx {BinomialDistribution} {@var{p} =} cdf (@var{pd}, @var{x}, @qcode{"upper"})
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
      p = binocdf (x, this.N, this.p);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        p(lb) = 0;
        p(ub) = 1;
        p(! (lb | ub)) -= binocdf (lx - 1, this.N, this.p);
        p(! (lb | ub)) /= diff (binocdf ([lx-1, ux], this.N, this.p));
      endif
      ## Apply uflag
      if (utail)
        p = 1 - p;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {BinomialDistribution} {@var{x} =} icdf (@var{pd}, @var{p})
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
      umax = binoinv (1, this.N, this.p);
      if (this.IsTruncated && this.Truncation(2) >= umax)
        ## Find out of range p values
        is_nan = p < 0 | p > 1;
        ## Get lower and upper boundaries
        lx = ceil (this.Truncation(1));
        ux = floor (this.Truncation(2));
        ux = min (ux, umax);
        lp = binocdf (lx - 1, this.N, this.p);
        up = binocdf (ux, this.N, this.p);
        p = lp + p * (up - lp);
        p(is_nan) = NaN;
      endif
      x = binoinv (p, this.N, this.p);
      if (this.IsTruncated)
        x(x < this.Truncation(1)) = this.Truncation(1);
        x(x > this.Truncation(2)) = this.Truncation(2);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {BinomialDistribution} {@var{r} =} iqr (@var{pd})
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
    ## @deftypefn  {BinomialDistribution} {@var{m} =} mean (@var{pd})
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
      m = binostat (this.N, this.p);
      if (this.IsTruncated)
        lx = ceil (this.Truncation(1));
        ux = floor (this.Truncation(2));
        ux = min (ux, binoinv (1, this.N, this.p));
        x = [lx:ux];
        m = sum (x .* pdf (this, x));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {BinomialDistribution} {@var{m} =} median (@var{pd})
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
        Fa_b = binocdf ([lx, ux], this.N, this.p);
        m = binoinv (sum (Fa_b) / 2, this.N, this.p);
      else
        if (! __traditional__() & this.p == 0.5 & rem (this.N, 2) == 1)
          m = this.mean ();
        else
          m = binoinv (0.5, this.N, this.p);
        endif
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {BinomialDistribution} {@var{nlogL} =} negloglik (@var{pd})
    ##
    ## Compute the negative loglikelihood of a probability distribution.
    ##
    ## @code{@var{nlogL} = negloglik (@var{pd})} computes the negative loglikelihood
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
      nlogL = - binolike ([this.N, this.p], this.InputData.data);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {BinomialDistribution} {@var{ci} =} paramci (@var{pd})
    ## @deftypefnx {BinomialDistribution} {@var{ci} =} paramci (@var{pd}, @var{Name}, @var{Value})
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
    ## @deftypefn  {BinomialDistribution} {@var{y} =} pdf (@var{pd}, @var{x})
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
      y = binopdf (x, this.N, this.p);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        y(lb | ub) = 0;
        y(! (lb | ub)) /= diff (binocdf ([lx-1, ux], this.N, this.p));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {BinomialDistribution} {} plot (@var{pd})
    ## @deftypefnx {BinomialDistribution} {} plot (@var{pd}, @var{Name}, @var{Value})
    ## @deftypefnx {BinomialDistribution} {@var{h} =} plot (@dots{})
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
    ## @deftypefn  {BinomialDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum})
    ## @deftypefnx {BinomialDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @qcode{"Display"}, @var{display})
    ## @deftypefnx {BinomialDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @var{setparam})
    ## @deftypefnx {BinomialDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @var{setparam}, @qcode{"Display"}, @var{display})
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
    ## For the binomial distribution, @qcode{@var{pnum} = 1} selects the
    ## parameter @qcode{N} and @qcode{@var{pnum} = 2} selects the parameter
    ## @qcode{p}.
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
    ## @deftypefn  {BinomialDistribution} {@var{y} =} random (@var{pd})
    ## @deftypefnx {BinomialDistribution} {@var{y} =} random (@var{pd}, @var{rows})
    ## @deftypefnx {BinomialDistribution} {@var{y} =} random (@var{pd}, @var{rows}, @var{cols}, @dots{})
    ## @deftypefnx {BinomialDistribution} {@var{y} =} random (@var{pd}, [@var{sz}])
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
        ratio = 1 / diff (binocdf ([lx-1, ux], this.N, this.p));
        nsize = fix (2 * ratio * ps);       # times 2 to be on the safe side
        ## Generate the numbers and remove out-of-bound random samples
        r = binornd (this.N, this.p, nsize, 1);
        r(r < lx | r > ux) = [];
        ## Randomly select the required size and reshape to requested dimensions
        idx = randperm (numel (r), ps);
        r = reshape (r(idx), sz);
      else
        r = binornd (this.N, this.p, varargin{:});
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {BinomialDistribution} {@var{s} =} std (@var{pd})
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
    ## @deftypefn  {BinomialDistribution} {@var{t} =} truncate (@var{pd}, @var{lower}, @var{upper})
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
    ## @deftypefn  {BinomialDistribution} {@var{v} =} var (@var{pd})
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
        ## Calculate untruncated mean and variance
        [um, uv] = binostat (this.N, this.p);
        ## Calculate truncated mean
        m = mean (this);
        ## Get lower and upper boundaries
        lx = ceil (this.Truncation(1));
        ux = floor (this.Truncation(2));
        ux = min (ux, binoinv (1, this.N, this.p));
        ## Handle infinite support on the right
        if (isequal (ux, Inf))
          ratio = 1 / diff (binocdf ([lx-1, ux], this.N, this.p));
          x = 0:lx-1;
          v = ratio * (uv + (um - m) ^ 2 - sum (((x - m) .^ 2) .* ...
                                                binopdf (x, this.N, this.p)));
        else
          x = lx:ux;
          v = sum (((x - m) .^ 2) .* pdf (this, x));
        endif
      else
        [~, v] = binostat (this.N, this.p);
      endif
    endfunction

  endmethods

  methods (Static, Hidden)

    function pd = fit (x, ntrials, varargin)
      ## Check input arguments
      if (nargin < 3)
        alpha = 0.05;
      else
        alpha = varargin{1};
      endif
      if (nargin < 4)
        freq = ones (size (x));
      else
        freq = varargin{2};
      endif
      ## Fit data
      [pshat, psci] = mle (x, "distribution", "binomial", "alpha", alpha, ...
                              "ntrials", ntrials, "frequency", freq);
      phat = [ntrials, pshat];
      pci = [[ntrials; ntrials], psci(:)];
      [~, acov] = binolike (phat, x, freq);
      ## Create fitted distribution object
      pd = BinomialDistribution.makeFitted (phat, pci, acov, x, freq);
    endfunction

    function pd = makeFitted (phat, pci, acov, x, freq)
      N = phat(1);
      p = phat(2);
      pd = BinomialDistribution (N, p);
      pd.ParameterCI = pci;
      pd.ParameterIsFixed = [true, false];
      pd.ParameterCovariance = acov;
      pd.InputData = struct ("data", x, "cens", [], "freq", freq);
    endfunction

  endmethods

endclassdef

function checkparams (N, p)
  if (! (isscalar (N) && isnumeric (N) && isreal (N) && isfinite (N)
                      && N > 0 && fix (N) == N))
    error ("BinomialDistribution: N must be a positive integer scalar.")
  endif
  if (! (isscalar (p) && isnumeric (p) && isreal (p) && isfinite (p)
                      && p >= 0 && p <= 1))
    error (strcat (["BinomialDistribution: p must be a real"], ...
                   [" scalar bounded in the range [0, 1]."]))
  endif
endfunction

## Test output
%!shared pd, t, t_inf
%! pd = BinomialDistribution (5, 0.5);
%! t = truncate (pd, 2, 4);
%! t_inf = truncate (pd, 2, Inf);
%!assert (cdf (pd, [0:5]), [0.0312, 0.1875, 0.5, 0.8125, 0.9688, 1], 1e-4);
%!assert (cdf (t, [0:5]), [0, 0, 0.4, 0.8, 1, 1], 1e-4);
%!assert (cdf (t_inf, [0:5]), [0, 0, 0.3846, 0.7692, 0.9615, 1], 1e-4);
%!assert (cdf (pd, [1.5, 2, 3, 4, NaN]), [0.1875, 0.5, 0.8125, 0.9688, NaN], 1e-4);
%!assert (cdf (t, [1.5, 2, 3, 4, NaN]), [0, 0.4, 0.8, 1, NaN], 1e-4);
%!assert (icdf (pd, [0:0.2:1]), [0, 2, 2, 3, 3, 5], 1e-4);
%!assert (icdf (t, [0:0.2:1]), [2, 2, 2, 3, 3, 4], 1e-4);
%!assert (icdf (t_inf, [0:0.2:1]), [2, 2, 3, 3, 4, 5], 1e-4);
%!assert (icdf (pd, [-1, 0.4:0.2:1, NaN]), [NaN, 2, 3, 3, 5, NaN], 1e-4);
%!assert (icdf (t, [-1, 0.4:0.2:1, NaN]), [NaN, 2, 3, 3, 4, NaN], 1e-4);
%!assert (iqr (pd), 1);
%!assert (iqr (t), 1);
%!assert (mean (pd), 2.5, 1e-10);
%!assert (mean (t), 2.8, 1e-10);
%!assert (mean (t_inf), 2.8846, 1e-4);
%!assert (median (pd), 2.5);
%!assert (median (t), 3);
%!assert (pdf (pd, [0:5]), [0.0312, 0.1562, 0.3125, 0.3125, 0.1562, 0.0312], 1e-4);
%!assert (pdf (t, [0:5]), [0, 0, 0.4, 0.4, 0.2, 0], 1e-4);
%!assert (pdf (t_inf, [0:5]), [0, 0, 0.3846, 0.3846, 0.1923, 0.0385], 1e-4);
%!assert (pdf (pd, [-1, 1.5, NaN]), [0, 0, NaN], 1e-4);
%!assert (pdf (t, [-1, 1.5, NaN]), [0, 0, NaN], 1e-4);
%!assert (isequal (size (random (pd, 100, 50)), [100, 50]))
%!assert (any (random (t, 1000, 1) < 2), false);
%!assert (any (random (t, 1000, 1) > 4), false);
%!assert (std (pd), 1.1180, 1e-4);
%!assert (std (t), 0.7483, 1e-4);
%!assert (std (t_inf), 0.8470, 1e-4);
%!assert (var (pd), 1.2500, 1e-4);
%!assert (var (t), 0.5600, 1e-4);
%!assert (var (t_inf), 0.7175, 1e-4);

## Test input validation
## 'BinomialDistribution' constructor
%!error <BinomialDistribution: N must be a positive integer scalar.> ...
%! BinomialDistribution(Inf, 0.5)
%!error <BinomialDistribution: N must be a positive integer scalar.> ...
%! BinomialDistribution(i, 0.5)
%!error <BinomialDistribution: N must be a positive integer scalar.> ...
%! BinomialDistribution("a", 0.5)
%!error <BinomialDistribution: N must be a positive integer scalar.> ...
%! BinomialDistribution([1, 2], 0.5)
%!error <BinomialDistribution: N must be a positive integer scalar.> ...
%! BinomialDistribution(NaN, 0.5)
%!error <BinomialDistribution: p must be a real scalar bounded in the range> ...
%! BinomialDistribution(1, 1.01)
%!error <BinomialDistribution: p must be a real scalar bounded in the range> ...
%! BinomialDistribution(1, -0.01)
%!error <BinomialDistribution: p must be a real scalar bounded in the range> ...
%! BinomialDistribution(1, Inf)
%!error <BinomialDistribution: p must be a real scalar bounded in the range> ...
%! BinomialDistribution(1, i)
%!error <BinomialDistribution: p must be a real scalar bounded in the range> ...
%! BinomialDistribution(1, "a")
%!error <BinomialDistribution: p must be a real scalar bounded in the range> ...
%! BinomialDistribution(1, [1, 2])
%!error <BinomialDistribution: p must be a real scalar bounded in the range> ...
%! BinomialDistribution(1, NaN)

## 'cdf' method
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (BinomialDistribution, 2, "uper")
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (BinomialDistribution, 2, 3)

## 'paramci' method
%!shared x
%! rand ("seed", 2);
%! x = binornd (5, 0.5, [1, 100]);
%!error <paramci: optional arguments must be in NAME-VALUE pairs.> ...
%! paramci (BinomialDistribution.fit (x, 6), "alpha")
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (BinomialDistribution.fit (x, 6), "alpha", 0)
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (BinomialDistribution.fit (x, 6), "alpha", 1)
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (BinomialDistribution.fit (x, 6), "alpha", [0.5 2])
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (BinomialDistribution.fit (x, 6), "alpha", "")
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (BinomialDistribution.fit (x, 6), "alpha", {0.05})
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (BinomialDistribution.fit (x, 6), "parameter", "p", ...
%!          "alpha", {0.05})
%!error <paramci: invalid VALUE size for 'Parameter' argument.> ...
%! paramci (BinomialDistribution.fit (x, 6), ...
%!          "parameter", {"N", "p", "param"})
%!error <paramci: invalid VALUE size for 'Parameter' argument.> ...
%! paramci (BinomialDistribution.fit (x, 6), "alpha", 0.01, ...
%!          "parameter", {"N", "p", "param"})
%!error <paramci: unknown distribution parameter.> ...
%! paramci (BinomialDistribution.fit (x, 6), "parameter", "param")
%!error <paramci: unknown distribution parameter.> ...
%! paramci (BinomialDistribution.fit (x, 6), "parameter", "N")
%!error <paramci: unknown distribution parameter.> ...
%! paramci (BinomialDistribution.fit (x, 6), "alpha", 0.01, ...
%!          "parameter", "param")
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (BinomialDistribution.fit (x, 6), "NAME", "value")
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (BinomialDistribution.fit (x, 6), "alpha", 0.01, ...
%!          "NAME", "value")
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (BinomialDistribution.fit (x, 6), "alpha", 0.01, ...
%!          "parameter", "p", "NAME", "value")

## 'plot' method
%!error <plot: optional arguments must be in NAME-VALUE pairs.> ...
%! plot (BinomialDistribution, "Parent")
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (BinomialDistribution, "PlotType", 12)
%!error <plot: invalid VALUE size for 'Parameter' argument.> ...
%! plot (BinomialDistribution, "PlotType", {"pdf", "cdf"})
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (BinomialDistribution, "PlotType", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (BinomialDistribution, "Discrete", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (BinomialDistribution, "Discrete", [1, 0])
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (BinomialDistribution, "Discrete", {true})
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (BinomialDistribution, "Parent", 12)
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (BinomialDistribution, "Parent", "hax")
%!error <plot: invalid NAME for optional argument.> ...
%! plot (BinomialDistribution, "invalidNAME", "pdf")
%!error <plot: no fitted DATA to plot a probability plot.> ...
%! plot (BinomialDistribution, "PlotType", "probability")

## 'proflik' method
%!error <proflik: no fitted data available.> ...
%! proflik (BinomialDistribution, 2)
%!error <proflik: PNUM must be a scalar number indexing a non-fixed parameter.> ...
%! proflik (BinomialDistribution.fit (x, 6), 3)
%!error <proflik: PNUM must be a scalar number indexing a non-fixed parameter.> ...
%! proflik (BinomialDistribution.fit (x, 6), [1, 2])
%!error <proflik: PNUM must be a scalar number indexing a non-fixed parameter.> ...
%! proflik (BinomialDistribution.fit (x, 6), {1})
%!error <proflik: SETPARAM must be a numeric vector.> ...
%! proflik (BinomialDistribution.fit (x, 6), 2, ones (2))
%!error <proflik: missing VALUE for 'Display' argument.> ...
%! proflik (BinomialDistribution.fit (x, 6), 2, "Display")
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (BinomialDistribution.fit (x, 6), 2, "Display", 1)
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (BinomialDistribution.fit (x, 6), 2, "Display", {1})
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (BinomialDistribution.fit (x, 6), 2, "Display", {"on"})
%!error <proflik: invalid VALUE size for 'Display' argument.> ...
%! proflik (BinomialDistribution.fit (x, 6), 2, "Display", ["on"; "on"])
%!error <proflik: invalid VALUE for 'Display' argument.> ...
%! proflik (BinomialDistribution.fit (x, 6), 2, "Display", "onnn")
%!error <proflik: invalid NAME for optional arguments.> ...
%! proflik (BinomialDistribution.fit (x, 6), 2, "NAME", "on")
%!error <proflik: invalid optional argument.> ...
%! proflik (BinomialDistribution.fit (x, 6), 2, {"NAME"}, "on")
%!error <proflik: invalid optional argument.> ...
%! proflik (BinomialDistribution.fit (x, 6), 2, {[1 2 3]}, "Display", "on")

## 'truncate' method
%!error <truncate: missing input argument.> ...
%! truncate (BinomialDistribution)
%!error <truncate: missing input argument.> ...
%! truncate (BinomialDistribution, 2)
%!error <truncate: invalid lower upper limits.> ...
%! truncate (BinomialDistribution, 4, 2)

## Catch errors when using array of probability objects with available methods
%!shared pd
%! pd = BinomialDistribution(1, 0.5);
%! pd(2) = BinomialDistribution(1, 0.6);
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
