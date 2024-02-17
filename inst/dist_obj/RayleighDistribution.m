## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

classdef RayleighDistribution
  ## -*- texinfo -*-
  ## @deftypefn {statistics} RayleighDistribution
  ##
  ## Rayleigh probability distribution object.
  ##
  ## A @code{RayleighDistribution} object consists of parameters, a model
  ## description, and sample data for a Rayleigh probability distribution.
  ##
  ## The Rayleigh distribution uses the following parameters.
  ##
  ## @multitable @columnfractions 0.25 0.48 0.27
  ## @headitem @var{Parameter} @tab @var{Description} @tab @var{Support}
  ##
  ## @item @qcode{sigma} @tab Scale parameter @tab @math{sigma > 0}
  ## @end multitable
  ##
  ## There are several ways to create a @code{RayleighDistribution} object.
  ##
  ## @itemize
  ## @item Fit a distribution to data using the @code{fitdist} function.
  ## @item Create a distribution with specified parameter values using the
  ## @code{makedist} function.
  ## @item Use the constructor @qcode{RayleighDistribution (@var{sigma})}
  ## to create a Rayleigh distribution with specified parameter values.
  ## @item Use the static method @qcode{RayleighDistribution.fit (@var{x},
  ## @var{censor}, @var{freq}, @var{options})} to a distribution to data @var{x}.
  ## @end itemize
  ##
  ## It is highly recommended to use @code{fitdist} and @code{makedist}
  ## functions to create probability distribution objects, instead of the
  ## constructor and the aforementioned static method.
  ##
  ## A @code{RayleighDistribution} object contains the following properties,
  ## which can be accessed using dot notation.
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
  ## A @code{RayleighDistribution} object contains the following methods:
  ## @code{cdf}, @code{icdf}, @code{iqr}, @code{mean}, @code{median},
  ## @code{negloglik}, @code{paramci}, @code{pdf}, @code{plot}, @code{proflik},
  ## @code{random}, @code{std}, @code{truncate}, @code{var}.
  ##
  ## Further information about the Rayleigh distribution can be found at
  ## @url{https://en.wikipedia.org/wiki/Rayleigh_distribution}
  ##
  ## @seealso{fitdist, makedist, raylcdf, raylinv, raylpdf, raylrnd, raylfit,
  ## rayllike, raylstat}
  ## @end deftypefn

  properties (Dependent = true)
    sigma
  endproperties

  properties (GetAccess = public, Constant = true)
    CensoringAllowed = true;
    DistributionName = "RayleighDistribution";
    DistributionCode = "rayl";
    NumParameters = 1;
    ParameterNames = {"sigma"};
    ParameterDescription = {"Scale"};
  endproperties

  properties (GetAccess = private, Constant = true)
    ParameterRange = [realmin; Inf];
    ParameterLogCI = true;
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

    function this = RayleighDistribution (sigma)
      if (nargin == 0)
        sigma = 1;
      endif
      checkparams (sigma)
      this.InputData = [];
      this.IsTruncated = false;
      this.ParameterValues = sigma;
      this.ParameterIsFixed = true;
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function display (this)
      fprintf ("%s =\n", inputname(1));
      __disp__ (this, "Rayleigh distribution");
    endfunction

    function disp (this)
      __disp__ (this, "Rayleigh distribution");
    endfunction

    function this = set.sigma (this, sigma)
      checkparams (sigma)
      this.InputData = [];
      this.ParameterValues(1) = sigma;
      this.ParameterIsFixed = true;
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function sigma = get.sigma (this)
      sigma = this.ParameterValues(1);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {RayleighDistribution} {@var{p} =} cdf (@var{pd}, @var{x})
    ## @deftypefnx {RayleighDistribution} {@var{p} =} cdf (@var{pd}, @var{x}, @qcode{"upper"})
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
      p = raylcdf (x, this.sigma);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        p(lb) = 0;
        p(ub) = 1;
        p(! (lb | ub)) -= raylcdf (lx, this.sigma);
        p(! (lb | ub)) /= diff (raylcdf ([lx, ux], this.sigma));
      endif
      ## Apply uflag
      if (utail)
        p = 1 - p;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RayleighDistribution} {@var{p} =} icdf (@var{pd}, @var{p})
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
        lp = raylcdf (this.Truncation(1), this.sigma);
        up = raylcdf (this.Truncation(2), this.sigma);
        ## Adjust p values within range of p @ lower limit and p @ upper limit
        np = lp + (up - lp) .* p;
        x = raylinv (np, this.sigma);
      else
        x = raylinv (p, this.sigma);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RayleighDistribution} {@var{r} =} iqr (@var{pd})
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
    ## @deftypefn  {RayleighDistribution} {@var{m} =} mean (@var{pd})
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
        m = raylstat (this.sigma);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RayleighDistribution} {@var{m} =} median (@var{pd})
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
        Fa_b = raylcdf ([lx, ux], this.sigma);
        m = raylinv (sum (Fa_b) / 2, this.sigma);
      else
        m = this.sigma .* sqrt (2 * log (2));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RayleighDistribution} {@var{nlogL} =} negloglik (@var{pd})
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
      nlogL = - rayllike (this.sigma, this.InputData.data, ...
                          this.InputData.cens, this.InputData.freq);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RayleighDistribution} {@var{ci} =} paramci (@var{pd})
    ## @deftypefnx {RayleighDistribution} {@var{ci} =} paramci (@var{pd}, @var{Name}, @var{Value})
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
    ## @deftypefn  {RayleighDistribution} {@var{y} =} pdf (@var{pd}, @var{x})
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
      y = raylpdf (x, this.sigma);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        y(lb | ub) = 0;
        y(! (lb | ub)) /= diff (raylcdf ([lx, ux], this.sigma));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RayleighDistribution} {} plot (@var{pd})
    ## @deftypefnx {RayleighDistribution} {} plot (@var{pd}, @var{Name}, @var{Value})
    ## @deftypefnx {RayleighDistribution} {@var{h} =} plot (@dots{})
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
    ## @deftypefn  {RayleighDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum})
    ## @deftypefnx {RayleighDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @qcode{"Display"}, @var{display})
    ## @deftypefnx {RayleighDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @var{setparam})
    ## @deftypefnx {RayleighDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @var{setparam}, @qcode{"Display"}, @var{display})
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
    ## For the Rayleigh distribution, @qcode{@var{pnum} = 1} selects the
    ## parameter @var{sigma}.
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
    ## @deftypefn  {RayleighDistribution} {@var{y} =} random (@var{pd})
    ## @deftypefnx {RayleighDistribution} {@var{y} =} random (@var{pd}, @var{rows})
    ## @deftypefnx {RayleighDistribution} {@var{y} =} random (@var{pd}, @var{rows}, @var{cols}, @dots{})
    ## @deftypefnx {RayleighDistribution} {@var{y} =} random (@var{pd}, [@var{sz}])
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
        lp = raylcdf (this.Truncation(1), this.sigma);
        up = raylcdf (this.Truncation(2), this.sigma);
        u = unifrnd (lp, up, varargin{:});
        r = sqrt (-2 .* log (1 - u) .* sigma .^ 2);
      else
        r = raylrnd (this.sigma, varargin{:});
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RayleighDistribution} {@var{s} =} std (@var{pd})
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
      v = var (this.sigma);
      s = sqrt (v);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RayleighDistribution} {@var{t} =} truncate (@var{pd}, @var{lower}, @var{upper})
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
      ## Check boundaries and constrain within support [0, Inf)
      lower(lower < 0) = 0;
      this.Truncation = [lower, upper];
      this.IsTruncated = true;
      this.InputData = [];
      this.ParameterIsFixed = true;
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RayleighDistribution} {@var{v} =} var (@var{pd})
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
        [~, v] = raylstat (this.sigma);
      endif
    endfunction

  endmethods

  methods (Static, Hidden)

    function pd = fit (x, varargin)
      ## Check input arguments
      if (nargin < 2)
        censor = [];
      endif
      if (nargin < 3)
        freq = [];
      endif
      ## Fit data
      [phat, pci] = raylfit (x, 0.05, censor, freq);
      [~, acov] = rayllike (phat, x, censor, freq);
      ## Create fitted distribution object
      pd = RayleighDistribution.makeFitted ...
           (phat, pci, acov, x, censor, freq);
    endfunction

    function pd = makeFitted (phat, pci, acov, x, censor, freq)
      sigma = phat(1);
      pd = RayleighDistribution (sigma);
      pd.ParameterCI = pci;
      pd.ParameterIsFixed = false;
      pd.ParameterCovariance = acov;
      pd.InputData = struct ("data", x, "cens", censor, "freq", freq);
    endfunction

  endmethods

endclassdef

function checkparams (sigma)
  if (! (isscalar (sigma) && isnumeric (sigma) && isreal (sigma)
                          && isfinite (sigma) && sigma > 0))
    error ("RayleighDistribution: SIGMA must be a positive real scalar.")
  endif
endfunction

## Test input validation
## 'RayleighDistribution' constructor
%!error <RayleighDistribution: SIGMA must be a positive real scalar.> ...
%! RayleighDistribution(0)
%!error <RayleighDistribution: SIGMA must be a positive real scalar.> ...
%! RayleighDistribution(-1)
%!error <RayleighDistribution: SIGMA must be a positive real scalar.> ...
%! RayleighDistribution(Inf)
%!error <RayleighDistribution: SIGMA must be a positive real scalar.> ...
%! RayleighDistribution(i)
%!error <RayleighDistribution: SIGMA must be a positive real scalar.> ...
%! RayleighDistribution("a")
%!error <RayleighDistribution: SIGMA must be a positive real scalar.> ...
%! RayleighDistribution([1, 2])
%!error <RayleighDistribution: SIGMA must be a positive real scalar.> ...
%! RayleighDistribution(NaN)

## 'cdf' method
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (RayleighDistribution, 2, "uper")
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (RayleighDistribution, 2, 3)

## 'paramci' method
%!error <paramci: optional arguments must be in NAME-VALUE pairs.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), "alpha")
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), "alpha", 0)
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), "alpha", 1)
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), "alpha", [0.5 2])
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), "alpha", "")
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), "alpha", {0.05})
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), ...
%! "parameter", "sigma", "alpha", {0.05})
%!error <paramci: invalid VALUE size for 'Parameter' argument.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), ...
%! "parameter", {"sigma", "param"})
%!error <paramci: invalid VALUE size for 'Parameter' argument.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), ...
%! "alpha", 0.01, "parameter", {"sigma", "param"})
%!error <paramci: unknown distribution parameter.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), ...
%! "parameter", "param")
%!error <paramci: unknown distribution parameter.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), ...
%! "alpha", 0.01, "parameter", "param")
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), ...
%! "NAME", "value")
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), ...
%! "alpha", 0.01, "NAME", "value")
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (RayleighDistribution.fit (raylrnd (1, [1, 100])), ...
%! "alpha", 0.01, "parameter", "sigma", "NAME", "value")

## 'plot' method
%!error <plot: optional arguments must be in NAME-VALUE pairs.> ...
%! plot (RayleighDistribution, "Parent")
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (RayleighDistribution, "PlotType", 12)
%!error <plot: invalid VALUE size for 'Parameter' argument.> ...
%! plot (RayleighDistribution, "PlotType", {"pdf", "cdf"})
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (RayleighDistribution, "PlotType", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (RayleighDistribution, "Discrete", "pdfcdf")
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (RayleighDistribution, "Discrete", [1, 0])
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (RayleighDistribution, "Discrete", {true})
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (RayleighDistribution, "Parent", 12)
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (RayleighDistribution, "Parent", "hax")

## 'proflik' method
%!error <proflik: no fitted data available.> ...
%! proflik (RayleighDistribution, 2)
%!error <proflik: PNUM must be a scalar number indexing a valid parameter.> ...
%! proflik (RayleighDistribution.fit (raylrnd (1, [1, 100])), 3)
%!error <proflik: PNUM must be a scalar number indexing a valid parameter.> ...
%! proflik (RayleighDistribution.fit (raylrnd (1, [1, 100])), [1, 2])
%!error <proflik: PNUM must be a scalar number indexing a valid parameter.> ...
%! proflik (RayleighDistribution.fit (raylrnd (1, [1, 100])), {1})
%!error <proflik: SETPARAM must be a numeric vector.> ...
%! proflik (RayleighDistribution.fit (raylrnd (1, [1, 100])), 1, ones (2))
%!error <proflik: missing VALUE for 'Display' argument.> ...
%! proflik (RayleighDistribution.fit (raylrnd (1, [1, 100])), 1, ...
%! "Display")
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (RayleighDistribution.fit (raylrnd (1, [1, 100])), 1, ...
%! "Display", 1)
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (RayleighDistribution.fit (raylrnd (1, [1, 100])), 1, ...
%! "Display", {1})
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (RayleighDistribution.fit (raylrnd (1, [1, 100])), 1, ...
%! "Display", {"on"})
%!error <proflik: invalid VALUE size for 'Display' argument.> ...
%! proflik (RayleighDistribution.fit (raylrnd (1, [1, 100])), 1, ...
%! "Display", ["on"; "on"])
%!error <proflik: invalid VALUE for 'Display' argument.> ...
%! proflik (RayleighDistribution.fit (raylrnd (1, [1, 100])), 1, ...
%! "Display", "onnn")
%!error <proflik: invalid NAME for optional arguments.> ...
%! proflik (RayleighDistribution.fit (raylrnd (1, [1, 100])), 1, ...
%! "NAME", "on")
%!error <proflik: invalid optional argument.> ...
%! proflik (RayleighDistribution.fit (raylrnd (1, [1, 100])), 1, ...
%! {"NAME"}, "on")
%!error <proflik: invalid optional argument.> ...
%! proflik (RayleighDistribution.fit (raylrnd (1, [1, 100])), 1, ...
%! {[1 2 3 4]}, "Display", "on")

## 'truncate' method
%!error <truncate: missing input argument.> ...
%! truncate (RayleighDistribution)
%!error <truncate: missing input argument.> ...
%! truncate (RayleighDistribution, 2)
%!error <truncate: invalid lower upper limits.> ...
%! truncate (RayleighDistribution, 4, 2)

## Catch errors when using array of probability objects with available methods
%!shared pd
%! pd = RayleighDistribution(1);
%! pd(2) = RayleighDistribution(3);
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
