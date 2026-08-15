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

classdef GeneralizedParetoDistribution < prob.ProbabilityDistribution
  ## -*- texinfo -*-
  ## @deftp {statistics} prob.GeneralizedParetoDistribution
  ##
  ## Generalized Pareto probability distribution object.
  ##
  ## A @code{prob.GeneralizedParetoDistribution} object consists of parameters, a
  ## model description, and sample data for a Generalized Pareto probability
  ## distribution.
  ##
  ## The Generalized Pareto distribution is a continuous probability
  ## distribution that models the tail behavior of other distributions, commonly
  ## used for extreme value analysis.  It is defined by shape parameter @var{k},
  ## scale parameter @var{sigma}, and location parameter @var{theta}.
  ##
  ## There are several ways to create a @code{prob.GeneralizedParetoDistribution}
  ## object.
  ##
  ## @itemize
  ## @item Fit a distribution to data using the @code{fitdist} function.
  ## @item Create a distribution with fixed parameter values using the
  ## @code{makedist} function.
  ## @item Use the constructor @qcode{prob.GeneralizedParetoDistribution (@var{k},
  ## @var{sigma}, @var{theta})} to create a Generalized Pareto distribution with
  ## fixed parameter values @var{k}, @var{sigma}, and @var{theta}.
  ## @item Use the static method @qcode{prob.GeneralizedParetoDistribution.fit
  ## (@var{x}, @var{theta}, @var{alpha}, @var{freq}, @var{options})} to fit a
  ## distribution to the data in @var{x} using the same input arguments as the
  ## @code{gpfit} function.
  ## @end itemize
  ##
  ## It is highly recommended to use @code{fitdist} and @code{makedist}
  ## functions to create probability distribution objects, instead of the class
  ## constructor or the aforementioned static method.
  ##
  ## Further information about the Generalized Pareto distribution can be found
  ## at
  ## @url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
  ##
  ## @seealso{fitdist, makedist, gpcdf, gpinv, gppdf, gprnd, gpfit,
  ## gplike, gpstat}
  ## @end deftp

  properties(Dependent = true)
    ## -*- texinfo -*-
    ## @deftp {prob.GeneralizedParetoDistribution} {property} k
    ##
    ## Shape parameter
    ##
    ## A scalar value characterizing the shape of the Generalized Pareto
    ## distribution. You can access the @qcode{k} property using dot name
    ## assignment.
    ##
    ## @end deftp
    k

    ## -*- texinfo -*-
    ## @deftp {prob.GeneralizedParetoDistribution} {property} sigma
    ##
    ## Scale parameter
    ##
    ## A positive scalar value characterizing the scale of the Generalized
    ## Pareto distribution. You can access the @qcode{sigma} property using dot
    ## name assignment.
    ##
    ## @end deftp
    sigma

    ## -*- texinfo -*-
    ## @deftp {prob.GeneralizedParetoDistribution} {property} theta
    ##
    ## Location parameter
    ##
    ## A scalar value characterizing the location of the Generalized Pareto
    ## distribution. You can access the @qcode{theta} property using dot name
    ## assignment.
    ##
    ## @end deftp
    theta
  endproperties

  properties(GetAccess = public, Constant = true)
    ## -*- texinfo -*-
    ## @deftp {prob.GeneralizedParetoDistribution} {property} DistributionName
    ##
    ## Probability distribution name
    ##
    ## A character vector specifying the name of the probability distribution
    ## object. This property is read-only.
    ##
    ## @end deftp
    DistributionName = 'Generalized Pareto';

    ## -*- texinfo -*-
    ## @deftp {prob.GeneralizedParetoDistribution} {property} NumParameters
    ##
    ## Number of parameters
    ##
    ## A scalar integer value specifying the number of parameters characterizing
    ## the probability distribution. This property is read-only.
    ##
    ## @end deftp
    NumParameters = 3;

    ## -*- texinfo -*-
    ## @deftp {prob.GeneralizedParetoDistribution} {property} ParameterNames
    ##
    ## Names of parameters
    ##
    ## A @math{3*1} cell array of character vectors with each element containing
    ## the name of a distribution parameter. This property is read-only.
    ##
    ## @end deftp
    ParameterNames = {'k', 'sigma', 'theta'};

    ## -*- texinfo -*-
    ## @deftp {prob.GeneralizedParetoDistribution} {property} ParameterDescription
    ##
    ## Description of parameters
    ##
    ## A @math{3*1} cell array of character vectors with each element containing
    ## a short description of a distribution parameter. This property is
    ## read-only.
    ##
    ## @end deftp
    ParameterDescription = {'Shape', 'Scale', 'Location'};
  endproperties

  properties(GetAccess = public, Constant = true, Hidden)
    CensoringAllowed = false;
    DistributionCode = 'gp';
    ParameterRange = [-Inf, realmin, -Inf; Inf, Inf, Inf];
    ParameterLogCI = [false, true, false];
  endproperties

  properties(GetAccess = public, SetAccess = protected)
    ## -*- texinfo -*-
    ## @deftp {prob.GeneralizedParetoDistribution} {property} ParameterValues
    ##
    ## Distribution parameter values
    ##
    ## A @math{3*1} numeric vector containing the values of the distribution
    ## parameters. This property is read-only. You can change the distribution
    ## parameters by assigning new values to the @qcode{k}, @qcode{sigma}, and
    ## @qcode{theta} properties.
    ##
    ## @end deftp
    ParameterValues

    ## -*- texinfo -*-
    ## @deftp {prob.GeneralizedParetoDistribution} {property} ParameterCovariance
    ##
    ## Covariance matrix of the parameter estimates
    ##
    ## A @math{3*3} numeric matrix containing the variance-covariance of the
    ## parameter estimates. Diagonal elements contain the variance of each
    ## estimated parameter, and non-diagonal elements contain the covariance
    ## between the parameter estimates. The covariance matrix is only meaningful
    ## when the distribution was fitted to data. If the distribution object was
    ## created with fixed parameters, or a parameter of a fitted distribution is
    ## modified, then all elements of the variance-covariance are zero. This
    ## property is read-only.
    ##
    ## @end deftp
    ParameterCovariance

    ## -*- texinfo -*-
    ## @deftp {prob.GeneralizedParetoDistribution} {property} ParameterIsFixed
    ##
    ## Flag for fixed parameters
    ##
    ## A @math{1*3} logical vector specifying which parameters are fixed and
    ## which are estimated. @qcode{true} values correspond to fixed parameters,
    ## @qcode{false} values correspond to parameter estimates. This property is
    ## read-only.
    ##
    ## @end deftp
    ParameterIsFixed

    ## -*- texinfo -*-
    ## @deftp {prob.GeneralizedParetoDistribution} {property} Truncation
    ##
    ## Truncation interval
    ##
    ## A @math{1*2} numeric vector specifying the truncation interval for the
    ## probability distribution. First element contains the lower boundary,
    ## second element contains the upper boundary. This property is read-only.
    ## You can only truncate a probability distribution with the
    ## @qcode{truncate} method.
    ##
    ## @end deftp
    Truncation

    ## -*- texinfo -*-
    ## @deftp {prob.GeneralizedParetoDistribution} {property} IsTruncated
    ##
    ## Flag for truncated probability distribution
    ##
    ## A logical scalar value specifying whether a probability distribution is
    ## truncated or not. This property is read-only.
    ##
    ## @end deftp
    IsTruncated

    ## -*- texinfo -*-
    ## @deftp {prob.GeneralizedParetoDistribution} {property} InputData
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

  properties(GetAccess = public, SetAccess = protected, Hidden)
    ParameterCI
  endproperties

  methods(Hidden)

    function this = GeneralizedParetoDistribution (k, sigma, theta)
      if (nargin == 0)
        k = 1;
        sigma = 1;
        theta = 1;
      endif
      checkparams (k, sigma, theta);
      this.InputData = [];
      this.IsTruncated = false;
      this.ParameterValues = [k, sigma, theta];
      this.ParameterIsFixed = [true, true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function display (this)
      fprintf ("%s =\n", inputname (1));
      __disp__ (this, 'Generalized Pareto distribution');
    endfunction

    function disp (this)
      __disp__ (this, 'Generalized Pareto distribution');
    endfunction

    function this = set.k (this, k)
      checkparams (k, this.sigma, this.theta);
      this.InputData = [];
      this.ParameterValues(1) = k;
      this.ParameterIsFixed = [true, true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function k = get.k (this)
      k = this.ParameterValues(1);
    endfunction

    function this = set.sigma (this, sigma)
      checkparams (this.k, sigma, this.theta);
      this.InputData = [];
      this.ParameterValues(2) = sigma;
      this.ParameterIsFixed = [true, true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function sigma = get.sigma (this)
      sigma = this.ParameterValues(2);
    endfunction

    function this = set.theta (this, theta)
      checkparams (this.k, this.sigma, theta);
      this.InputData = [];
      this.ParameterValues(3) = theta;
      this.ParameterIsFixed = [true, true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    function theta = get.theta (this)
      theta = this.ParameterValues(3);
    endfunction

  endmethods

  methods(Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {@var{p} =} cdf (@var{pd}, @var{x})
    ## @deftypefnx {prob.GeneralizedParetoDistribution} {@var{p} =} cdf (@var{pd}, @var{x}, @qcode{'upper'})
    ##
    ## Compute the cumulative distribution function (CDF).
    ##
    ## @code{@var{p} = cdf (@var{pd}, @var{x})} computes the CDF of the
    ## probability distribution object, @var{pd}, evaluated at the values in
    ## @var{x}.
    ##
    ## @code{@var{p} = cdf (@dots{}, @qcode{'upper'})} returns the complement of
    ## the CDF of the probability distribution object, @var{pd}, evaluated at
    ## the values in @var{x}.
    ##
    ## @end deftypefn
    function p = cdf (this, x, uflag)
      if (! isscalar (this))
        error ("cdf: requires a scalar probability distribution.");
      endif
      ## Check for "upper" flag
      if (nargin > 2 && strcmpi (uflag, 'upper'))
        utail = true;
      elseif (nargin > 2 && ! strcmpi (uflag, 'upper'))
        error ("cdf: invalid argument for upper tail.");
      else
        utail = false;
      endif
      ## Do the computations
      p = gpcdf (x, this.k, this.sigma, this.theta);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        p(lb) = 0;
        p(ub) = 1;
        p(! (lb | ub)) -= gpcdf (lx, this.k, this.sigma, this.theta);
        p(! (lb | ub)) /= diff (gpcdf ([lx, ux], this.k, this.sigma, this.theta));
      endif
      ## Apply uflag
      if (utail)
        p = 1 - p;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {@var{x} =} icdf (@var{pd}, @var{p})
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
        lp = gpcdf (this.Truncation(1), this.k, this.sigma, this.theta);
        up = gpcdf (this.Truncation(2), this.k, this.sigma, this.theta);
        ## Adjust p values within range of p @ lower limit and p @ upper limit
        is_nan = p < 0 | p > 1;
        p(is_nan) = NaN;
        np = lp + (up - lp) .* p;
        x = gpinv (np, this.k, this.sigma, this.theta);
        x(x < this.Truncation(1)) = this.Truncation(1);
        x(x > this.Truncation(2)) = this.Truncation(2);
      else
        x = gpinv (p, this.k, this.sigma, this.theta);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {@var{r} =} iqr (@var{pd})
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
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {@var{m} =} mean (@var{pd})
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
        m = gpstat (this.k, this.sigma, this.theta);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {@var{m} =} median (@var{pd})
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
        Fa_b = gpcdf ([lx, ux], this.k, this.sigma, this.theta);
        m = gpinv (sum (Fa_b) / 2, this.k, this.sigma, this.theta);
      else
        m = gpinv (0.5, this.k, this.sigma, this.theta);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {@var{nlogL} =} negloglik (@var{pd})
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
      nlogL = gplike ([this.k, this.sigma, this.theta], ...
                        this.InputData.data, this.InputData.freq);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {@var{ci} =} paramci (@var{pd})
    ## @deftypefnx {prob.GeneralizedParetoDistribution} {@var{ci} =} paramci (@var{pd}, @var{Name}, @var{Value})
    ##
    ## Compute the confidence intervals for probability distribution parameters.
    ##
    ## @code{@var{ci} = paramci (@var{pd})} computes the lower and upper
    ## boundaries of the 95% confidence interval for each parameter of the
    ## probability distribution object, @var{pd}.
    ##
    ## @code{@var{ci} = paramci (@var{pd}, @var{Name}, @var{Value})} computes
    ## the
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
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {@var{y} =} pdf (@var{pd}, @var{x})
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
      y = gppdf (x, this.k, this.sigma, this.theta);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        y(lb | ub) = 0;
        y(! (lb | ub)) /= diff (gpcdf ([lx, ux], this.k, this.sigma, this.theta));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {} plot (@var{pd})
    ## @deftypefnx {prob.GeneralizedParetoDistribution} {} plot (@var{pd}, @var{Name}, @var{Value})
    ## @deftypefnx {prob.GeneralizedParetoDistribution} {@var{h} =} plot (@dots{})
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
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'PlotType'} @tab A character vector specifying the plot
    ## type.  @qcode{'pdf'} plots the probability density function (PDF).  When
    ## @var{pd} is fit to data, the PDF is superimposed on a histogram of the
    ## data.  @qcode{'cdf'} plots the cumulative density function (CDF).  When
    ## @var{pd} is fit to data, the CDF is superimposed over an empirical CDF.
    ## @qcode{'probability'} plots a probability plot using a CDF of the data
    ## and a CDF of the fitted probability distribution.  This option is
    ## available only when @var{pd} is fitted to data.
    ##
    ## @item @qcode{'Discrete'} @tab A logical scalar to specify whether to
    ## plot the PDF or CDF of a discrete distribution object as a line plot or a
    ## stem plot, by specifying @qcode{false} or @qcode{true}, respectively.  By
    ## default, it is @qcode{true} for discrete distributions and @qcode{false}
    ## for continuous distributions.  When @var{pd} is a continuous distribution
    ## object, option is ignored.
    ##
    ## @item @qcode{'Parent'} @tab An axes graphics object for plot.  If
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
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum})
    ## @deftypefnx {prob.GeneralizedParetoDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @qcode{'Display'}, @var{display})
    ## @deftypefnx {prob.GeneralizedParetoDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @var{setparam})
    ## @deftypefnx {prob.GeneralizedParetoDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd}, @var{pnum}, @var{setparam}, @qcode{'Display'}, @var{display})
    ## @deftypefnx {prob.GeneralizedParetoDistribution} {[@var{nlogL}, @var{param}] =} proflik (@var{pd})
    ## @deftypefnx {prob.GeneralizedParetoDistribution} {[@var{nlogL}, @var{param}, @var{other}] =} proflik (@dots{})
    ##
    ## Profile likelihood function for a probability distribution object.
    ##
    ## @code{[@var{nlogL}, @var{param}] = proflik (@var{pd}, @var{pnum})}
    ## returns a vector @var{nlogL} of negative loglikelihood values and a
    ## vector @var{param} of corresponding parameter values for the parameter in
    ## the position indicated by @var{pnum}.  By default, @code{proflik} uses
    ## the lower and upper bounds of the 98% confidence interval and computes
    ## 101 equispaced values for the selected parameter when it is the only one
    ## being estimated, and 21 values otherwise.  @var{pd} must be fitted to
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
    ## @code{[@var{nlogL}, @var{param}] = proflik (@var{pd})} selects the
    ## first parameter that is not fixed.
    ##
    ## @code{[@var{nlogL}, @var{param}, @var{other}] = proflik (@dots{})} also
    ## returns a matrix @var{other} holding, in each row, the values of the
    ## remaining parameters that maximize the likelihood at the corresponding
    ## value of @var{param}.  A fixed parameter keeps its own value.
    ##
    ## For the Generalized Pareto distribution, @qcode{@var{pnum} = 1} selects
    ## the parameter @qcode{k}, @qcode{@var{pnum} = 2} selects the parameter
    ## @qcode{sigma}, and @qcode{@var{pnum} = 3} selects the parameter
    ## @qcode{theta}.
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
      if (nargin < 2)
        pnum = [];
      endif
      [varargout{1:nargout}] = __proflik__ (this, pnum, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {@var{r} =} random (@var{pd})
    ## @deftypefnx {prob.GeneralizedParetoDistribution} {@var{r} =} random (@var{pd}, @var{rows})
    ## @deftypefnx {prob.GeneralizedParetoDistribution} {@var{r} =} random (@var{pd}, @var{rows}, @var{cols}, @dots{})
    ## @deftypefnx {prob.GeneralizedParetoDistribution} {@var{r} =} random (@var{pd}, [@var{sz}])
    ##
    ## Generate random arrays from the probability distribution object.
    ##
    ## @code{@var{r} = random (@var{pd})} returns a random number from the
    ## distribution object @var{pd}.
    ##
    ## When called with a single size argument, @code{random} returns a square
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
        ratio = 1 / diff (gpcdf ([lx, ux], this.k, this.sigma, this.theta));
        nsize = fix (2 * ratio * ps);       # times 2 to be on the safe side
        ## Generate the numbers and remove out-of-bound random samples
        r = gprnd (this.k, this.sigma, this.theta, nsize, 1);
        r(r < lx | r > ux) = [];
        ## Randomly select the required size and reshape to requested dimensions
        idx = randperm (numel (r), ps);
        r = reshape (r(idx), sz);
      else
        r = gprnd (this.k, this.sigma, this.theta, varargin{:});
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {@var{s} =} std (@var{pd})
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
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {@var{t} =} truncate (@var{pd}, @var{lower}, @var{upper})
    ##
    ## Truncate a probability distribution.
    ##
    ## @code{@var{t} = truncate (@var{pd}, @var{lower}, @var{upper})} returns a
    ## probability distribution @var{t}, which is the probability distribution
    ## @var{pd} truncated to the specified interval with lower limit,
    ## @var{lower},
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
      this.Truncation = [lower, upper];
      this.IsTruncated = true;
      this.InputData = [];
      this.ParameterIsFixed = [true, true, true];
      this.ParameterCovariance = zeros (this.NumParameters);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {prob.GeneralizedParetoDistribution} {@var{v} =} var (@var{pd})
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
        [~, v] = gpstat (this.k, this.sigma, this.theta);
      endif
    endfunction

  endmethods

  methods(Static, Hidden)

    function pd = fit (x, theta, varargin)
      ## Check input arguments
      if (nargin < 3)
        alpha = 0.05;
      else
        alpha = varargin{1};
      endif
      if (nargin < 4)
        freq = [];
      else
        freq = varargin{2};
      endif
      if (nargin < 5)
        options.Display = 'off';
        options.MaxFunEvals = 400;
        options.MaxIter = 200;
        options.TolX = 1e-6;
      else
        options = varargin{3};
      endif
      ## Fit data
      [phat, pci] = gpfit (x, theta, alpha, freq, options);
      [~, acov] = gplike (phat, x, freq);
      ## Create fitted distribution object
      pd = prob.GeneralizedParetoDistribution.makeFitted (phat, pci, acov, x, freq);
    endfunction

    function pd = makeFitted (phat, pci, acov, x, freq)
      k = phat(1);
      sigma = phat(2);
      theta = phat(3);
      pd = prob.GeneralizedParetoDistribution (k, sigma, theta);
      pd.ParameterCI = pci;
      pd.ParameterIsFixed = [false, false, true];
      pd.ParameterCovariance = acov;
      pd.InputData = struct ('data', x, 'cens', [], 'freq', freq);
    endfunction

  endmethods

endclassdef

function checkparams (k, sigma, theta)
  if (! (isscalar (k) && isnumeric (k) && isreal (k) && isfinite (k)))
    error ("GeneralizedParetoDistribution: K must be a real scalar.")
  endif
  if (! (isscalar (sigma) && isnumeric (sigma) && isreal (sigma)
                          && isfinite (sigma) && sigma > 0))
    error ("GeneralizedParetoDistribution: SIGMA must be a positive real scalar.")
  endif
  if (! (isscalar (theta) && isnumeric (theta) && isreal (theta)
                          && isfinite (theta)))
    error ("GeneralizedParetoDistribution: THETA must be a real scalar.")
  endif
endfunction

## Test output
%!shared pd, t
%! pd = prob.GeneralizedParetoDistribution (1, 1, 1);
%! t = truncate (pd, 2, 4);
%!assert_equal (cdf (pd, [0:5]), [0, 0, 0.5, 0.6667, 0.75, 0.8], 1e-4);
%!assert_equal (cdf (t, [0:5]), [0, 0, 0, 0.6667, 1, 1], 1e-4);
%!assert_equal (cdf (pd, [1.5, 2, 3, 4]), [0.3333, 0.5, 0.6667, 0.75], 1e-4);
%!assert_equal (cdf (t, [1.5, 2, 3, 4]), [0, 0, 0.6667, 1], 1e-4);
%!assert_equal (icdf (pd, [0:0.2:1]), [1, 1.25, 1.6667, 2.5, 5, Inf], 1e-4);
%!assert_equal (icdf (t, [0:0.2:1]), [2, 2.2222, 2.5, 2.8571, 3.3333, 4], 1e-4);
%!assert_equal (icdf (pd, [-1, 0.4:0.2:1, NaN]), [NaN, 1.6667, 2.5, 5, Inf, NaN], 1e-4);
%!assert_equal (icdf (t, [-1, 0.4:0.2:1, NaN]), [NaN, 2.5, 2.8571, 3.3333, 4, NaN], 1e-4);
%!assert_equal (iqr (pd), 2.6667, 1e-4);
%!assert_equal (iqr (t), 0.9143, 1e-4);
%!assert_equal (mean (pd), Inf);
%!assert_equal (mean (t), 2.7726, 1e-4);
%!assert_equal (median (pd), 2);
%!assert_equal (median (t), 2.6667, 1e-4);
%!assert_equal (pdf (pd, [0:5]), [0, 1, 0.25, 0.1111, 0.0625, 0.04], 1e-4);
%!assert_equal (pdf (t, [0:5]), [0, 0, 1, 0.4444, 0.25, 0], 1e-4);
%!assert_equal (pdf (pd, [-1, 1:4, NaN]), [0, 1, 0.25, 0.1111, 0.0625, NaN], 1e-4);
%!assert_equal (pdf (t, [-1, 1:4, NaN]), [0, 0, 1, 0.4444, 0.25, NaN], 1e-4);
%!assert_equal (isequal (size (random (pd, 100, 50)), [100, 50]), true)
%!assert_equal (any (random (t, 1000, 1) < 2), false);
%!assert_equal (any (random (t, 1000, 1) > 4), false);
%!assert_equal (std (pd), Inf);
%!assert_equal (std (t), 0.5592, 1e-4);
%!assert_equal (var (pd), Inf);
%!assert_equal (var (t), 0.3128, 1e-4);

%!test
%! ## The profile over the first free parameter: 21 grid values, one row of
%! ## OTHER per value, and the likelihood peaking at the fitted estimate.
%! x = [1.2; 0.4; 3.1; 0.7; 2.5; 1.8; 0.3; 4.2; 1.1; 0.9; ...
%!      2.2; 0.6; 1.5; 3.7; 0.8; 2.9; 1.3; 0.5; 2.0; 1.6] - 0.3 + 1e-8;
%! pd = fitdist (x, 'GeneralizedPareto', 'theta', 0);
%! [nlogL, param, other] = proflik (pd, 1);
%! assert_equal (size (param), [1, 21]);
%! assert_equal (size (other), [21, 2]);
%! assert_equal (proflik (pd), nlogL);
%! [~, imax] = max (nlogL);
%! assert_equal (abs (param(imax) - pd.ParameterValues(1)) <= param(2) - param(1), true);

## Test input validation
## 'prob.GeneralizedParetoDistribution' constructor
%!error <GeneralizedParetoDistribution: K must be a real scalar.> ...
%! prob.GeneralizedParetoDistribution (Inf, 1, 1)
%!error <GeneralizedParetoDistribution: K must be a real scalar.> ...
%! prob.GeneralizedParetoDistribution (i, 1, 1)
%!error <GeneralizedParetoDistribution: K must be a real scalar.> ...
%! prob.GeneralizedParetoDistribution ('a', 1, 1)
%!error <GeneralizedParetoDistribution: K must be a real scalar.> ...
%! prob.GeneralizedParetoDistribution ([1, 2], 1, 1)
%!error <GeneralizedParetoDistribution: K must be a real scalar.> ...
%! prob.GeneralizedParetoDistribution (NaN, 1, 1)
%!error <GeneralizedParetoDistribution: SIGMA must be a positive real scalar.> ...
%! prob.GeneralizedParetoDistribution (1, 0, 1)
%!error <GeneralizedParetoDistribution: SIGMA must be a positive real scalar.> ...
%! prob.GeneralizedParetoDistribution (1, -1, 1)
%!error <GeneralizedParetoDistribution: SIGMA must be a positive real scalar.> ...
%! prob.GeneralizedParetoDistribution (1, Inf, 1)
%!error <GeneralizedParetoDistribution: SIGMA must be a positive real scalar.> ...
%! prob.GeneralizedParetoDistribution (1, i, 1)
%!error <GeneralizedParetoDistribution: SIGMA must be a positive real scalar.> ...
%! prob.GeneralizedParetoDistribution (1, 'a', 1)
%!error <GeneralizedParetoDistribution: SIGMA must be a positive real scalar.> ...
%! prob.GeneralizedParetoDistribution (1, [1, 2], 1)
%!error <GeneralizedParetoDistribution: SIGMA must be a positive real scalar.> ...
%! prob.GeneralizedParetoDistribution (1, NaN, 1)
%!error <GeneralizedParetoDistribution: THETA must be a real scalar.> ...
%! prob.GeneralizedParetoDistribution (1, 1, Inf)
%!error <GeneralizedParetoDistribution: THETA must be a real scalar.> ...
%! prob.GeneralizedParetoDistribution (1, 1, i)
%!error <GeneralizedParetoDistribution: THETA must be a real scalar.> ...
%! prob.GeneralizedParetoDistribution (1, 1, 'a')
%!error <GeneralizedParetoDistribution: THETA must be a real scalar.> ...
%! prob.GeneralizedParetoDistribution (1, 1, [1, 2])
%!error <GeneralizedParetoDistribution: THETA must be a real scalar.> ...
%! prob.GeneralizedParetoDistribution (1, 1, NaN)

## 'cdf' method
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (prob.GeneralizedParetoDistribution, 2, 'uper')
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (prob.GeneralizedParetoDistribution, 2, 3)

## 'paramci' method
%!shared x
%! x = gprnd (1, 1, 1, [1, 100]);
%!error <paramci: optional arguments must be in NAME-VALUE pairs.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), 'alpha')
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), 'alpha', 0)
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), 'alpha', 1)
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), 'alpha', [0.5 2])
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), 'alpha', '')
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), 'alpha', {0.05})
%!error <paramci: invalid VALUE for 'Alpha' argument.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), ...
%!          'parameter', 'sigma', 'alpha', {0.05})
%!error <paramci: invalid VALUE size for 'Parameter' argument.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), ...
%!          'parameter', {'k', 'sigma', 'theta', 'param'})
%!error <paramci: invalid VALUE size for 'Parameter' argument.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), 'alpha', 0.01, ...
%!          'parameter', {'k', 'sigma', 'theta', 'param'})
%!error <paramci: unknown distribution parameter.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), 'parameter', 'param')
%!error <paramci: unknown distribution parameter.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), 'alpha', 0.01, ...
%!          'parameter', 'param')
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), 'NAME', 'value')
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), 'alpha', 0.01, ...
%!          'NAME', 'value')
%!error <paramci: invalid NAME for optional argument.> ...
%! paramci (prob.GeneralizedParetoDistribution.fit (x, 1), 'alpha', 0.01, ...
%!          'parameter', 'sigma', 'NAME', 'value')

## 'plot' method
%!error <plot: optional arguments must be in NAME-VALUE pairs.> ...
%! plot (prob.GeneralizedParetoDistribution, 'Parent')
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (prob.GeneralizedParetoDistribution, 'PlotType', 12)
%!error <plot: invalid VALUE size for 'Parameter' argument.> ...
%! plot (prob.GeneralizedParetoDistribution, 'PlotType', {'pdf', 'cdf'})
%!error <plot: invalid VALUE for 'PlotType' argument.> ...
%! plot (prob.GeneralizedParetoDistribution, 'PlotType', 'pdfcdf')
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (prob.GeneralizedParetoDistribution, 'Discrete', 'pdfcdf')
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (prob.GeneralizedParetoDistribution, 'Discrete', [1, 0])
%!error <plot: invalid VALUE for 'Discrete' argument.> ...
%! plot (prob.GeneralizedParetoDistribution, 'Discrete', {true})
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (prob.GeneralizedParetoDistribution, 'Parent', 12)
%!error <plot: invalid VALUE for 'Parent' argument.> ...
%! plot (prob.GeneralizedParetoDistribution, 'Parent', 'hax')
%!error <plot: invalid NAME for optional argument.> ...
%! plot (prob.GeneralizedParetoDistribution, 'invalidNAME', 'pdf')
%!error <plot: no fitted DATA to plot a probability plot.> ...
%! plot (prob.GeneralizedParetoDistribution, 'PlotType', 'probability')

## 'proflik' method
%!error <proflik: no fitted data available.> ...
%! proflik (prob.GeneralizedParetoDistribution, 2)
%!error <proflik: PNUM must be a scalar number indexing a non-fixed parameter.> ...
%! proflik (prob.GeneralizedParetoDistribution.fit (x, 1), 3)
%!error <proflik: PNUM must be a scalar number indexing a non-fixed parameter.> ...
%! proflik (prob.GeneralizedParetoDistribution.fit (x, 1), [1, 2])
%!error <proflik: PNUM must be a scalar number indexing a non-fixed parameter.> ...
%! proflik (prob.GeneralizedParetoDistribution.fit (x, 1), {1})
%!error <proflik: SETPARAM must be a numeric vector.> ...
%! proflik (prob.GeneralizedParetoDistribution.fit (x, 1), 1, ones (2))
%!error <proflik: missing VALUE for 'Display' argument.> ...
%! proflik (prob.GeneralizedParetoDistribution.fit (x, 1), 1, 'Display')
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (prob.GeneralizedParetoDistribution.fit (x, 1), 1, 'Display', 1)
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (prob.GeneralizedParetoDistribution.fit (x, 1), 1, 'Display', {1})
%!error <proflik: invalid VALUE type for 'Display' argument.> ...
%! proflik (prob.GeneralizedParetoDistribution.fit (x, 1), 1, 'Display', {'on'})
%!error <proflik: invalid VALUE size for 'Display' argument.> ...
%! proflik (prob.GeneralizedParetoDistribution.fit (x, 1), 1, ...
%!          'Display', ['on'; 'on'])
%!error <proflik: invalid VALUE for 'Display' argument.> ...
%! proflik (prob.GeneralizedParetoDistribution.fit (x, 1), 1, 'Display', 'onnn')
%!error <proflik: invalid NAME for optional arguments.> ...
%! proflik (prob.GeneralizedParetoDistribution.fit (x, 1), 1, 'NAME', 'on')
%!error <proflik: invalid optional argument.> ...
%! proflik (prob.GeneralizedParetoDistribution.fit (x, 1), 1, {'NAME'}, 'on')
%!error <proflik: invalid optional argument.> ...
%! proflik (prob.GeneralizedParetoDistribution.fit (x, 1), 1, {[1 2 3 4]}, ...
%!          'Display', 'on')

## 'truncate' method
%!error <truncate: missing input argument.> ...
%! truncate (prob.GeneralizedParetoDistribution)
%!error <truncate: missing input argument.> ...
%! truncate (prob.GeneralizedParetoDistribution, 2)
%!error <truncate: invalid lower upper limits.> ...
%! truncate (prob.GeneralizedParetoDistribution, 4, 2)

## Catch errors when using array of probability objects with available methods
%!shared pd
%! pd = prob.GeneralizedParetoDistribution (1, 1, 1);
%! pd(2) = prob.GeneralizedParetoDistribution (1, 3, 1);
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
