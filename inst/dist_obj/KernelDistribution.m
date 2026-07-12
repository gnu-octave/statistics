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

classdef KernelDistribution
  ## -*- texinfo -*-
  ## @deftp {statistics} KernelDistribution
  ##
  ## Kernel probability distribution object.
  ##
  ## A @code{KernelDistribution} object consists of a nonparametric kernel
  ## smoothing density estimate fitted to sample data, together with a model
  ## description.  Unlike the parametric distribution objects, it has no
  ## estimated parameters; the fitted distribution is defined entirely by the
  ## data, the smoothing kernel, and the bandwidth.
  ##
  ## A @code{KernelDistribution} object can only be created by fitting a kernel
  ## smoothing distribution to data with the @code{fitdist} function.  Unlike
  ## the parametric distributions, it cannot be created with the @code{makedist}
  ## function, since it is not parametric and requires data.
  ##
  ## Further information about the kernel density estimation can be found at
  ## @url{https://en.wikipedia.org/wiki/Kernel_density_estimation}
  ##
  ## @seealso{fitdist, ksdensity, mvksdensity}
  ## @end deftp

  properties (Dependent = true)
    ## -*- texinfo -*-
    ## @deftp {KernelDistribution} {property} Kernel
    ##
    ## Kernel smoothing function
    ##
    ## A character vector specifying the type of smoothing kernel used for the
    ## density estimate.  It is one of @qcode{'normal'}, @qcode{'box'},
    ## @qcode{'triangle'}, or @qcode{'epanechnikov'}.  You can access the
    ## @qcode{Kernel} property using dot name assignment.
    ##
    ## @end deftp
    Kernel

    ## -*- texinfo -*-
    ## @deftp {KernelDistribution} {property} Bandwidth
    ##
    ## Bandwidth of the smoothing kernel
    ##
    ## A positive scalar value specifying the bandwidth of the smoothing kernel.
    ## You can access the @qcode{Bandwidth} property using dot name assignment.
    ##
    ## @end deftp
    Bandwidth
  endproperties

  properties (GetAccess = public, Constant = true)
    ## -*- texinfo -*-
    ## @deftp {KernelDistribution} {property} DistributionName
    ##
    ## Probability distribution name
    ##
    ## A character vector specifying the name of the probability distribution
    ## object.  This property is read-only.
    ##
    ## @end deftp
    DistributionName = 'KernelDistribution';

    ## -*- texinfo -*-
    ## @deftp {KernelDistribution} {property} NumParameters
    ##
    ## Number of parameters
    ##
    ## A scalar integer value specifying the number of parameters characterizing
    ## the probability distribution.  A kernel distribution is nonparametric, so
    ## this value is always @qcode{0}.  This property is read-only.
    ##
    ## @end deftp
    NumParameters = 0;

    ## -*- texinfo -*-
    ## @deftp {KernelDistribution} {property} ParameterNames
    ##
    ## Names of parameters
    ##
    ## An empty cell array, since a kernel distribution has no estimated
    ## parameters.  This property is read-only.
    ##
    ## @end deftp
    ParameterNames = {};

    ## -*- texinfo -*-
    ## @deftp {KernelDistribution} {property} ParameterDescription
    ##
    ## Description of parameters
    ##
    ## An empty cell array, since a kernel distribution has no estimated
    ## parameters.  This property is read-only.
    ##
    ## @end deftp
    ParameterDescription = {};
  endproperties

  properties (GetAccess = public, Constant = true, Hidden)
    CensoringAllowed = false;
    DistributionCode = 'kernel';
  endproperties

  properties (GetAccess = public, SetAccess = protected)
    ## -*- texinfo -*-
    ## @deftp {KernelDistribution} {property} Support
    ##
    ## Support of the probability distribution
    ##
    ## A scalar structure containing the following fields:
    ## @itemize
    ## @item @qcode{range}: either the character vector @qcode{'unbounded'} or
    ## @qcode{'positive'}, or a two-element numeric vector @math{[L, U]} with the
    ## lower and upper bounds of the support.
    ## @item @qcode{closedbound}: a two-element logical vector specifying whether
    ## each bound is closed.
    ## @item @qcode{iscontinuous}: a logical scalar, always @qcode{true} for a
    ## kernel distribution.
    ## @end itemize
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Support

    ## -*- texinfo -*-
    ## @deftp {KernelDistribution} {property} Truncation
    ##
    ## Truncation interval
    ##
    ## A @math{1*2} numeric vector specifying the truncation interval for the
    ## probability distribution.  First element contains the lower boundary,
    ## second element contains the upper boundary.  This property is read-only.
    ## You can only truncate a probability distribution with the
    ## @qcode{truncate} method.
    ##
    ## @end deftp
    Truncation

    ## -*- texinfo -*-
    ## @deftp {KernelDistribution} {property} IsTruncated
    ##
    ## Flag for truncated probability distribution
    ##
    ## A logical scalar value specifying whether a probability distribution is
    ## truncated or not.  This property is read-only.
    ##
    ## @end deftp
    IsTruncated

    ## -*- texinfo -*-
    ## @deftp {KernelDistribution} {property} InputData
    ##
    ## Data used for fitting a probability distribution
    ##
    ## A scalar structure containing the following fields:
    ## @itemize
    ## @item @qcode{data}: a numeric vector containing the data used for
    ## distribution fitting.
    ## @item @qcode{cens}: an empty array, since censoring is not supported for
    ## a kernel distribution.
    ## @item @qcode{freq}: a numeric vector of non-negative integer values
    ## containing the frequency information corresponding to the elements of the
    ## data used for distribution fitting.  If no frequency vector was used for
    ## distribution fitting, then this field defaults to an empty array.
    ## @end itemize
    ##
    ## @end deftp
    InputData
  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    KernelName
    BandwidthValue
  endproperties

  methods (Hidden)

    function this = KernelDistribution (data, kernel, bw, support, freq)
      if (nargin == 0)
        data = [0; 1];
        kernel = 'normal';
        support = make_support ('unbounded');
        freq = [];
        bw = default_bandwidth (data, kernel, support, freq);
      endif
      this.KernelName = kernel;
      this.BandwidthValue = bw;
      this.Support = support;
      this.IsTruncated = false;
      this.Truncation = [];
      this.InputData = struct ('data', data(:), 'cens', [], 'freq', freq);
    endfunction

    function display (this)
      fprintf ("%s =\n", inputname (1));
      __disp__ (this, 'kernel distribution');
    endfunction

    function disp (this)
      __disp__ (this, 'kernel distribution');
    endfunction

    function this = set.Kernel (this, kernel)
      if (! (ischar (kernel) && isrow (kernel) && any (strcmpi (kernel, ...
             {'normal', 'box', 'triangle', 'epanechnikov'}))))
        error (strcat ("KernelDistribution: 'Kernel' must be 'normal',", ...
                       " 'box', 'triangle', or 'epanechnikov'."));
      endif
      this.KernelName = lower (kernel);
    endfunction

    function kernel = get.Kernel (this)
      kernel = this.KernelName;
    endfunction

    function this = set.Bandwidth (this, bw)
      if (! (isnumeric (bw) && isscalar (bw) && isreal (bw) && bw > 0))
        error ("KernelDistribution: 'Bandwidth' must be a positive scalar.");
      endif
      this.BandwidthValue = bw;
    endfunction

    function bw = get.Bandwidth (this)
      bw = this.BandwidthValue;
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {KernelDistribution} {@var{p} =} cdf (@var{pd}, @var{x})
    ## @deftypefnx {KernelDistribution} {@var{p} =} cdf (@var{pd}, @var{x}, @qcode{'upper'})
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
      p = basecdf (this, x);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        p(lb) = 0;
        p(ub) = 1;
        Fab = basecdf (this, [lx, ux]);
        p(! (lb | ub)) -= Fab(1);
        p(! (lb | ub)) /= diff (Fab);
      endif
      ## Apply uflag
      if (utail)
        p = 1 - p;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {KernelDistribution} {@var{x} =} icdf (@var{pd}, @var{p})
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
        Fab = basecdf (this, this.Truncation);
        lp = Fab(1);
        up = Fab(2);
        ## Adjust p values within range of p @ lower limit and p @ upper limit
        is_nan = p < 0 | p > 1;
        p(is_nan) = NaN;
        np = lp + (up - lp) .* p;
        x = baseicdf (this, np);
        x(x < this.Truncation(1)) = this.Truncation(1);
        x(x > this.Truncation(2)) = this.Truncation(2);
      else
        x = baseicdf (this, p);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {KernelDistribution} {@var{r} =} iqr (@var{pd})
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
    ## @deftypefn  {KernelDistribution} {@var{m} =} mean (@var{pd})
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
      if (this.IsTruncated || ! is_unbounded (this))
        [lb, ub] = integration_limits (this);
        fm = @(x) x .* pdf (this, x);
        m = integral (fm, lb, ub, 'ArrayValued', 1);
      else
        w = weights (this);
        m = sum (w .* this.InputData.data(:)');
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {KernelDistribution} {@var{m} =} median (@var{pd})
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
    ## @deftypefn  {KernelDistribution} {@var{nlogL} =} negloglik (@var{pd})
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
      data = this.InputData.data(:);
      f = pdf (this, data);
      if (isempty (this.InputData.freq))
        nlogL = - sum (log (f));
      else
        nlogL = - sum (this.InputData.freq(:) .* log (f));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {KernelDistribution} {@var{y} =} pdf (@var{pd}, @var{x})
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
      y = basepdf (this, x);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        y(lb | ub) = 0;
        y(! (lb | ub)) /= diff (basecdf (this, [lx, ux]));
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {KernelDistribution} {} plot (@var{pd})
    ## @deftypefnx {KernelDistribution} {} plot (@var{pd}, @var{Name}, @var{Value})
    ## @deftypefnx {KernelDistribution} {@var{h} =} plot (@dots{})
    ##
    ## Plot a probability distribution object.
    ##
    ## @code{plot (@var{pd})} plots a probability density function (PDF) of the
    ## probability distribution object @var{pd}, superimposed over a histogram
    ## of the data used to fit it.
    ##
    ## @code{plot (@var{pd}, @var{Name}, @var{Value})} specifies additional
    ## options with the @qcode{Name-Value} pair arguments listed below.
    ##
    ## @multitable @columnfractions 0.18 0.82
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'PlotType'} @tab A character vector specifying the plot
    ## type.  @qcode{'pdf'} plots the probability density function (PDF)
    ## superimposed on a histogram of the data.  @qcode{'cdf'} plots the
    ## cumulative distribution function (CDF) superimposed over an empirical
    ## CDF.  @qcode{'probability'} plots a probability plot using a CDF of the
    ## data and a CDF of the fitted probability distribution.
    ##
    ## @item @qcode{'Parent'} @tab An axes graphics object for plot.  If not
    ## specified, the @code{plot} function plots into the current axes or
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
    ## @deftypefn  {KernelDistribution} {@var{r} =} random (@var{pd})
    ## @deftypefnx {KernelDistribution} {@var{r} =} random (@var{pd}, @var{rows})
    ## @deftypefnx {KernelDistribution} {@var{r} =} random (@var{pd}, @var{rows}, @var{cols}, @dots{})
    ## @deftypefnx {KernelDistribution} {@var{r} =} random (@var{pd}, [@var{sz}])
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
      u = unifrnd (0, 1, varargin{:});
      r = icdf (this, u);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {KernelDistribution} {@var{s} =} std (@var{pd})
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
      s = sqrt (var (this));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {KernelDistribution} {@var{t} =} truncate (@var{pd}, @var{lower}, @var{upper})
    ##
    ## Truncate a probability distribution.
    ##
    ## @code{@var{t} = truncate (@var{pd}, @var{lower}, @var{upper})} returns a
    ## probability distribution @var{t}, which is the probability distribution
    ## @var{pd} truncated to the specified interval with lower limit,
    ## @var{lower}, and upper limit, @var{upper}.
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
    ## @deftypefn  {KernelDistribution} {@var{v} =} var (@var{pd})
    ##
    ## Compute the variance of a probability distribution.
    ##
    ## @code{@var{v} = var (@var{pd})} computes the variance of the probability
    ## distribution object, @var{pd}.
    ##
    ## @end deftypefn
    function v = var (this)
      if (! isscalar (this))
        error ("var: requires a scalar probability distribution.");
      endif
      if (this.IsTruncated || ! is_unbounded (this))
        [lb, ub] = integration_limits (this);
        m = mean (this);
        fv = @(x) ((x - m) .^ 2) .* pdf (this, x);
        v = integral (fv, lb, ub, 'ArrayValued', 1);
      else
        w = weights (this);
        data = this.InputData.data(:)';
        mu = sum (w .* data);
        v = sum (w .* (data - mu) .^ 2) + this.Bandwidth ^ 2;
      endif
    endfunction

  endmethods

  methods (Access = private)

    ## Untruncated PDF, delegating to the shipped ksdensity engine.  The query
    ## values are flattened to a vector for ksdensity and reshaped back, since
    ## ksdensity only accepts vector query points.
    function y = basepdf (this, x)
      y = reshape (ksdensity (this.InputData.data, x(:), ksargs (this){:}, ...
                              'Function', 'pdf'), size (x));
    endfunction

    ## Untruncated CDF, delegating to the shipped ksdensity engine.
    function p = basecdf (this, x)
      p = reshape (ksdensity (this.InputData.data, x(:), ksargs (this){:}, ...
                              'Function', 'cdf'), size (x));
    endfunction

    ## Untruncated inverse CDF, delegating to the shipped ksdensity engine.
    function x = baseicdf (this, p)
      x = reshape (ksdensity (this.InputData.data, p(:), ksargs (this){:}, ...
                              'Function', 'icdf'), size (p));
    endfunction

    ## Common Name-Value arguments forwarded to ksdensity.
    function args = ksargs (this)
      args = {'Kernel', this.Kernel, 'Bandwidth', this.Bandwidth, ...
              'Support', support_arg(this.Support)};
      if (! isempty (this.InputData.freq))
        args = [args, {'Weights', this.InputData.freq}];
      endif
    endfunction

    ## Normalised weights of the sample points.
    function w = weights (this)
      n = numel (this.InputData.data);
      if (isempty (this.InputData.freq))
        w = ones (1, n) / n;
      else
        f = this.InputData.freq(:)';
        w = f / sum (f);
      endif
    endfunction

    ## Finite integration limits for numeric moment computation.
    function [lb, ub] = integration_limits (this)
      if (this.IsTruncated)
        lb = this.Truncation(1);
        ub = this.Truncation(2);
      else
        lb = baseicdf (this, 1e-10);
        ub = baseicdf (this, 1 - 1e-10);
      endif
    endfunction

    ## True for an unbounded (whole real line) support.
    function tf = is_unbounded (this)
      tf = ischar (this.Support.range) ...
           && strcmp (this.Support.range, 'unbounded');
    endfunction

  endmethods

  methods (Static, Hidden)

    function pd = fit (x, kernel, support, width, freq)
      if (nargin < 2 || isempty (kernel))
        kernel = 'normal';
      endif
      if (nargin < 3 || isempty (support))
        support = 'unbounded';
      endif
      if (nargin < 4)
        width = [];
      endif
      if (nargin < 5)
        freq = [];
      endif
      ## A trivial (all-ones) frequency vector is treated as no weighting
      if (! isempty (freq) && all (freq(:) == 1))
        freq = [];
      endif
      ## Validate the kernel
      if (! (ischar (kernel) && isrow (kernel) && any (strcmpi (kernel, ...
             {'normal', 'box', 'triangle', 'epanechnikov'}))))
        error (strcat ("KernelDistribution: 'Kernel' must be 'normal',", ...
                       " 'box', 'triangle', or 'epanechnikov'."));
      endif
      kernel = lower (kernel);
      S = make_support (support);
      ## Bandwidth: use the supplied value or the ksdensity default rule
      if (isempty (width))
        bw = default_bandwidth (x, kernel, S, freq);
      elseif (isnumeric (width) && isscalar (width) && isreal (width)
                                                    && width > 0)
        bw = width;
      else
        error ("KernelDistribution: 'Width' must be a positive scalar.");
      endif
      ## Create fitted distribution object
      pd = KernelDistribution (x, kernel, bw, S, freq);
    endfunction

  endmethods

endclassdef

## Build the Support structure from a keyword or a two-element vector.
function S = make_support (support)
  if (ischar (support) && isrow (support))
    switch (lower (support))
      case 'unbounded'
        S = struct ('range', 'unbounded', 'closedbound', ...
                    [false, false], 'iscontinuous', true);
      case 'positive'
        S = struct ('range', 'positive', 'closedbound', ...
                    [false, false], 'iscontinuous', true);
      otherwise
        error (strcat ("KernelDistribution: 'Support' must be", ...
                       " 'unbounded', 'positive', or [L U]."));
    endswitch
  elseif (isnumeric (support) && isreal (support) && numel (support) == 2
                                          && support(1) < support(2))
    S = struct ('range', [support(1), support(2)], 'closedbound', ...
                [false, false], 'iscontinuous', true);
  else
    error (strcat ("KernelDistribution: 'Support' must be", ...
                   " 'unbounded', 'positive', or [L U]."));
  endif
endfunction

## Convert a Support structure to the argument expected by ksdensity.
function arg = support_arg (S)
  if (ischar (S.range))
    arg = S.range;
  else
    arg = S.range;
  endif
endfunction

## Default bandwidth via the ksdensity normal-reference rule.
function bw = default_bandwidth (x, kernel, S, freq)
  args = {'Kernel', kernel, 'Support', support_arg(S)};
  if (! isempty (freq))
    args = [args, {'Weights', freq}];
  endif
  [~, ~, bw] = ksdensity (x, args{:});
endfunction

%!demo
%! ## Fit a kernel distribution to a sample and plot its PDF over a histogram.
%! load patients
%! pd = fitdist (Weight, 'Kernel');
%! plot (pd)
%! title ('Kernel distribution fitted to patient weights')

## Test output
%!shared x, pd, pdbox, t
%! x = [2.1 0.3 1.2 -0.7 0.9 1.5 2.8 0.1 0.4 1.1 3.2 0.6 2.0 0.9 1.7]';
%! pd = fitdist (x, 'Kernel');
%! pdbox = fitdist (x, 'Kernel', 'Kernel', 'box', 'Width', 0.5);
%! t = truncate (pd, 0, 3);
%!assert_equal (pd.DistributionName, 'KernelDistribution');
%!assert_equal (pd.Kernel, 'normal');
%!assert_equal (pd.NumParameters, 0);
%!assert_equal (pd.ParameterNames, {});
%!assert_equal (pd.IsTruncated, false);
%!assert_equal (pd.Bandwidth, 0.639566, 1e-4);
%!assert_equal (pd.Support.range, 'unbounded');
%!assert_equal (pdf (pd, [-1 0 0.5 1 1.5 2 2.5 3]), ...
%!    [0.0589 0.2141 0.3051 0.3394 0.3061 0.2375 0.1697 0.1166], 2e-3);
%!assert_equal (cdf (pd, [-1 0 0.5 1 1.5 2 2.5 3]), ...
%!    [0.0272 0.1538 0.2850 0.4492 0.6128 0.7494 0.8506 0.9217], 2e-3);
%!assert_equal (icdf (pd, [0.1 0.25 0.5 0.75 0.9]), ...
%!    [-0.2909 0.3820 1.1504 2.0027 2.8265], 5e-3);
%!assert_equal (mean (pd), 1.2067, 1e-4);
%!assert_equal (std (pd), 1.1918, 1e-3);
%!assert_equal (var (pd), 1.4203, 1e-3);
%!assert_equal (median (pd), 1.1504, 5e-3);
%!assert_equal (iqr (pd), 1.6207, 5e-3);
%!assert_equal (negloglik (pd), 21.5835, 1e-3);
%!assert_equal (pdbox.Kernel, 'box');
%!assert_equal (pdbox.Bandwidth, 0.5);
%!assert_equal (pdf (pdbox, [0 1 2]), [0.1925 0.3464 0.2309], 2e-3);
%!assert_equal (pdf (t, [-1 0 1 2 3 4]), ...
%!    [0 0.2788 0.4420 0.3093 0.1518 0], 2e-3);
%!assert_equal (cdf (t, [-1 0 1 2 3 4]), ...
%!    [0 0 0.3846 0.7755 1 1], 2e-3);
%!assert_equal (mean (t), 1.3293, 1e-3);

%!test  ## positive support (log boundary correction)
%! y = [0.2 0.5 0.7 1.1 1.4 2 2.6 3.3 4.1 5.5]';
%! pdpos = fitdist (y, 'Kernel', 'Support', 'positive');
%! assert_equal (pdpos.Bandwidth, 0.768201, 1e-4);
%! assert_equal (pdpos.Support.range, 'positive');
%! assert_equal (pdf (pdpos, [0.1 0.5 1 2 3 5]), ...
%!    [0.4302 0.3919 0.2738 0.1572 0.0999 0.0463], 2e-3);

%!test  ## random values honour truncation bounds
%! x = [2.1 0.3 1.2 -0.7 0.9 1.5 2.8 0.1 0.4 1.1 3.2 0.6 2.0 0.9 1.7]';
%! pd = fitdist (x, 'Kernel');
%! t = truncate (pd, 0, 3);
%! r = random (t, 1000, 1);
%! assert_equal (all (r >= 0 & r <= 3), true);
%! assert_equal (size (random (pd, 10, 5)), [10, 5]);

## Test input validation
%!error <cdf: invalid argument for upper tail.> ...
%! cdf (fitdist ([1 2 3 4 5]', 'Kernel'), 2, 'uper')
%!error <truncate: missing input argument.> ...
%! truncate (fitdist ([1 2 3 4 5]', 'Kernel'))
%!error <truncate: invalid lower upper limits.> ...
%! truncate (fitdist ([1 2 3 4 5]', 'Kernel'), 4, 2)
%!error <KernelDistribution: 'Kernel' must be 'normal', 'box', 'triangle', or 'epanechnikov'.> ...
%! fitdist ([1 2 3 4 5]', 'Kernel', 'Kernel', 'cosine')
%!error <KernelDistribution: 'Support' must be 'unbounded', 'positive', or .L U..> ...
%! fitdist ([1 2 3 4 5]', 'Kernel', 'Support', 'half')
%!error <KernelDistribution: 'Width' must be a positive scalar.> ...
%! fitdist ([1 2 3 4 5]', 'Kernel', 'Width', -1)

## Catch errors when using array of probability objects with available methods
%!shared pd
%! pd = KernelDistribution ();
%! pd(2) = KernelDistribution ();
%!error <cdf: requires a scalar probability distribution.> cdf (pd, 1)
%!error <icdf: requires a scalar probability distribution.> icdf (pd, 0.5)
%!error <iqr: requires a scalar probability distribution.> iqr (pd)
%!error <mean: requires a scalar probability distribution.> mean (pd)
%!error <median: requires a scalar probability distribution.> median (pd)
%!error <negloglik: requires a scalar probability distribution.> negloglik (pd)
%!error <pdf: requires a scalar probability distribution.> pdf (pd, 1)
%!error <plot: requires a scalar probability distribution.> plot (pd)
%!error <random: requires a scalar probability distribution.> random (pd)
%!error <std: requires a scalar probability distribution.> std (pd)
%!error <truncate: requires a scalar probability distribution.> ...
%! truncate (pd, 2, 4)
%!error <var: requires a scalar probability distribution.> var (pd)
