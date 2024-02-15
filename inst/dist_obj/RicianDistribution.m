classdef RicianDistribution

  properties (Dependent = true)
    nu
    sigma
  endproperties

  properties (GetAccess = public, Constant = true)
    DistributionName = "RicianDistribution";
    DistributionCode = "rice";
    CensoringAllowed = true;
    NumParameters = 2;
    ParameterNames = {"nu", "sigma"};
    ParameterRange = [0, realmin; Inf, Inf];
    ParameterLogCI = [true, true];
    ParameterDescription = {"Non-centrality Distance", "Scale"};
  endproperties

  properties (GetAccess = public , SetAccess = protected)
    ParameterValues
    ParameterCI
    NegativeLogLikelihood
    ParameterCovariance
    ParameterIsFixed
    Truncation
    IsTruncated
    InputData
  endproperties

  methods (Hidden)

    function pd = tLocationScaleDistribution (nu, sigma)
      if (nargin == 0)
        nu = 0;
        sigma = 1;
      endif
      checkparams (nu, sigma)
      pd.ParameterValues = [nu, sigma];
      pd.ParameterIsFixed = [false, false];
      pd.ParameterCovariance = zeros (pd.NumParameters);
      pd.IsTruncated = false;
    endfunction

    function display (this)
      fprintf ("%s =\n", inputname(1));
      __disp__ (this, "Rician distribution");
    endfunction

    function disp (this)
      __disp__ (this, "Rician distribution");
    endfunction

    function this = set.nu (this, nu)
        checkparams (nu, this.sigma)
        this.ParameterValues(1) = nu;
        this.ParameterCovariance = zeros (this.NumParameters);
        this.InputData = [];
    endfunction

    function nu = get.nu (this)
        nu = this.ParameterValues(1);
    endfunction

    function this = set.sigma (this, sigma)
        checkparams (this.nu, sigma)
        this.ParameterValues(2) = sigma;
        this.ParameterCovariance = zeros (this.NumParameters);
        this.InputData = [];
    endfunction

    function sigma = get.sigma (this)
        sigma = this.ParameterValues(2);
    endfunction

  endmethods

  methods (Access = public)

    function p = cdf (this, x, uflag)
      p = ricecdf (x, this.nu, this.sigma, uflag);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        if (strcmpi (uflag, "upper"))
          p(lb) = 1;
          p(ub) = 0;
          p(! (lb | ub)) += ricecdf (lx, this.nu, this.sigma, uflag);
          p(! (lb | ub)) /= diff (ricecdf ([lx, ux], this.nu, this.sigma, uflag));
        else
          p(lb) = 0;
          p(ub) = 1;
          p(! (lb | ub)) -= ricecdf (lx, this.nu, this.sigma);
          p(! (lb | ub)) /= diff (ricecdf ([ux, lx], this.nu, this.sigma));
        endif
      endif
    endfunction

    function x = icdf (this, p)
      x = riceinv (p, this.nu, this.sigma);
      if (this.IsTruncated)
        x(x < this.Truncation(1)) = this.Truncation(1);
        x(x > this.Truncation(2)) = this.Truncation(2);
      endif
    endfunction

    function r = irq (this)
      if (this.IsTruncated)
        lp = ricecdf (this.Truncation(1), this.nu, this.sigma);
        up = ricecdf (this.Truncation(2), this.nu, this.sigma);
        r = diff (riceinv ([0.75-up, 0.25+lb], this.nu, this.sigma));
      else
        r = diff (riceinv ([0.75, 0.25], this.nu, this.sigma));
      endif
    endfunction

    function m = mean (this)
      m = ricestat (this.nu, this.sigma);
    endfunction

    function m = median (this)
      m = ricestat (this.nu, this.sigma);
    endfunction

    function nlogL = negloglik (this)
      if (isempty (this.InputData))
        error ("negloglik: no data available");
      endif
      params = [this.nu, this.sigma];
      nlogL = ricelike (params, this.InputData.x, ...
                       this.InputData.cens, this.InputData.freq);
    endfunction

    function ci = paramci (this, varargin)
      if (isempty (this.InputData))
        ci = [this.ParameterValues; this.ParameterValues];
      else
        ci = __paramci__ (this, varargin);
      endif
    endfunction

    function y = pdf (this, x)
      y = ricepdf (x, this.nu, this.sigma);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        y(lb | ub) = 0;
        y(! (lb | ub)) /= diff (ricepdf ([ux, lx], this.nu, this.sigma));
      endif
    endfunction

    function [varargout] = plot (this, varargin)
      h = __plot__ (this, false, varargin{:});
      if (nargout > 0)
        varargout{1} = h;
      endif
    endfunction

    function [varargout] = proflik (this, pnum, varargin)
      [varargout{1:nargout}] = __proflik__ (this, pnum, varargin);
    endfunction

    function r = random (this, varargin)
      if (this.IsTruncated)
        sz = [varargin{:}];
        ps = prod (sz);
        ## Get an estimate of how many more random numbers we need to randomly
        ## pick the appropriate size from
        lx = this.Truncation(1);
        ux = this.Truncation(2);
        ratio = 1 / diff (ricecdf ([ux, lx], this.nu, this.sigma));
        nsize = 2 * ratio * ps;       # times 2 to be on the safe side
        ## Generate the numbers and remove out-of-bound random samples
        r = ricernd (this.nu, this.sigma, nsize, 1);
        r(r < lx | r > ux) = [];
        ## Randomly select the required size and reshape to requested dimensions
        r = randperm (r, ps);
        r = reshape (r, sz);
      else
        r = ricernd (this.nu, this.sigma, varargin{:});
      endif
    endfunction

    function s = std (this)
      [~, v] = ricestat (this.nu, this.sigma);
      s = sqrt (v);
    endfunction

    function t = truncate (this, lower, upper)
      if (nargin < 3)
        error ("truncate: missing input argument.");
      elseif (lower >= upper)
        error ("truncate: invalid lower upper limits.");
      endif
      this.Truncation = [lower, upper];
      this.IsTruncated = true;
      this.InputData = [];
      t = this;
    endfunction

    function v = var (this)
      [~, v] = ricestat (this.nu, this.sigma);
    endfunction

  endmethods

  methods (Static, Hidden)

    function pd = fit(varargin)
      [paramhat, paramci] = ricefit (x, 0.05, censor, freq, opt);
      [nlogL, acov] = ricelike (paramhat, x, censor, freq);
      pd = RicianDistribution.makeFitted ...
           (paramhat, paramci, nlogL, acov, x, censor, freq);
    endfunction

    function pd = makeFitted (paramhat, paramci, nlogL, acov, x, censor, freq)
      nu = paramhat(1);
      sigma = paramhat(2);
      pd = RicianDistribution (nu, sigma);
      pd.ParameterCI = paramci;
      pd.NegativeLogLikelihood = nlogL;
      pd.ParameterCovariance = acov;
      pd.InputData = struct ("data", x, "cens", censor, "freq", freq);
    endfunction

  endmethods

endclassdef

function checkparams (nu, sigma)
  if (! (isscalar (nu) && isnumeric (nu) && isreal (nu) && isfinite (nu)
                       && nu >= 0 ))
    error ("tLocationScaleDistribution: NU must be a non-negative real scalar.")
  endif
  if (! (isscalar (sigma) && isnumeric (sigma) && isreal (sigma)
                          && isfinite (sigma) && sigma > 0))
    error ("tLocationScaleDistribution: SIGMA must be a positive real scalar.")
  endif
endfunction
