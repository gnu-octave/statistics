classdef tLocationScaleDistribution

  properties (Dependent = true)
    mu
    sigma
    df
  endproperties

  properties (GetAccess = public, Constant = true)
    DistributionName = "tLocationScaleDistribution";
    DistributionCode = "tls";
    CensoringAllowed = true;
    NumParameters = 3;
    ParameterNames = {"mu", "sigma", "df"};
    ParameterRange = [-Inf, realmin, realmin; Inf, Inf, Inf];
    ParameterLogCI = [true, true, true];
    ParameterDescription = {"Location", "Scale", "Degrees of Freedom"};
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

    function pd = tLocationScaleDistribution (mu, sigma, df)
      if (nargin == 0)
        mu = 0;
        sigma = 1;
        df = 5;
      endif
      checkparams (mu, sigma, df)
      pd.ParameterValues = [mu, sigma, df];
      pd.ParameterIsFixed = [false, false, false];
      pd.ParameterCovariance = zeros (pd.NumParameters);
      pd.IsTruncated = false;
    endfunction

    function display (this)
      fprintf ("%s =\n", inputname(1));
      __disp__ (this, "t Location-Scale distribution");
    endfunction

    function disp (this)
      __disp__ (this, "t Location-Scale distribution");
    endfunction

    function this = set.mu (this, mu)
        checkparams (mu, this.sigma, this.df)
        this.ParameterValues(1) = mu;
        this.ParameterCovariance = zeros (this.NumParameters);
        this.InputData = [];
    endfunction

    function mu = get.mu (this)
        mu = this.ParameterValues(1);
    endfunction

    function this = set.sigma (this, sigma)
        checkparams (this.mu, sigma, this.df)
        this.ParameterValues(2) = sigma;
        this.ParameterCovariance = zeros (this.NumParameters);
        this.InputData = [];
    endfunction

    function sigma = get.sigma (this)
        sigma = this.ParameterValues(2);
    endfunction

    function this = set.df (this, df)
        checkparams (this.mu, this.sigma, df)
        this.ParameterValues(3) = df;
        this.ParameterCovariance = zeros (this.NumParameters);
        this.InputData = [];
    endfunction

    function df = get.df (this)
        df = this.ParameterValues(3);
    endfunction

  endmethods

  methods (Access = public)

    function p = cdf (this, x, uflag)
      p = tlscdf (x, this.mu, this.sigma, this.df, uflag);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        if (strcmpi (uflag, "upper"))
          p(lb) = 1;
          p(ub) = 0;
          p(! (lb | ub)) += tlscdf (lx, this.mu, this.sigma, this.df, uflag);
          p(! (lb | ub)) /= diff (tlscdf ([lx, ux], this.mu, ...
                                          this.sigma, this.df, uflag));
        else
          p(lb) = 0;
          p(ub) = 1;
          p(! (lb | ub)) -= tlscdf (lx, this.mu, this.sigma, this.df);
          p(! (lb | ub)) /= diff (tlscdf ([ux, lx], this.mu, this.sigma, this.df));
        endif
      endif
    endfunction

    function x = icdf (this, p)
      x = tlsinv (p, this.mu, this.sigma, this.df);
      if (this.IsTruncated)
        x(x < this.Truncation(1)) = this.Truncation(1);
        x(x > this.Truncation(2)) = this.Truncation(2);
      endif
    endfunction

    function r = irq (this)
      if (this.IsTruncated)
        lp = ricecdf (this.Truncation(1), this.nu, this.sigma);
        up = ricecdf (this.Truncation(2), this.nu, this.sigma);
        r = diff (tlsinv ([0.75-up, 0.25+lb], this.mu, this.sigma, this.df));
      else
        r = diff (tlsinv ([0.75, 0.25], this.mu, this.sigma, this.df));
      endif
    endfunction

    function m = mean (this)
      m = tlsstat (this.mu, this.sigma, this.df);
    endfunction

    function m = median (this)
      m = tlsstat (this.mu, this.sigma, this.df);
    endfunction

    function nlogL = negloglik (this)
      if (isempty (this.InputData))
        error ("negloglik: no data available");
      endif
      params = [this.mu, this.sigma, this.df];
      nlogL = tlslike (params, this.InputData.x, ...
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
      y = tlspdf (x, this.mu, this.sigma, this.df);
      if (this.IsTruncated)
        lx = this.Truncation(1);
        lb = x < lx;
        ux = this.Truncation(2);
        ub = x > ux;
        y(lb | ub) = 0;
        y(! (lb | ub)) /= diff (tlspdf ([ux, lx], this.mu, this.sigma, this.df));
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
        nsize = 2 * ratio * ps;     # times 1.5 to be on the safe side
        ## Generate the numbers and remove out-of-bound random samples
        r = tlsrnd (this.mu, this.sigma, this.df, nsize, 1);
        r(r < lx | r > ux) = [];
        ## Randomly select the required size and reshape to requested dimensions
        r = randperm (r, ps);
        r = reshape (r, sz);
      else
        r = tlsrnd (this.mu, this.sigma, this.df, varargin{:});
      endif
    endfunction

    function s = std (this)
      [~, v] = tlsstat (this.mu, this.sigma, this.df);
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
      [~, v] = tlsstat (this.mu, this.sigma, this.df);
    endfunction

  endmethods

  methods (Static, Hidden)

    function pd = fit(varargin)
      [paramhat, paramci] = tlsfit (x, 0.05, censor, freq, opt);
      [nlogL, acov] = tlslike (paramhat, x, censor, freq);
      pd = tLocationScaleDistribution.makeFitted ...
           (paramhat, paramci, nlogL, acov, x, censor, freq);
    endfunction

    function pd = makeFitted (paramhat, paramci, nlogL, acov, x, censor, freq)
      mu = paramhat(1);
      sigma = paramhat(2);
      df = paramhat(3);
      pd = tLocationScaleDistribution (mu, sigma, df);
      pd.ParameterCI = paramci;
      pd.NegativeLogLikelihood = nlogL;
      pd.ParameterCovariance = acov;
      pd.InputData = struct ("data", x, "cens", censor, "freq", freq);
    endfunction

  endmethods

endclassdef

function checkparams (mu, sigma, df)
  if (! (isscalar (mu) && isnumeric (mu) && isreal (mu) && isfinite (mu)))
    error ("tLocationScaleDistribution: MU must be a real scalar.")
  endif
  if (! (isscalar (sigma) && isnumeric (sigma) && isreal (sigma)
                          && isfinite (sigma) && sigma > 0))
    error ("tLocationScaleDistribution: SIGMA must be a positive real scalar.")
  endif
  if (! (isscalar (df) && isnumeric (df) && isreal (df)
                       && isfinite (df) && df > 0))
    error ("tLocationScaleDistribution: DF must be a positive real scalar.")
  endif
endfunction
