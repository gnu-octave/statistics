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

classdef paretotails
  ## -*- texinfo -*-
  ## @deftp {statistics} paretotails
  ##
  ## Piecewise distribution with generalized Pareto tails.
  ##
  ## A @code{paretotails} object is a piecewise probability distribution fit to
  ## sample data.  A generalized Pareto distribution (GPD) is fit to each tail
  ## of the data, below a lower quantile and above an upper quantile, while the
  ## middle of the distribution is described by the empirical cumulative
  ## distribution function of the data.  This gives a smooth model for the
  ## tails, useful for extreme value analysis, together with a nonparametric
  ## description of the central region.
  ##
  ## Create a @code{paretotails} object with the constructor
  ## @code{@var{pt} = paretotails (@var{x}, @var{pl}, @var{pu})}, where @var{x}
  ## is the sample data and @var{pl} and @var{pu} are the cumulative
  ## probabilities at the lower and upper tail boundaries.  Data at or below the
  ## @var{pl} quantile form the lower tail, data at or above the @var{pu}
  ## quantile form the upper tail, and the rest form the middle segment.
  ##
  ## Query the fitted object with the methods @code{cdf}, @code{pdf},
  ## @code{icdf}, @code{random}, @code{boundary}, @code{nsegments},
  ## @code{segment}, @code{lowerparams}, and @code{upperparams}.
  ##
  ## @strong{Note:} the kernel-smoothed middle option of @sc{matlab}
  ## (@code{paretotails (@var{x}, @var{pl}, @var{pu}, "kernel")}) is not yet
  ## supported; only the default empirical (@qcode{"ecdf"}) middle is available.
  ##
  ## @seealso{gpfit, gpcdf, gppdf, gpinv, ecdf, fitdist,
  ## GeneralizedParetoDistribution}
  ## @end deftp

  properties (SetAccess = private)
    ## -*- texinfo -*-
    ## @deftp {paretotails} {property} NumSegments
    ##
    ## Number of segments in the piecewise distribution (a lower tail, a middle,
    ## and an upper tail give three).
    ##
    ## @end deftp
    NumSegments = 3

    ## -*- texinfo -*-
    ## @deftp {paretotails} {property} NumParameters
    ##
    ## Number of estimated parameters (two per fitted generalized Pareto tail).
    ##
    ## @end deftp
    NumParameters = 4
  endproperties

  properties (Access = private)
    lowerP = []          # [k, sigma] of the lower-tail GPD, or [] if no tail
    upperP = []          # [k, sigma] of the upper-tail GPD, or [] if no tail
    pl = 0               # cumulative probability at the lower boundary
    pu = 1               # cumulative probability at the upper boundary
    ql = -Inf            # lower boundary quantile
    qu = Inf             # upper boundary quantile
    xknot = []           # middle-segment interpolation knots (x)
    pknot = []           # middle-segment interpolation knots (cumulative prob)
    nobs = 0             # number of (non-NaN) observations
    method = "ecdf"      # description of the middle segment
  endproperties

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {paretotails} {@var{pt} =} paretotails (@var{x}, @var{pl}, @var{pu})
    ## @deftypefnx {paretotails} {@var{pt} =} paretotails (@var{x}, @var{pl}, @var{pu}, @var{cdffun})
    ##
    ## Fit a piecewise distribution with generalized Pareto tails to @var{x}.
    ##
    ## @var{pl} and @var{pu} are the cumulative probabilities of the lower and
    ## upper tail boundaries, with @code{0 <= @var{pl} < @var{pu} <= 1}.  A
    ## generalized Pareto distribution is fit by maximum likelihood to the
    ## exceedances in each tail; the middle segment uses the empirical
    ## cumulative distribution of @var{x}.
    ##
    ## @var{cdffun} selects the middle segment and defaults to @qcode{"ecdf"};
    ## the @qcode{"kernel"} option of @sc{matlab} is not yet supported.
    ##
    ## @end deftypefn
    function this = paretotails (x, pl, pu, cdffun)

      if (nargin < 3)
        print_usage ();
      endif

      if (! (isnumeric (x) && isreal (x) && isvector (x)))
        error ("paretotails: X must be a numeric vector.");
      endif
      x = x(:);
      x(isnan (x)) = [];
      if (numel (x) < 2)
        error ("paretotails: X must contain at least two non-NaN values.");
      endif

      if (! (isnumeric (pl) && isscalar (pl) && isreal (pl) ...
             && pl >= 0 && pl <= 1))
        error ("paretotails: PL must be a scalar in the range [0, 1].");
      endif
      if (! (isnumeric (pu) && isscalar (pu) && isreal (pu) ...
             && pu >= 0 && pu <= 1))
        error ("paretotails: PU must be a scalar in the range [0, 1].");
      endif
      if (pl >= pu)
        error ("paretotails: PL must be less than PU.");
      endif

      if (nargin > 3)
        if (! ischar (cdffun))
          error ("paretotails: CDFFUN must be 'ecdf' or 'kernel'.");
        endif
        switch (lower (cdffun))
          case "ecdf"
            ## default
          case "kernel"
            error (strcat ("paretotails: the 'kernel' middle segment is", ...
                           " not yet supported; use the default 'ecdf'."));
          otherwise
            error ("paretotails: CDFFUN must be 'ecdf' or 'kernel'.");
        endswitch
      endif

      xs = sort (x);
      n = numel (xs);
      pp = ((1:n)' - 0.5) ./ n;

      this.pl = pl;
      this.pu = pu;
      this.nobs = n;
      this.ql = paretotails.quantile_pp (xs, pp, pl);
      this.qu = paretotails.quantile_pp (xs, pp, pu);

      ## Fit the generalized Pareto tails to the exceedances (gpfit returns
      ## [k, sigma, theta]; theta is the fixed threshold 0, so keep [k, sigma]).
      if (pl > 0)
        lp = gpfit (this.ql - xs(xs < this.ql), 0);
        this.lowerP = lp(1:2);
      endif
      if (pu < 1)
        up = gpfit (xs(xs > this.qu) - this.qu, 0);
        this.upperP = up(1:2);
      endif

      ## Middle-segment interpolation knots: the interior data points, with the
      ## two boundary points prepended and appended so the piecewise-linear cdf
      ## maps [ql, qu] onto [pl, pu] exactly.
      mid = xs > this.ql & xs < this.qu;
      xk = xs(mid);
      pk = pp(mid);
      if (pl > 0)
        xk = [this.ql; xk];
        pk = [pl; pk];
      endif
      if (pu < 1)
        xk = [xk; this.qu];
        pk = [pk; pu];
      endif
      ## Drop any repeated abscissae so interp1 sees a strictly increasing grid
      keep = [true; diff(xk) > 0];
      this.xknot = xk(keep);
      this.pknot = pk(keep);

      this.NumSegments = 1 + (pl > 0) + (pu < 1);
      this.NumParameters = 2 .* (pl > 0) + 2 .* (pu < 1);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {paretotails} {@var{p} =} cdf (@var{pt}, @var{x})
    ##
    ## Cumulative distribution function of the @code{paretotails} object
    ## @var{pt} evaluated at the values in @var{x}.
    ##
    ## @end deftypefn
    function p = cdf (this, x)
      if (nargin != 2)
        print_usage ();
      endif
      p = nan (size (x));
      lo = x < this.ql;
      hi = x > this.qu;
      mid = ! lo & ! hi;
      if (any (lo(:)))
        p(lo) = this.pl .* (1 - gpcdf (this.ql - x(lo), ...
                                       this.lowerP(1), this.lowerP(2), 0));
      endif
      if (any (hi(:)))
        p(hi) = this.pu + (1 - this.pu) .* gpcdf (x(hi) - this.qu, ...
                                       this.upperP(1), this.upperP(2), 0);
      endif
      if (any (mid(:)))
        p(mid) = interp1 (this.xknot, this.pknot, x(mid), "linear");
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {paretotails} {@var{y} =} pdf (@var{pt}, @var{x})
    ##
    ## Probability density function of the @code{paretotails} object @var{pt}
    ## evaluated at the values in @var{x}.
    ##
    ## @end deftypefn
    function y = pdf (this, x)
      if (nargin != 2)
        print_usage ();
      endif
      y = nan (size (x));
      lo = x < this.ql;
      hi = x > this.qu;
      mid = ! lo & ! hi;
      if (any (lo(:)))
        y(lo) = this.pl .* gppdf (this.ql - x(lo), ...
                                  this.lowerP(1), this.lowerP(2), 0);
      endif
      if (any (hi(:)))
        y(hi) = (1 - this.pu) .* gppdf (x(hi) - this.qu, ...
                                        this.upperP(1), this.upperP(2), 0);
      endif
      if (any (mid(:)))
        ## The middle cdf is piecewise linear, so its density is the piecewise
        ## constant slope of the segment each point falls in.
        slope = diff (this.pknot) ./ diff (this.xknot);
        idx = paretotails.bin_index (this.xknot, x(mid));
        y(mid) = slope(idx);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {paretotails} {@var{x} =} icdf (@var{pt}, @var{p})
    ##
    ## Inverse cumulative distribution function (quantile function) of the
    ## @code{paretotails} object @var{pt} evaluated at the probabilities
    ## @var{p}.
    ##
    ## @end deftypefn
    function x = icdf (this, p)
      if (nargin != 2)
        print_usage ();
      endif
      x = nan (size (p));
      valid = p >= 0 & p <= 1;
      lo = valid & p < this.pl;
      hi = valid & p > this.pu;
      mid = valid & ! lo & ! hi;
      if (any (lo(:)))
        x(lo) = this.ql - gpinv (1 - p(lo) ./ this.pl, ...
                                 this.lowerP(1), this.lowerP(2), 0);
      endif
      if (any (hi(:)))
        x(hi) = this.qu + gpinv ((p(hi) - this.pu) ./ (1 - this.pu), ...
                                 this.upperP(1), this.upperP(2), 0);
      endif
      if (any (mid(:)))
        x(mid) = interp1 (this.pknot, this.xknot, p(mid), "linear");
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {paretotails} {@var{r} =} random (@var{pt})
    ## @deftypefnx {paretotails} {@var{r} =} random (@var{pt}, @var{sz})
    ## @deftypefnx {paretotails} {@var{r} =} random (@var{pt}, @var{m}, @var{n}, @dots{})
    ##
    ## Random values drawn from the @code{paretotails} object @var{pt}, by
    ## inverse transform sampling.  The size arguments follow @code{rand}.
    ##
    ## @end deftypefn
    function r = random (this, varargin)
      r = this.icdf (rand (varargin{:}));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {paretotails} {[@var{p}, @var{q}] =} boundary (@var{pt})
    ##
    ## Boundary probabilities @var{p} and quantiles @var{q} of the segments of
    ## the @code{paretotails} object @var{pt}, as column vectors.
    ##
    ## @end deftypefn
    function [p, q] = boundary (this)
      p = [];
      q = [];
      if (this.pl > 0)
        p = [p; this.pl];
        q = [q; this.ql];
      endif
      if (this.pu < 1)
        p = [p; this.pu];
        q = [q; this.qu];
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {paretotails} {@var{n} =} nsegments (@var{pt})
    ##
    ## Number of segments in the @code{paretotails} object @var{pt}.
    ##
    ## @end deftypefn
    function n = nsegments (this)
      n = this.NumSegments;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {paretotails} {@var{params} =} lowerparams (@var{pt})
    ##
    ## Shape and scale parameters @code{[@var{k}, @var{sigma}]} of the
    ## generalized Pareto distribution fit to the lower tail of @var{pt}.
    ##
    ## @end deftypefn
    function params = lowerparams (this)
      if (isempty (this.lowerP))
        error ("paretotails: the distribution has no lower tail.");
      endif
      params = this.lowerP(:)';
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {paretotails} {@var{params} =} upperparams (@var{pt})
    ##
    ## Shape and scale parameters @code{[@var{k}, @var{sigma}]} of the
    ## generalized Pareto distribution fit to the upper tail of @var{pt}.
    ##
    ## @end deftypefn
    function params = upperparams (this)
      if (isempty (this.upperP))
        error ("paretotails: the distribution has no upper tail.");
      endif
      params = this.upperP(:)';
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {paretotails} {@var{s} =} segment (@var{pt}, @var{x}, @var{p})
    ##
    ## Segment indices for the @code{paretotails} object @var{pt}.  Supply the
    ## data values in @var{x} (with @var{p} empty) or the cumulative
    ## probabilities in @var{p} (with @var{x} empty).  Segment @code{1} is the
    ## lower tail, @code{2} the middle, and @code{3} the upper tail.
    ##
    ## @end deftypefn
    function s = segment (this, x, p)
      if (nargin < 3)
        p = [];
      endif
      if (nargin < 2)
        x = [];
      endif
      if (isempty (x) && ! isempty (p))
        s = 2 .* ones (size (p));
        s(p < this.pl) = 1;
        s(p > this.pu) = 3;
      else
        s = 2 .* ones (size (x));
        s(x < this.ql) = 1;
        s(x > this.qu) = 3;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {paretotails} {} disp (@var{pt})
    ##
    ## Display a summary of the @code{paretotails} object @var{pt}.
    ##
    ## @end deftypefn
    function disp (this)
      printf ("  Piecewise distribution with %d segments fit to %d ", ...
              this.NumSegments, this.nobs);
      printf ("observations:\n\n");
      if (this.pl > 0)
        printf ("    P = [0, %g]: lower tail, generalized Pareto ", this.pl);
        printf ("(k = %g, sigma = %g)\n", this.lowerP(1), this.lowerP(2));
      endif
      printf ("    P = [%g, %g]: middle, %s\n", this.pl, this.pu, this.method);
      if (this.pu < 1)
        printf ("    P = [%g, 1]: upper tail, generalized Pareto ", this.pu);
        printf ("(k = %g, sigma = %g)\n", this.upperP(1), this.upperP(2));
      endif
    endfunction

  endmethods

  methods (Static, Access = private)

    ## Empirical quantile by linear interpolation of the plotting positions
    ## (i-0.5)/n, clamped to the data range (the default of MATLAB's quantile).
    function q = quantile_pp (xs, pp, p)
      if (p <= pp(1))
        q = xs(1);
      elseif (p >= pp(end))
        q = xs(end);
      else
        q = interp1 (pp, xs, p, "linear");
      endif
    endfunction

    ## Index of the interpolation bin containing each value of x, clamped to the
    ## valid range of the knot vector xk.
    function idx = bin_index (xk, x)
      idx = zeros (size (x));
      for i = 1:numel (x)
        j = find (xk <= x(i), 1, "last");
        if (isempty (j))
          j = 1;
        elseif (j >= numel (xk))
          j = numel (xk) - 1;
        endif
        idx(i) = j;
      endfor
    endfunction

  endmethods

endclassdef

%!demo
%! ## Fit Pareto tails to a normal sample and compare the tail cdf to the data
%! x = norminv (((1:100) - 0.5) / 100);
%! pt = paretotails (x, 0.1, 0.9);
%! lowerparams (pt)
%! [p, q] = boundary (pt)
%! cdf (pt, [-2.5, 0, 2.5])

## Shared probe data and MATLAB reference values (x = norminv of 100 plotting
## positions; pl = 0.1, pu = 0.9).
%!shared x, pt
%! x = norminv (((1:100) - 0.5) / 100);
%! pt = paretotails (x, 0.1, 0.9);

%!test
%! assert_equal (lowerparams (pt), ...
%!               [-0.381277950146653, 0.652296030248444], 1e-5);
%! assert_equal (upperparams (pt), ...
%!               [-0.381277950146658, 0.652296030248448], 1e-5);
%!test
%! [p, q] = boundary (pt);
%! assert_equal (p, [0.1; 0.9], 1e-12);
%! assert_equal (q, [-1.28207227531929; 1.28207227531929], 1e-10);
%!test
%! assert_equal (nsegments (pt), 3);
%! assert_equal (pt.NumSegments, 3);

%!test
%! xq = [-2.8, -2.3, -1.8, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 1.8, 2.3, 2.8];
%! ref = [0.000326503983434112, 0.00934245750239852, 0.0388388580428329, ...
%!        0.0699512617277107, 0.158702922884361, 0.308553672872447, 0.5, ...
%!        0.691446327127553, 0.84129707711564, 0.930048738272289, ...
%!        0.961161141957167, 0.990657542497602, 0.999673496016566];
%! assert_equal (cdf (pt, xq), ref, 1e-5);
%!test
%! xq = [-2.8, -2.3, -1.8, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 1.8, 2.3, 2.8];
%! ref = [0.00443959395367343, 0.0353636322020254, 0.0853936076680005, ...
%!        0.122892916352678, 0.243260728154123, 0.352775902406303, ...
%!        0.398931835816165, 0.352775902406307, 0.243260728154123, ...
%!        0.122892916352677, 0.0853936076680006, 0.0353636322020256, ...
%!        0.00443959395367328];
%! assert_equal (pdf (pt, xq), ref, 1e-5);
%!test
%! pq = [0.005, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.98, ...
%!       0.99, 0.995];
%! ref = [-2.44694213045116, -2.28179640154802, -2.06669489811039, ...
%!        -1.67939673030492, -1.28207227531929, -0.674573258334612, 0, ...
%!        0.67457325833461, 1.28207227531929, 1.67939673030492, ...
%!        2.06669489811039, 2.28179640154802, 2.44694213045116];
%! assert_equal (icdf (pt, pq), ref, 1e-5);
%!test
%! xq = [-2.8, -2.3, -1.8, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 1.8, 2.3, 2.8];
%! ref = [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3];
%! assert_equal (segment (pt, xq, []), ref);

## cdf and icdf are inverses in the middle and tails
%!test
%! p = [0.02, 0.2, 0.5, 0.8, 0.98];
%! assert_equal (cdf (pt, icdf (pt, p)), p, 1e-6);

## random returns values of the requested size within the support
%!test
%! r = random (pt, 3, 4);
%! assert_equal (size (r), [3, 4]);

## Boundary probabilities outside [0,1] give NaN from icdf
%!test
%! assert_equal (icdf (pt, [-0.1, 1.1]), [NaN, NaN]);

## Test input validation
%!error <Invalid call to paretotails> paretotails (1)
%!error <paretotails: X must be a numeric vector.> paretotails ("a", 0.1, 0.9)
%!error <paretotails: X must contain at least two non-NaN values.> ...
%! paretotails (1, 0.1, 0.9)
%!error <paretotails: PL must be a scalar in the range> ...
%! paretotails (1:10, -0.1, 0.9)
%!error <paretotails: PU must be a scalar in the range> ...
%! paretotails (1:10, 0.1, 1.5)
%!error <paretotails: PL must be less than PU.> paretotails (1:10, 0.9, 0.1)
%!error <paretotails: the 'kernel' middle segment is not yet supported> ...
%! paretotails (1:10, 0.1, 0.9, "kernel")
%!error <paretotails: CDFFUN must be 'ecdf' or 'kernel'.> ...
%! paretotails (1:10, 0.1, 0.9, "foo")
