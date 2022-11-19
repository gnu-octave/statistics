## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

## -*- texinfo -*-
## @deftypefn {Function File} @var{p} = nctcdf (@var{x}, @var{df}, @var{delta})
## @deftypefnx {Function File} @var{p} = nctcdf (@var{x}, @var{df}, @var{delta}, @var{uflag})
##
## Noncentral T cumulative distribution function (cdf).
##
## @code{@var{p} = nctcdf (@var{x}, @var{df}, @var{delta})} returns the
## noncentral T cdf with @var{df} degrees of freedom and noncentrality parameter
## @var{delta} at the values of @var{X}.
##
## The size of @var{p} is the common size of the input arguments. Scalar input
## arguments @var{x}, @var{df}, @var{delta} are regarded as constant matrices of
## the same size as the other inputs.
##
## @code{@var{p} = nctcdf (@var{x}, @var{df}, @var{delta}, "upper"} returns the
## upper tail probability of the noncentral T distribution with @var{df} degrees
## of freedom and noncentrality parameter @var{delta} at the values in @var{x}.
##
## @seealso{nctinv, nctpdf, nctrnd, nctstat}
## @end deftypefn

function p = nctcdf (x, df, delta, uflag)

  ## Check for valid input arguments
  if (nargin <  3)
    error ("nctcdf: too few imputs.");
  endif

  ## Check and fix size of input arguments
  [err, x, df, delta] = common_size (x, df, delta);
  if (err > 0)
    error ("nctcdf: input size mismatch.");
  endif

  ## Check for upper tail option
  if (nargin > 3)
    if (! strcmpi (uflag, "upper"))
      error ("nctcdf: improper definition of upper tail option.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Initialize p
  if (isa (x, "single") || isa (df, "single") || isa (delta, "single"))
    p = zeros (size (x), "single");
    c_eps = eps ("single");
  else
    p = zeros (size (x));
    c_eps = eps;
  endif

  ## Find NaNs in input arguments (if any) and propagate them to p
  is_nan = isnan (x) | isnan (df) | isnan (delta);
  p(is_nan) = NaN;

  ## Find special cases for delta==0 and x<0; and x = Inf.
  case_Dinf = (df <= 0 | isinf(delta)) & ! is_nan;
  case_Dzero = delta == 0 & ! case_Dinf & ! is_nan;
  case_Xzero = x < 0 & ! case_Dzero & ! case_Dinf & ! is_nan;
  case_Xinf = x == Inf & ! case_Dzero & ! case_Dinf & ! is_nan;
  case_DFbig = df > 2e6 & ! case_Dzero & ! case_Dinf & ! case_Xinf & ! is_nan;
  flag_Dinf = any (case_Dinf(:));
  flag_Dzero = any (case_Dzero(:));
  flag_Xzero = any (case_Xzero(:));
  flag_Xinf = any (case_Xinf(:));
  flag_DFbig = any (case_DFbig(:));

  ## Handle special cases
  if (flag_Dinf || flag_Dzero || flag_Xzero || flag_Xinf || flag_DFbig)
    if (flag_Dinf)
      p(case_Dinf) = NaN;
    endif
    if (flag_Dzero)
      if (uflag)
        p(case_Dzero) = tcdf (x(case_Dzero), df(case_Dzero), "upper");
      else
        p(case_Dzero) = tcdf (x(case_Dzero), df(case_Dzero));
      endif
    endif
    if (flag_Xinf)
      if (uflag)
        p(case_Xinf) = 0;
      else
        p(case_Xinf) = 1;
      endif
    endif
    if (flag_DFbig)
      s = 1 - 1 ./ (4 * df);
      d = sqrt (1 + x .^ 2 ./ (2 * df));
      if (uflag)
        p(case_DFbig) = normcdf (x(case_DFbig) .* s(case_DFbig), ...
                                 delta(case_DFbig), d(case_DFbig), "upper");
      else
        p(case_DFbig) = normcdf (x(case_DFbig) .* s(case_DFbig), ...
                                 delta(case_DFbig), d(case_DFbig));
      endif
    endif
    fp = ! (case_Dinf | case_Dzero | case_Xzero | case_Xinf | case_DFbig);
    if (any (fp(:)))
      if (uflag)
        p(fp) = nctcdf (x(fp), df(fp), delta(fp), "upper");
      else
        p(fp) = nctcdf (x(fp), df(fp), delta(fp));
      endif
    endif
    if (flag_Xzero)
      if (uflag)
        p(case_Xzero) = nctcdf (-x(case_Xzero), df(case_Xzero), ...
                                -delta(case_Xzero));
      else
        p(case_Xzero) = nctcdf (-x(case_Xzero), df(case_Xzero), ...
                                -delta(case_Xzero), "upper");
      endif
    endif
    return
  endif

  ## Compute value for betainc function.
  x_square = x .^ 2;
  denom = df + x_square;
  P = x_square ./ denom;
  Q = df ./ denom;
  ## Initialize infinite sum.
  d_square = delta .^ 2;

  ## Compute probability P[TD<0] (first term)
  if (uflag)
    x_zero = x == 0 & ! is_nan;
    if (any (x_zero(:)))
      fx = normcdf (- delta, 0, 1, "upper");
      p(x_zero)= fx(x_zero);
    endif
  else
    p(! is_nan) = normcdf (- delta(! is_nan), 0, 1);
  endif

  ## Compute probability P[0<TD<x] (second term)
  x_notzero = find (x != 0 & ! is_nan);
  if (any (x_notzero(:)))
    P = P(x_notzero);
    Q = Q(x_notzero);
    df = df(x_notzero);
    d_square = d_square(x_notzero);
    d_sign = sign (delta(x_notzero));
    subtotal = zeros (size (x_notzero));

    ## Start looping over term jj and higher, this should be near the
    ## peak of the E part of the term (see below)
    jj = 2 * floor(d_square/2);

    ## Compute an infinite sum using Johnson & Kotz eq 9, or new
    ## edition eq 31.16, each term having this form:
    ##  B  = betainc(P,(j+1)/2,df/2);
    ##  E  = (exp(0.5*j*log(0.5*delta^2) - gammaln(j/2+1)));
    ##  term = E .* B;
    ##
    ## We'll compute betainc at the beginning, and then update using
    ## recurrence formulas (Abramowitz & Stegun 26.5.16).  We'll sum the
    ## series two terms at a time to make the recurrence work out.
    jj = 2 * floor (d_square / 2);
    E1 = exp (0.5 * jj .* log (0.5 * d_square) - ...
              d_square / 2 - gammaln (jj / 2 + 1));
    E2 = d_sign .* exp (0.5 * (jj + 1) .* log (0.5 * d_square) - ...
                        d_square / 2 - gammaln ((jj + 1) / 2 + 1));

    ## Use either P or Q, whichever is more accurately computed
    TD = (P < 0.5);
    B1 = zeros (size (P));
    B2 = zeros (size (P));
    if (uflag)
      if any(TD)
        B1(TD) = betainc (P(TD), (jj(TD) + 1) / 2, df(TD) / 2, "upper");
        B2(TD) = betainc (P(TD), (jj(TD) + 2) / 2, df(TD) / 2, "upper");
      endif
      TD = ! TD;
      if (any (TD))
        B1(TD) = betainc (Q(TD), df(TD) / 2, (jj(TD) + 1) / 2, "lower");
        B2(TD) = betainc (Q(TD), df(TD) / 2, (jj(TD) + 2) / 2, "lower");
      endif
    else
      if (any (TD))
        B1(TD) = betainc (P(TD), (jj(TD) + 1) / 2, df(TD) / 2, "lower");
        B2(TD) = betainc (P(TD), (jj(TD) + 2) / 2, df(TD) / 2, "lower");
      endif
      TD = ! TD;
      if (any (TD))
        B1(TD) = betainc (Q(TD), df(TD) / 2, (jj(TD) + 1) / 2, "upper");
        B2(TD) = betainc (Q(TD), df(TD) / 2, (jj(TD) + 2) / 2, "upper");
      endif
    endif
    R1 = exp (gammaln ((jj + 1) / 2 + df / 2) - gammaln ((jj + 3) / 2) - ...
             gammaln (df / 2) + ((jj + 1) / 2) .* log (P) + (df / 2) .* log(Q));
    R2 = exp (gammaln ((jj + 2) / 2 + df / 2) - gammaln ((jj + 4) / 2) - ...
             gammaln (df / 2) + ((jj + 2) / 2) .* log (P) + (df / 2) .* log(Q));

    ## Keep terms
    E10 = E1; E20 = E2; B10 = B1; B20 = B2; R10 = R1; R20 = R2; j0 = jj;
    TD = true (size (d_square));
    while(true)
      ## Probability that TD lies between 0 and x
      twoterms = E1(TD) .* B1(TD) + E2(TD) .* B2(TD);
      subtotal(TD) = subtotal(TD) + twoterms;
      ## Convergence test.
      TD(TD) = (abs(twoterms) > (abs (subtotal(TD)) + c_eps) * c_eps);
      if (! any (TD))
        break;
      end
      ## Update for next iteration
      jj = jj+2;
      E1(TD) = E1(TD) .* d_square(TD) ./ (jj(TD));
      E2(TD) = E2(TD) .* d_square(TD) ./ (jj(TD) + 1);
      if (uflag)
        B1(TD) = betainc (P(TD), (jj(TD) + 1) / 2, df(TD) / 2, "upper");
        B2(TD) = betainc (P(TD), (jj(TD) + 2) / 2, df(TD) / 2, "upper");
      else
        B1(TD) = B1(TD) - R1(TD);
        B2(TD) = B2(TD) - R2(TD);
        R1(TD) = R1(TD) .* P(TD) .* (jj(TD)+df(TD)-1) ./ (jj(TD)+1);
        R2(TD) = R2(TD) .* P(TD) .* (jj(TD)+df(TD)  ) ./ (jj(TD)+2);
      endif
    endwhile

    ## Go back to the peak and start looping downward as far as necessary.
    E1 = E10; E2 = E20; B1 = B10; B2 = B20; R1 = R10; R2 = R20;
    jj = j0;
    TD = (jj > 0);
    while (any (TD))
      JJ = jj(TD);
      E1(TD) = E1(TD) .* (JJ  ) ./ d_square(TD);
      E2(TD) = E2(TD) .* (JJ+1) ./ d_square(TD);
      R1(TD) = R1(TD) .* (JJ+1) ./ ((JJ+df(TD)-1) .* P(TD));
      R2(TD) = R2(TD) .* (JJ+2) ./ ((JJ+df(TD))   .* P(TD));
      if (uflag)
        B1(TD) = betainc (P(TD), (JJ - 1) / 2, df(TD) / 2, "upper");
        B2(TD) = betainc (P(TD), JJ / 2, df(TD) / 2, "upper");
      else
        B1(TD) = B1(TD) + R1(TD);
        B2(TD) = B2(TD) + R2(TD);
      end
      twoterms = E1(TD) .* B1(TD) + E2(TD) .* B2(TD);
      subtotal(TD) = subtotal(TD) + twoterms;
      jj = jj - 2;
      TD(TD) = (abs (twoterms) > (abs (subtotal(TD)) + c_eps) * c_eps) & ...
               (jj(TD) > 0);
    endwhile
    p(x_notzero) = min (1, max (0, p(x_notzero) + subtotal / 2));
  endif

endfunction

%!demo
%! ## Compare the noncentral t cdf with DELTA = 1 to the t cdf
%! ## with the same number of degrees of freedom (10).
%!
%! x = (-5:0.1:5)';
%! p1 = nctcdf (x, 10, 1);
%! p = tcdf (x, 10);
%! plot (x, p, "-", x, p1, ":")

## Input validation tests
%!error<nctcdf: too few imputs.> p = nctcdf (2, 4);
%!error<nctcdf: input size mismatch.> p = nctcdf (2, [4, 3], [3, 4, 5]);
%!error<nctcdf: improper definition of upper tail option.> ...
%! p = nctcdf (2, 4, 2, "lower");

## Output validation tests
%!test
%! x = (-2:0.1:2)';
%! p = nctcdf (x, 10, 1);
%! assert (p(1), 0.003302485766631558, 1e-14);
%! assert (p(2), 0.004084668193532631, 1e-14);
%! assert (p(3), 0.005052800319478737, 1e-14);
%! assert (p(41), 0.8076115625303751, 1e-14);
%!test
%! p = nctcdf (12, 10, 3);
%! assert (p, 0.9997719343243797, 1e-14);
%!test
%! p = nctcdf (2, 3, 2);
%! assert (p, 0.4430757822176028, 1e-14);
%!test
%! p = nctcdf (2, 3, 2, "upper");
%! assert (p, 0.5569242177823971, 1e-14);
%!test
%! p = nctcdf ([3, 6], 3, 2, "upper");
%! assert (p, [0.3199728259444777, 0.07064855592441913], 1e-14);


