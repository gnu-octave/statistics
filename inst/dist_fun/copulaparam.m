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

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{param} =} copulaparam (@var{family}, @var{r})
## @deftypefnx {statistics} {@var{param} =} copulaparam (@dots{}, @qcode{"type"}, @var{type})
##
## Copula parameter as a function of rank correlation.
##
## @code{@var{param} = copulaparam (@var{family}, @var{r})} returns the linear
## or copula parameter @var{param} corresponding to a copula of the family
## @var{family} that has Kendall's rank correlation @var{r}.  It is the inverse
## of @code{copulastat}.
##
## @var{family} is the copula family name.  It can be @qcode{"Gaussian"} for the
## Gaussian family, @qcode{"t"} for the Student's t family, @qcode{"Clayton"}
## for the Clayton family, @qcode{"Gumbel"} for the Gumbel-Hougaard family, or
## @qcode{"Frank"} for the Frank family.
##
## For the Gaussian and Student's t families, @var{r} is a scalar rank
## correlation or a matrix of pairwise rank correlations, and @var{param} is the
## corresponding linear correlation of the same size.  For the Clayton,
## Gumbel-Hougaard, and Frank families, @var{r} is a scalar rank correlation and
## @var{param} is the scalar copula parameter.  The Gumbel-Hougaard family
## models positive dependence only, so @var{r} must be non-negative for that
## family.
##
## @code{@var{param} = copulaparam (@dots{}, @qcode{"type"}, @var{type})}
## selects the measure of rank correlation given in @var{r}.  @var{type} can be
## @qcode{"Kendall"} (the default) for Kendall's tau, or @qcode{"Spearman"} for
## Spearman's rho.
##
## The Gaussian and Student's t families and the Kendall's tau of the Clayton
## and Gumbel-Hougaard families are inverted in closed form.  The remaining
## cases are inverted numerically from @code{copulastat}.
##
## @strong{Note:} for the Archimedean families with @qcode{"Spearman"}, the
## underlying relationship is computed by exact numerical integration rather
## than the interpolated table used by @sc{matlab}, so results may differ from
## @sc{matlab} by up to about @math{10^{-4}}.  See @code{copulastat}.
##
## @seealso{copulastat, copulafit, copulacdf, copulapdf, copularnd}
## @end deftypefn

function param = copulaparam (family, r, varargin)

  ## Check arguments
  if (nargin < 2)
    print_usage ();
  endif

  if (! ischar (family))
    error (strcat ("copulaparam: FAMILY must be one of 'Gaussian',", ...
                   " 't', 'Clayton', 'Gumbel', and 'Frank'."));
  endif

  if (! isnumeric (r) || ! isreal (r))
    error ("copulaparam: R must be real.");
  endif

  ## Parse the 'type' option
  type = 'kendall';
  if (numel (varargin) > 0)
    if (numel (varargin) != 2 || ! ischar (varargin{1}) || ...
        ! strcmpi (varargin{1}, 'type'))
      error ("copulaparam: invalid optional argument.");
    endif
    if (! ischar (varargin{2}) || ...
        ! any (strcmpi (varargin{2}, {'kendall', 'spearman'})))
      error ("copulaparam: TYPE must be either 'Kendall' or 'Spearman'.");
    endif
    type = lower (varargin{2});
  endif

  lower_family = lower (family);

  switch (lower_family)

    case {'gaussian', 't'}
      ## Elliptical families: closed-form inverse, applied elementwise
      if (any (abs (r(:)) > 1))
        error ("copulaparam: R must be a correlation in the range -1 to 1.");
      endif
      if (strcmp (type, 'kendall'))
        param = sin (pi .* r ./ 2);
      else
        param = 2 .* sin (pi .* r ./ 6);
      endif

    case {'clayton', 'gumbel', 'frank'}
      ## Archimedean families: scalar rank correlation
      if (! isscalar (r))
        error ("copulaparam: R must be a scalar for the %s family.", family);
      endif
      if (abs (r) >= 1)
        error ("copulaparam: R must be in the range -1 to 1.");
      endif
      if (strcmp (lower_family, 'gumbel') && r < 0)
        error (strcat ("copulaparam: R must be non-negative for the", ...
                       " Gumbel family."));
      endif
      switch (lower_family)
        case 'clayton'
          if (strcmp (type, 'kendall'))
            param = 2 .* r ./ (1 - r);
          else
            param = invert_stat (lower_family, r, type);
          endif
        case 'gumbel'
          if (strcmp (type, 'kendall'))
            param = 1 ./ (1 - r);
          else
            param = invert_stat (lower_family, r, type);
          endif
        case 'frank'
          param = invert_stat (lower_family, r, type);
      endswitch

    otherwise
      error ("copulaparam: unknown copula family '%s'.", family);

  endswitch

endfunction

## Invert a rank correlation numerically from copulastat for an Archimedean
## family.  The dependence measure is monotone increasing in the parameter, so
## we bracket the root and solve.
function a = invert_stat (family, r, type)
  g = @(aa) copulastat (family, aa, 'type', type) - r;
  ## Parameter value at which the family is the independence copula
  if (strcmp (family, 'gumbel'))
    s0 = 1;
  else
    s0 = 0;
  endif
  if (r == 0)
    a = s0;
    return;
  elseif (r > 0)
    lo = s0;
    hi = s0 + 1;
    while (g (hi) < 0 && hi < 1e8)
      hi = s0 + (hi - s0) .* 2;
    endwhile
  else
    ## r < 0 (Clayton or Frank only)
    hi = s0;
    if (strcmp (family, 'clayton'))
      lo = -1 + 1e-8;
    else
      lo = s0 - 1;
      while (g (lo) > 0 && lo > -1e8)
        lo = s0 - (s0 - lo) .* 2;
      endwhile
    endif
  endif
  a = fzero (g, [lo, hi]);
endfunction

%!demo
%! ## Copula parameter of a Gaussian copula with Kendall's tau 0.3
%! rho = copulaparam ("Gaussian", 0.3)

%!demo
%! ## copulaparam inverts copulastat
%! alpha = copulaparam ("Clayton", 0.5)
%! tau = copulastat ("Clayton", alpha)

## Test output against MATLAB
%!test
%! assert_equal (copulaparam ("Gaussian", 0.3), 0.453990499739547, 1e-14);
%! assert_equal (copulaparam ("Gaussian", 0.3, "type", "Spearman"), ...
%!               0.312868930080462, 1e-14);
%! assert_equal (copulaparam ("t", 0.3), 0.453990499739547, 1e-14);
%! assert_equal (copulaparam ("Clayton", 0.3), 0.857142857142857, 1e-14);
%! assert_equal (copulaparam ("Frank", 0.3), 2.91743444592452, 1e-8);
%! assert_equal (copulaparam ("Gumbel", 0.3), 1.42857142857143, 1e-13);

## copulaparam is the inverse of copulastat (round trip)
%!test
%! for tau = [0.1, 0.25, 0.5, 0.7]
%!   for fam = {"Clayton", "Gumbel", "Frank"}
%!     a = copulaparam (fam{1}, tau);
%!     assert_equal (copulastat (fam{1}, a), tau, 1e-8);
%!   endfor
%! endfor
%!test
%! for rs = [0.1, 0.3, 0.6]
%!   for fam = {"Clayton", "Gumbel", "Frank"}
%!     a = copulaparam (fam{1}, rs, "type", "Spearman");
%!     assert_equal (copulastat (fam{1}, a, "type", "Spearman"), rs, 1e-7);
%!   endfor
%! endfor

## Negative dependence for the signed families
%!test
%! a = copulaparam ("Frank", -0.3);
%! assert_equal (copulastat ("Frank", a), -0.3, 1e-8);
%! a = copulaparam ("Clayton", -0.2);
%! assert_equal (copulastat ("Clayton", a), -0.2, 1e-8);

## Elliptical families accept a correlation matrix and act elementwise
%!test
%! tau = [1, 0.3; 0.3, 1];
%! assert_equal (copulaparam ("Gaussian", tau), sin (pi .* tau ./ 2), 1e-14);

## Test input validation
%!error <copulaparam: FAMILY must be one of 'Gaussian', 't', 'Clayton', 'Gumbel', and 'Frank'.> ...
%! copulaparam (5, 0.3)
%!error <copulaparam: R must be real.> copulaparam ("Gaussian", 2i)
%!error <copulaparam: TYPE must be either 'Kendall' or 'Spearman'.> ...
%! copulaparam ("Gaussian", 0.3, "type", "Pearson")
%!error <copulaparam: invalid optional argument.> ...
%! copulaparam ("Gaussian", 0.3, "foo")
%!error <copulaparam: R must be a correlation in the range -1 to 1.> ...
%! copulaparam ("Gaussian", 1.5)
%!error <copulaparam: R must be a scalar for the Clayton family.> ...
%! copulaparam ("Clayton", [0.1, 0.2])
%!error <copulaparam: R must be in the range -1 to 1.> ...
%! copulaparam ("Clayton", 1)
%!error <copulaparam: R must be non-negative for the Gumbel family.> ...
%! copulaparam ("Gumbel", -0.3)
%!error <copulaparam: unknown copula family 'Foo'.> copulaparam ("Foo", 0.3)
