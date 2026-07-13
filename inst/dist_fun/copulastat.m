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
## @deftypefn  {statistics} {@var{r} =} copulastat (@var{family}, @var{param})
## @deftypefnx {statistics} {@var{r} =} copulastat (@dots{}, @qcode{"type"}, @var{type})
##
## Rank correlation for a copula family.
##
## @code{@var{r} = copulastat (@var{family}, @var{param})} returns Kendall's
## rank correlation @var{r} corresponding to a copula of the family
## @var{family} with linear or copula parameter @var{param}.
##
## @var{family} is the copula family name.  It can be @qcode{"Gaussian"} for the
## Gaussian family, @qcode{"t"} for the Student's t family, @qcode{"Clayton"}
## for the Clayton family, @qcode{"Gumbel"} for the Gumbel-Hougaard family, or
## @qcode{"Frank"} for the Frank family.
##
## For the Gaussian and Student's t families, @var{param} is a linear
## correlation coefficient @var{rho} in the range @math{[-1,1]}, or a
## @math{p}-by-@math{p} correlation matrix, in which case @var{r} has the same
## size and each element is computed elementwise.  For the Clayton,
## Gumbel-Hougaard, and Frank families, @var{param} is the scalar copula
## parameter.
##
## @code{@var{r} = copulastat (@dots{}, @qcode{"type"}, @var{type})} selects the
## measure of rank correlation.  @var{type} can be @qcode{"Kendall"} (the
## default) for Kendall's tau, or @qcode{"Spearman"} for Spearman's rho.
##
## The relationships are closed form for the Gaussian and Student's t families
## (@math{r = 2 \arcsin(\rho) / \pi} for Kendall's tau and
## @math{r = 6 \arcsin(\rho/2) / \pi} for Spearman's rho) and for the Kendall's
## tau of the Archimedean families.  Spearman's rho of the Archimedean families
## has no closed form and is computed by accurate numerical integration of the
## copula.
##
## @strong{Note:} @sc{matlab} returns the Archimedean Spearman's rho by
## interpolating an internal precomputed table, whose values deviate from the
## true relationship by up to about @math{10^{-4}}.  This implementation returns
## the mathematically exact value instead, so results for
## @code{copulastat (@var{family}, @var{param}, "type", "Spearman")} with an
## Archimedean @var{family} may differ from @sc{matlab} at that level.
##
## @seealso{copulaparam, copulafit, copulacdf, copulapdf, copularnd}
## @end deftypefn

function r = copulastat (family, param, varargin)

  ## Check arguments
  if (nargin < 2)
    print_usage ();
  endif

  if (! ischar (family))
    error (strcat ("copulastat: FAMILY must be one of 'Gaussian',", ...
                   " 't', 'Clayton', 'Gumbel', and 'Frank'."));
  endif

  if (! isnumeric (param) || ! isreal (param))
    error ("copulastat: PARAM must be real.");
  endif

  ## Parse the 'type' option
  type = 'kendall';
  if (numel (varargin) > 0)
    if (numel (varargin) != 2 || ! ischar (varargin{1}) || ...
        ! strcmpi (varargin{1}, 'type'))
      error ("copulastat: invalid optional argument.");
    endif
    if (! ischar (varargin{2}) || ...
        ! any (strcmpi (varargin{2}, {'kendall', 'spearman'})))
      error ("copulastat: TYPE must be either 'Kendall' or 'Spearman'.");
    endif
    type = lower (varargin{2});
  endif

  lower_family = lower (family);

  switch (lower_family)

    case {'gaussian', 't'}
      ## Elliptical families: closed form, applied elementwise
      if (any (abs (param(:)) > 1))
        error ("copulastat: PARAM must be a correlation in the range -1 to 1.");
      endif
      if (strcmp (type, 'kendall'))
        r = 2 ./ pi .* asin (param);
      else
        r = 6 ./ pi .* asin (param ./ 2);
      endif

    case {'clayton', 'gumbel', 'frank'}
      ## Archimedean families: scalar copula parameter
      if (! isscalar (param))
        error ("copulastat: PARAM must be a scalar for the %s family.", family);
      endif
      switch (lower_family)
        case 'clayton'
          if (param < -1)
            error (strcat ("copulastat: PARAM must be greater than or", ...
                           " equal to -1 for the Clayton family."));
          endif
        case 'gumbel'
          if (param < 1)
            error (strcat ("copulastat: PARAM must be greater than or", ...
                           " equal to 1 for the Gumbel family."));
          endif
      endswitch
      if (strcmp (type, 'kendall'))
        r = archimedean_kendall (lower_family, param);
      else
        r = archimedean_spearman (lower_family, param);
      endif

    otherwise
      error ("copulastat: unknown copula family '%s'.", family);

  endswitch

endfunction

## Kendall's tau for the Archimedean families (closed form)
function tau = archimedean_kendall (family, alpha)
  switch (family)
    case 'clayton'
      tau = alpha ./ (alpha + 2);
    case 'gumbel'
      tau = 1 - 1 ./ alpha;
    case 'frank'
      if (alpha == 0)
        tau = 0;
      else
        tau = 1 - 4 ./ alpha .* (1 - debye1 (alpha));
      endif
  endswitch
endfunction

## Spearman's rho for the Archimedean families by numerical integration of the
## copula: rho_s = 12 * integral over the unit square of C(u,v), minus 3.
function rho = archimedean_spearman (family, alpha)
  if (family_is_independent (family, alpha))
    rho = 0;
    return;
  endif
  [x, w] = gauss_legendre_unit ();
  [U, V] = meshgrid (x, x);
  switch (family)
    case 'clayton'
      C = (U .^ (-alpha) + V .^ (-alpha) - 1) .^ (-1 ./ alpha);
    case 'gumbel'
      C = exp (-((-log (U)) .^ alpha + (-log (V)) .^ alpha) .^ (1 ./ alpha));
    case 'frank'
      C = -log (1 + (expm1 (-alpha .* U) .* expm1 (-alpha .* V)) ...
                    ./ expm1 (-alpha)) ./ alpha;
  endswitch
  rho = 12 .* (w(:)' * C * w(:)) - 3;
endfunction

## True when the parameter reduces the family to the independence copula
function tf = family_is_independent (family, alpha)
  tf = ((strcmp (family, 'clayton') || strcmp (family, 'frank')) ...
        && alpha == 0) || (strcmp (family, 'gumbel') && alpha == 1);
endfunction

## Debye function of the first order, D1(x) = (1/x) * integral_0^x t/(e^t-1) dt,
## valid for any real x (the removable singularity at t=0 is handled by expm1).
function d = debye1 (x)
  d = integral (@(t) t ./ expm1 (t), 0, x) ./ x;
endfunction

## Nodes and weights of a 64-point Gauss-Legendre rule mapped to [0,1], cached
## across calls (Golub-Welsch).
function [x, w] = gauss_legendre_unit ()
  persistent nodes weights
  if (isempty (nodes))
    N = 64;
    beta = 0.5 ./ sqrt (1 - (2 .* (1:N-1)) .^ (-2));
    T = diag (beta, 1) + diag (beta, -1);
    [V, D] = eig (T);
    [nodes, idx] = sort (diag (D));
    weights = 2 .* (V(1, idx) .^ 2)';
    nodes = (nodes + 1) ./ 2;
    weights = weights ./ 2;
  endif
  x = nodes;
  w = weights;
endfunction

%!demo
%! ## Kendall's tau and Spearman's rho of a Gaussian copula with correlation 0.5
%! tau = copulastat ("Gaussian", 0.5)
%! rho = copulastat ("Gaussian", 0.5, "type", "Spearman")

%!demo
%! ## Kendall's tau of a Clayton copula as its parameter grows
%! alpha = [0.5, 1, 2, 5];
%! tau = arrayfun (@(a) copulastat ("Clayton", a), alpha)

## Test output against MATLAB
%!test
%! assert_equal (copulastat ("Gaussian", 0.5), 1/3, 1e-14);
%! assert_equal (copulastat ("Gaussian", 0.5, "type", "Spearman"), ...
%!               0.482583739530997, 1e-14);
%! assert_equal (copulastat ("t", 0.5), 1/3, 1e-14);
%! assert_equal (copulastat ("Clayton", 2), 0.5, 1e-14);
%! assert_equal (copulastat ("Frank", 3), 0.307246959430723, 1e-12);
%! assert_equal (copulastat ("Gumbel", 2), 0.5, 1e-14);

## Elliptical families accept a correlation matrix and act elementwise
%!test
%! rho = [1, 0.5; 0.5, 1];
%! assert_equal (copulastat ("Gaussian", rho), 2/pi .* asin (rho), 1e-14);

## Independence limits
%!test
%! assert_equal (copulastat ("Clayton", 0), 0, 1e-14);
%! assert_equal (copulastat ("Frank", 0), 0, 1e-14);
%! assert_equal (copulastat ("Gumbel", 1), 0, 1e-14);
%! assert_equal (copulastat ("Clayton", 0, "type", "Spearman"), 0, 1e-14);
%! assert_equal (copulastat ("Gumbel", 1, "type", "Spearman"), 0, 1e-14);

## Test input validation
%!error <copulastat: FAMILY must be one of 'Gaussian', 't', 'Clayton', 'Gumbel', and 'Frank'.> ...
%! copulastat (5, 0.5)
%!error <copulastat: PARAM must be real.> copulastat ("Gaussian", 2i)
%!error <copulastat: TYPE must be either 'Kendall' or 'Spearman'.> ...
%! copulastat ("Gaussian", 0.5, "type", "Pearson")
%!error <copulastat: invalid optional argument.> ...
%! copulastat ("Gaussian", 0.5, "foo")
%!error <copulastat: PARAM must be a correlation in the range -1 to 1.> ...
%! copulastat ("Gaussian", 1.5)
%!error <copulastat: PARAM must be a scalar for the Clayton family.> ...
%! copulastat ("Clayton", [1, 2])
%!error <copulastat: PARAM must be greater than or equal to 1 for the Gumbel family.> ...
%! copulastat ("Gumbel", 0.5)
%!error <copulastat: unknown copula family 'Foo'.> copulastat ("Foo", 0.5)
