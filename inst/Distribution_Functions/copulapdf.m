## Copyright (C) 2008 Arno Onken <asnelt@asnelt.org>
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
## @deftypefn  {statistics} {@var{y} =} copulapdf (@var{family}, @var{x}, @var{theta})
## @deftypefnx {statistics} {@var{y} =} copulapdf ('t', @var{x}, @var{theta}, @var{df})
##
## Copula family probability density functions (PDF).
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{family} is the copula family name. Currently, @var{family} can
## be @code{'Gaussian'} for the Gaussian family, @code{'t'} for the
## Student's t family, @code{'Clayton'} for the Clayton family,
## @code{'Gumbel'} for the Gumbel-Hougaard family, @code{'Frank'} for the
## Frank family, @code{'AMH'} for the Ali-Mikhail-Haq family, or
## @code{'FGM'} for the Farlie-Gumbel-Morgenstern family.  The last two are
## Octave extensions that MATLAB does not provide.
##
## @item
## @var{x} is the support where each row corresponds to an observation.
##
## @item
## @var{theta} is the parameter of the copula. For the Gaussian and
## Student's t families it is the linear correlation matrix, and a scalar
## is expanded to a bivariate one. For the remaining families the elements
## of @var{theta} must be greater than or equal to @code{-1} for the
## Clayton family, greater than or equal to @code{1} for the
## Gumbel-Hougaard family, arbitrary for the Frank family, and greater
## than or equal to @code{-1} and lower than @code{1} for the
## Ali-Mikhail-Haq family. Moreover, @var{theta} must be non-negative
## for dimensions greater than @code{2}. @var{theta} must be a column
## vector with the same number of rows as @var{x} or be scalar.  The
## Farlie-Gumbel-Morgenstern family instead takes one parameter for every
## subset of the variables of order two or more, so @var{theta} is a row
## vector of length @code{2^d-d-1} or a matrix with one such row per
## observation; parameter sets violating the family's linear constraints
## give @code{NaN}.
##
## @item
## @var{df} is the degrees of freedom of the Student's t family, and is
## required by it. It must be a vector with the same number of rows as
## @var{x} or be scalar.
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{y} is the probability density of the copula at each row of
## @var{x} and corresponding parameter @var{theta}.
## @end itemize
##
## @subheading Examples
##
## @example
## @group
## x = [0.2:0.2:0.6; 0.2:0.2:0.6];
## theta = [1; 2];
## y = copulapdf ("Clayton", x, theta)
## @end group
##
## @group
## y = copulapdf ("Gumbel", x, 2)
## @end group
## @end example
##
## @subheading References
##
## @enumerate
## @item
## Roger B. Nelsen. @cite{An Introduction to Copulas}. Springer,
## New York, second edition, 2006.
## @end enumerate
##
## Input arguments must be @qcode{double} or @qcode{single}; integer, logical,
## and character arrays are rejected.  MATLAB accepts a character array and
## evaluates it at the character codes, which Octave deliberately does not,
## since a character array is an integer type and integers are refused too.
##
## @seealso{copulacdf, copularnd}
## @end deftypefn

function y = copulapdf (family, x, theta, df)

  ## Check arguments
  if (nargin != 3 && (nargin != 4 || ! strcmpi (family, 't')))
    print_usage ();
  endif

  if (! ischar (family))
    error (strcat ("copulapdf: family must be one of 'Gaussian',", ...
                   " 't', 'Clayton', 'Gumbel', 'Frank', and 'AMH'."));
  endif

  ## Check for X and THETA being double or single
  if (! (isfloat (x) && isfloat (theta)))
    error ("copulapdf: X and THETA must be double or single.");
  endif

  if (! isempty (x) && ! ismatrix (x))
    error ("copulapdf: X must be a numeric matrix.");
  endif

  [n, d] = size (x);

  lower_family = lower (family);

  ## The Farlie-Gumbel-Morgenstern family carries one parameter per subset of
  ## the variables of order two or more, as it does in copulacdf.
  if (strcmp (lower_family, 'fgm'))
    if (! ismatrix (theta) || size (theta, 2) != (2 .^ d - d - 1) || ...
        (size (theta, 1) != 1 && size (theta, 1) != n))
      error (strcat ("copulapdf: THETA must be a row vector of length", ...
                     " 2^d-d-1 or a matrix of size N x (2^d-d-1)."));
    endif
    if (n > 1 && size (theta, 1) == 1)
      theta = repmat (theta, n, 1);
    endif
  endif

  ## The two elliptical families take a correlation matrix rather than a
  ## one-parameter THETA, and are validated the way copulacdf validates them.
  is_elliptical = any (strcmp (lower_family, {'gaussian', 't'}));
  if (is_elliptical)
    if (d == 2 && isscalar (theta))
      ## Expand a scalar to a correlation matrix
      theta = [1, theta; theta, 1];
    endif
    if (any (size (theta) != [d, d]) || any (diag (theta) != 1) || ...
        any (any (theta != theta')) || min (eig (theta)) <= 0)
      error ("copulapdf: THETA must be a correlation matrix.");
    endif
    if (nargin == 4)
      if (! isscalar (df) && (! isvector (df) || length (df) != n))
        error (strcat ("copulapdf: DF must be a vector with the same", ...
                       " number of rows as X or be scalar."));
      endif
      df = df(:);
    endif
  elseif (! strcmp (lower_family, 'fgm') && (! isvector (theta) ...
          || (! isscalar (theta) && size (theta, 1) != n)))
    error (strcat ("copulapdf: THETA must be a column vector with the", ...
                   " same number of rows as X or be scalar."));
  endif

  if (n == 0)
    ## Input is empty
    y = zeros (0, 1);
  else
    if (n > 1 && isscalar (theta) && ! is_elliptical ...
        && ! strcmp (lower_family, 'fgm'))
      theta = repmat (theta, n, 1);
    endif

    ## Truncate input to unit hypercube
    x(x < 0) = 0;
    x(x > 1) = 1;

    ## Compute the density according to family
    lowerarg = lower_family;

    if (strcmp (lowerarg, 'gaussian'))
      ## The Gaussian family: the density of the correlated normal relative to
      ## the independent one, at the normal quantiles of X.
      z = norminv (x);
      y = exp (-0.5 * sum ((z * (inv (theta) - eye (d))) .* z, 2)) ...
          ./ sqrt (det (theta));
      ## No parameter bounds check
      k = [];
    elseif (strcmp (lowerarg, 't'))
      ## The Student's t family: the multivariate t density at the t quantiles
      ## of X, divided by the univariate ones it would factor into.
      if (nargin < 4)
        error ("copulapdf: DF is required for the 't' copula family.");
      endif
      z = tinv (x, df);
      y = mvtpdf (z, theta, df) ./ prod (tpdf (z, df), 2);
      ## No parameter bounds check
      k = [];
    elseif (strcmp (lowerarg, 'clayton'))
      ## The Clayton family
      log_cdf = -log (max (sum (x .^ (repmat (-theta, 1, d)), 2) ...
                - d + 1, 0)) ./ theta;
      y = prod (repmat (theta, 1, d) .* repmat (0:(d - 1), n, 1) + 1, 2) ...
                                     .* exp ((1 + theta .* d) .* log_cdf - ...
                                            (theta + 1) .* sum (log (x), 2));
      ## Product copula at columns where theta == 0
      k = find (theta == 0);
      if (any (k))
        y(k) = 1;
      endif
      ## Check theta
      if (d > 2)
        k = find (! (theta >= 0) | ! (theta < inf));
      else
        k = find (! (theta >= -1) | ! (theta < inf));
      endif
    elseif (strcmp (lowerarg, 'gumbel'))
      ## The Gumbel-Hougaard family
      g = sum ((-log (x)) .^ repmat (theta, 1, d), 2);
      c = exp (-g .^ (1 ./ theta));
      y = ((prod (-log (x), 2)) .^ (theta - 1)) ./ prod (x, 2) .* c .* ...
          (g .^ (2 ./ theta - 2) + (theta - 1) .* g .^ (1 ./ theta - 2));
      ## Check theta
      k = find (! (theta >= 1) | ! (theta < inf));
    elseif (strcmp (lowerarg, 'frank'))
      ## The Frank family
      if (d != 2)
        error ("copulapdf: Frank copula PDF implemented as bivariate only.");
      endif
      y = (theta .* exp (theta .* (1 + sum (x, 2))) .* (exp (theta) - 1)) ./ ...
          (exp (theta) - exp (theta + theta .* x(:, 1)) + ...
           exp (theta .* sum (x, 2)) - exp (theta + theta .* x(:, 2))) .^ 2;
      ## Product copula at columns where theta == 0
      k = find (theta == 0);
      if (any (k))
        y(k) = 1;
      endif
      ## Check theta
      k = find (! (theta > -inf) | ! (theta < inf));
    elseif (strcmp (lowerarg, 'amh'))
      ## The Ali-Mikhail-Haq family
      if (d != 2)
        error (strcat ("copulapdf: Ali-Mikhail-Haq copula PDF", ...
                       " implemented as bivariate only."));
      endif
      z = theta .* prod (x - 1, 2) - 1;
      y = (theta .* (1 - sum (x, 2) - prod (x, 2) - z) - 1) ./ (z .^ 3);
      ## Check theta
      k = find (! (theta >= -1) | ! (theta < 1));
    elseif (strcmp (lowerarg, 'fgm'))
      ## The Farlie-Gumbel-Morgenstern family.  Differentiating the
      ## distribution once in every variable turns each u_i (1 - u_i) factor
      ## into (1 - 2 u_i) and leaves the leading product at one.
      bcomb = logical (floor (mod (((0:(2 .^ d - 1))' * 2 .^ ...
                       ((1 - d):0)), 2)));
      ecomb = ones (size (bcomb));
      ecomb(bcomb) = -1;
      ## Summation over all combinations of order >= 2
      bcomb = bcomb(sum (bcomb, 2) >= 2, end:-1:1);
      ## Linear constraints matrix
      ac = zeros (size (ecomb, 1), size (bcomb, 1));
      ## Matrix to compute y
      ap = zeros (n, size (bcomb, 1));
      for i = 1:size (bcomb, 1)
        ac(:, i) = -prod (ecomb(:, bcomb(i, :)), 2);
        ap(:, i) = prod (1 - 2 * x(:, bcomb(i, :)), 2);
      endfor
      y = 1 + sum (ap .* theta, 2);
      ## Check linear constraints
      k = false (n, 1);
      for i = 1:n
        k(i) = any (ac * theta(i, :)' > 1);
      endfor
    else
      error ("copulapdf: unknown copula family '%s'.", family);
    endif

    if (any (k))
      y(k) = NaN;
    endif

  endif

endfunction

## Test output
%!test
%! x = [0.2:0.2:0.6; 0.2:0.2:0.6];
%! theta = [1; 2];
%! y = copulapdf ('Clayton', x, theta);
%! expected_p = [0.9872; 0.7295];
%! assert_equal (y, expected_p, 0.001);
%!test
%! x = [0.2:0.2:0.6; 0.2:0.2:0.6];
%! y = copulapdf ('Gumbel', x, 2);
%! expected_p = [0.9468; 0.9468];
%! assert_equal (y, expected_p, 0.001);
%!test
%! x = [0.2, 0.6; 0.2, 0.6];
%! theta = [1; 2];
%! y = copulapdf ('Frank', x, theta);
%! expected_p = [0.9378; 0.8678];
%! assert_equal (y, expected_p, 0.001);
%!test
%! x = [0.2, 0.6; 0.2, 0.6];
%! theta = [0.3; 0.7];
%! y = copulapdf ('AMH', x, theta);
%! expected_p = [0.9540; 0.8577];
%! assert_equal (y, expected_p, 0.001);

## Test input validation
%!error<copulapdf: X and THETA must be double or single.> copulapdf ('Clayton', int32 ([0, 0]), 2)
%!error<copulapdf: X and THETA must be double or single.> copulapdf ('Clayton', [true, true], 2)
%!error<copulapdf: X and THETA must be double or single.> copulapdf ('Clayton', 'ab', 2)

## The Gaussian and Student's t families, verified against MATLAB R2024a.
%!test
%! x = [0.1, 0.2; 0.3, 0.6; 0.5, 0.4; 0.7, 0.9; 0.45, 0.55];
%! y = copulapdf ('Gaussian', x, 0.5);
%! assert_equal (y, [1.60177371945198; 0.998741486235102; 1.14241401106385; ...
%!                   1.31299420633171; 1.13661012971028], 1e-12);
%! y = copulapdf ('Gaussian', x, -0.3);
%! assert_equal (y, [0.653989319149703; 1.07700187314878; 1.04496288532534; ...
%!                   0.76398020459843; 1.05211178116014], 1e-12);

%!test  # an uncorrelated Gaussian copula is the independence copula
%! x = [0.1, 0.2; 0.3, 0.6; 0.5, 0.4; 0.7, 0.9; 0.45, 0.55];
%! assert_equal (copulapdf ('Gaussian', x, 0), ones (5, 1), 1e-12);

%!test
%! x = [0.1, 0.2; 0.3, 0.6; 0.5, 0.4; 0.7, 0.9; 0.45, 0.55];
%! y = copulapdf ('t', x, 0.5, 5);
%! assert_equal (y, [1.66488234707408; 1.00205894407007; 1.24574194276063; ...
%!                   1.25091343011063; 1.24054754220011], 1e-12);
%! y = copulapdf ('t', x, -0.3, 10);
%! assert_equal (y, [0.644342340201519; 1.12031977908833; 1.09385787202575; ...
%!                   0.730314426588563; 1.10518920332876], 1e-12);

%!test  # both elliptical families take a correlation matrix beyond two columns
%! R = [1, 0.4, 0.2; 0.4, 1, 0.3; 0.2, 0.3, 1];
%! x = [0.2, 0.4, 0.6; 0.5, 0.5, 0.5];
%! assert_equal (copulapdf ('Gaussian', x, R), ...
%!               [1.1162786093147; 1.14859096884849], 1e-12);
%! assert_equal (copulapdf ('t', x, R, 8), ...
%!               [1.19326414606185; 1.37528246652052], 1e-12);

%!test  # the density integrates to one over the unit square
%! g = ((1:40)' - 0.5) / 40;
%! [A, B] = meshgrid (g, g);
%! x = [A(:), B(:)];
%! assert_equal (sum (copulapdf ('Gaussian', x, 0.5)) / 1600, 1, 1e-3);
%! assert_equal (sum (copulapdf ('t', x, 0.5, 6)) / 1600, 1, 1e-2);

%!error<copulapdf: THETA must be a correlation matrix.> ...
%! copulapdf ('Gaussian', [0.2, 0.4], [1, 2; 2, 1])
%!error<copulapdf: DF is required for the 't' copula family.> ...
%! copulapdf ('t', [0.2, 0.4], 0.5)
%!error<Invalid call to copulapdf> copulapdf ('Gaussian', [0.2, 0.4], 0.5, 5)

## The Farlie-Gumbel-Morgenstern family, an Octave extension.  Its density is
## checked against a finite difference of copulacdf, which is independent of it.
%!test
%! x = [0.35, 0.62];
%! h = 1e-5;
%! for theta = [-1, -0.4, 0.7, 1]
%!   fd = (copulacdf ('FGM', [x(1)+h, x(2)+h], theta) ...
%!       - copulacdf ('FGM', [x(1)+h, x(2)-h], theta) ...
%!       - copulacdf ('FGM', [x(1)-h, x(2)+h], theta) ...
%!       + copulacdf ('FGM', [x(1)-h, x(2)-h], theta)) / (4 * h * h);
%!   assert_equal (copulapdf ('FGM', x, theta), fd, 1e-6);
%! endfor

%!test  # the bivariate density in closed form
%! x = [0.35, 0.62; 0.2, 0.3];
%! theta = 0.7;
%! assert_equal (copulapdf ('FGM', x, theta), ...
%!               1 + theta * (1 - 2 * x(:,1)) .* (1 - 2 * x(:,2)), 1e-14);

%!test  # a zero parameter gives the independence copula
%! x = [0.35, 0.62; 0.2, 0.3; 0.9, 0.1];
%! assert_equal (copulapdf ('FGM', x, 0), ones (3, 1), 1e-14);

%!test  # one parameter per subset of order two or more beyond two variables
%! x = [0.3, 0.5, 0.7];
%! theta = [0.1, 0.1, 0.1, 0.1];
%! assert_equal (copulapdf ('FGM', x, theta), 0.984, 1e-12);

%!test  # a parameter set outside the family's linear constraints gives NaN
%! assert_equal (copulapdf ('FGM', [0.3, 0.5], 2), NaN);

%!test  # the density integrates to one over the unit square
%! g = ((1:60)' - 0.5) / 60;
%! [A, B] = meshgrid (g, g);
%! x = [A(:), B(:)];
%! assert_equal (sum (copulapdf ('FGM', x, 0.7)) / 3600, 1, 1e-10);

%!error<copulapdf: THETA must be a row vector of length 2\^d-d-1> ...
%! copulapdf ('FGM', [0.3, 0.5, 0.7], [0.1, 0.1])
