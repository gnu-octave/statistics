## Copyright (C) 2012  Arno Onken
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
## @deftypefn  {Function File} {@var{r} =} copularnd (@var{family}, @var{theta}, @var{n})
## @deftypefnx {Function File} {@var{r} =} copularnd (@var{family}, @var{theta}, @var{n}, @var{d})
## @deftypefnx {Function File} {@var{r} =} copularnd ('t', @var{theta}, @var{df}, @var{n})
##
## Random arrays from the copula family distributions.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{family} is the copula family name. Currently, @var{family} can be
## @code{'Gaussian'} for the Gaussian family, @code{'t'} for the Student's t
## family, @code{'Clayton'} for the Clayton family, @code{'Frank'} for the
## Frank family, @code{'Gumbel'} for the Gumbel-Hougaard family, @code{'AMH'}
## for the Ali-Mikhail-Haq family, or @code{'FGM'} for the
## Farlie-Gumbel-Morgenstern family.  The last two are Octave extensions that
## MATLAB does not provide.  Every family but Clayton is generated as
## bivariate only.
##
## @item
## @var{theta} is the parameter of the copula. For the Gaussian and Student's t
## copula, @var{theta} must be a correlation matrix. For bivariate copulas
## @var{theta} can also be a correlation coefficient. For the Clayton, Frank
## and Gumbel-Hougaard families, @var{theta} must be a vector with the same
## number of elements as samples to be generated or be scalar.  Values outside
## a family's range give @code{NaN} rows: at or above @code{1} for the
## Gumbel-Hougaard family, at or above @code{-1} for the bivariate Clayton
## family, and any finite value for the Frank family.  The Ali-Mikhail-Haq
## family takes @var{theta} in @code{[-1, 1)} and the
## Farlie-Gumbel-Morgenstern family in @code{[-1, 1]}.
##
## @item
## @var{df} is the degrees of freedom for the Student's t family. @var{df} must
## be a vector with the same number of elements as samples to be generated or
## be scalar.
##
## @item
## @var{n} is the number of rows of the matrix to be generated. @var{n} must be
## a non-negative integer and corresponds to the number of samples to be
## generated.
##
## @item
## @var{d} is the number of columns of the matrix to be generated. @var{d} must
## be a positive integer and corresponds to the dimension of the copula.
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{r} is a matrix of random samples from the copula with @var{n} samples
## of distribution dimension @var{d}.
## @end itemize
##
## @subheading Examples
##
## @example
## @group
## theta = 0.5;
## r = copularnd ("Gaussian", theta);
## @end group
##
## @group
## theta = 0.5;
## df = 2;
## r = copularnd ("t", theta, df);
## @end group
##
## @group
## theta = 0.5;
## n = 2;
## r = copularnd ("Clayton", theta, n);
## @end group
## @end example
##
## @subheading References
##
## @enumerate
## @item
## Roger B. Nelsen. @cite{An Introduction to Copulas}. Springer, New York,
## second edition, 2006.
## @end enumerate
## @end deftypefn

function r = copularnd (family, theta, df, n)

  ## Check arguments
  if (nargin < 2)
    print_usage ();
  endif

  if (! ischar (family))
    error (strcat ("copularnd: family must be one of 'Gaussian',", ...
                   " 't', 'Clayton', 'Frank', 'Gumbel', 'AMH', and 'FGM'."));
  endif

  lower_family = lower (family);

  ## Check family and copula parameters
  switch (lower_family)

    case {'gaussian'}
      ## Gaussian family
      if (isscalar (theta))
        ## Expand a scalar to a correlation matrix
        theta = [1, theta; theta, 1];
      endif
      if (! ismatrix (theta) || any (diag (theta) != 1) || ...
          any (any (theta != theta')) || min (eig (theta)) <= 0)
        error ("copularnd: THETA must be a correlation matrix.");
      endif
      if (nargin > 3)
        d = n;
        if (! isscalar (d) || d != size (theta, 1))
          error ("copularnd: D must correspond to dimension of theta.");
        endif
      else
        d = size (theta, 1);
      endif
      if (nargin < 3)
        n = 1;
      else
        n = df;
        if (! isscalar (n) || (n < 0) || round (n) != n)
          error ("copularnd: N must be a non-negative integer.");
        endif
      endif

    case {'t'}
      ## Student's t family
      if (nargin < 3)
        print_usage ();
      endif
      if (isscalar (theta))
        ## Expand a scalar to a correlation matrix
        theta = [1, theta; theta, 1];
      endif
      if (! ismatrix (theta) || any (diag (theta) != 1) || ...
          any (any (theta != theta')) || min (eig (theta)) <= 0)
        error ("copularnd: THETA must be a correlation matrix.");
      endif
      if (! isscalar (df) && (! isvector (df) || length (df) != n))
        error (strcat ("copularnd: DF must be a vector with the same", ...
                       " number of rows as r or be scalar."));
      endif
      df = df(:);
      if (nargin < 4)
        n = 1;
      else
        if (! isscalar (n) || (n < 0) || round (n) != n)
          error ("copularnd: N must be a non-negative integer.");
        endif
      endif

    case {'clayton', 'frank', 'gumbel', 'amh', 'fgm'}
      ## Archimedian one parameter family
      if (nargin < 4)
        ## Default is bivariate
        d = 2;
      else
        d = n;
        if (! isscalar (d) || (d < 2) || round (d) != d)
          error ("copularnd: D must be an integer greater than 1.");
        endif
      endif
      ## Only the Clayton family is available for more than two dimensions.
      if (d != 2 && ! strcmp (lower_family, 'clayton'))
        error (strcat ("copularnd: the '%s' copula is implemented as", ...
                       " bivariate only."), family);
      endif
      if (nargin < 3)
        ## Default is one sample
        n = 1;
      else
        n = df;
        if (! isscalar (n) || (n < 0) || round (n) != n)
          error ("copularnd: N must be a non-negative integer.");
        endif
      endif
      if (! isvector (theta) || (! isscalar (theta) && size (theta, 1) != n))
        error (strcat ("copularnd: THETA must be a column vector with", ...
                       " the number of rows equal to N or be scalar."));
      endif
      if (n > 1 && isscalar (theta))
        theta = repmat (theta, n, 1);
      endif

    otherwise
      error ("copularnd: unknown copula family '%s'.", family);

  endswitch

  if (n == 0)
    ## Input is empty
    r = zeros (0, d);
  else

    ## Draw random samples according to family
    switch (lower_family)

      case {'gaussian'}
        ## The Gaussian family
        r = normcdf (mvnrnd (zeros (1, d), theta, n), 0, 1);
        ## No parameter bounds check
        k = [];

      case {'t'}
        ## The Student's t family
        r = tcdf (mvtrnd (theta, df, n), df);
        ## No parameter bounds check
        k = [];

      case {'clayton'}
        ## The Clayton family
        u = rand (n, d);
        if (d == 2)
          r = zeros (n, 2);
          ## Conditional distribution method for the bivariate case which also
          ## works for theta < 0
          r(:, 1) = u(:, 1);
          r(:, 2) = (1 + u(:, 1) .^ (-theta) .* (u(:, 2) .^ ...
                    (-theta ./ (1 + theta)) - 1)) .^ (-1 ./ theta);
        else
          ## Apply the algorithm by Marshall and Olkin:
          ## Frailty distribution for Clayton copula is gamma
          y = randg (1 ./ theta, n, 1);
          r = (1 - log (u) ./ repmat (y, 1, d)) .^ (-1 ./ repmat (theta, 1, d));
        endif
        k = find (theta == 0);
        if (any (k))
          ## Product copula at columns k
          r(k, :) = u(k, :);
        endif
        ## Continue argument check
        if (d == 2)
          k = find (! (theta >= -1) | ! (theta < inf));
        else
          k = find (! (theta >= 0) | ! (theta < inf));
        endif

      case {'frank'}
        ## The Frank family, by inverting the conditional distribution.  It has
        ## a closed form here, which the Gumbel family below does not.
        u = rand (n, 2);
        e = exp (-theta);
        eu = exp (-theta .* u(:, 1));
        w = 1 + u(:, 2) .* (1 - e) ./ (u(:, 2) .* (eu - 1) - eu);
        u2 = -log (w) ./ theta;
        r = [u(:, 1), u2];
        ## Product copula at rows where theta == 0
        k = find (theta == 0);
        if (any (k))
          r(k, :) = u(k, :);
        endif
        ## Check bounds
        k = find (! (theta > -inf) | ! (theta < inf));

      case {'gumbel'}
        ## The Gumbel-Hougaard family, by the algorithm of Marshall and Olkin.
        ## Its frailty is positive stable with index 1 / THETA, drawn by
        ## Kanter's method; the conditional distribution cannot be inverted in
        ## closed form as Clayton's and Frank's can.
        a = 1 ./ theta;
        U = pi * rand (n, 1);
        W = -log (rand (n, 1));
        aU = sin ((1 - a) .* U) .* (sin (a .* U) .^ (a ./ (1 - a))) ...
             ./ (sin (U) .^ (1 ./ (1 - a)));
        V = (aU ./ W) .^ ((1 - a) ./ a);
        E = -log (rand (n, 2));
        r = exp (- (E ./ V) .^ a);
        ## The independence copula at theta == 1, where the frailty degenerates
        k = find (theta == 1);
        if (any (k))
          r(k, :) = exp (-E(k, :));
        endif
        ## Check bounds
        k = find (! (theta >= 1) | ! (theta < inf));

      case {'amh'}
        ## The Ali-Mikhail-Haq family.  Its conditional distribution is a
        ## quadratic in 1 - V, so it inverts in closed form.
        u = rand (n, 2);
        w = u(:, 2);
        b = 1 - u(:, 1);
        qa = theta - theta .^ 2 .* b .^ 2 .* w;
        qb = 2 .* theta .* b .* w - (1 + theta);
        qc = 1 - w;
        z = w;
        lin = abs (qa) <= eps;
        z(! lin) = (-qb(! lin) - sqrt (qb(! lin) .^ 2 ...
                    - 4 .* qa(! lin) .* qc(! lin))) ./ (2 .* qa(! lin));
        z(lin) = -qc(lin) ./ qb(lin);
        r = [u(:, 1), 1 - z];
        ## Check bounds
        k = find (! (theta >= -1) | ! (theta < 1));

      case {'fgm'}
        ## The Farlie-Gumbel-Morgenstern family, likewise a quadratic.
        u = rand (n, 2);
        w = u(:, 2);
        A = theta .* (1 - 2 .* u(:, 1));
        v = w;
        nz = A != 0;
        v(nz) = ((1 + A(nz)) - sqrt ((1 + A(nz)) .^ 2 ...
                 - 4 .* A(nz) .* w(nz))) ./ (2 .* A(nz));
        r = [u(:, 1), v];
        ## Check bounds
        k = find (! (theta >= -1) | ! (theta <= 1));

    endswitch

    ## Out of bounds parameters
    if (any (k))
      r(k, :) = NaN;
    endif

  endif

endfunction

## Test output
%!test
%! theta = 0.5;
%! r = copularnd ('Gaussian', theta);
%! assert_equal (size (r), [1, 2]);
%! assert_equal (all ((all ((r >= 0) & (r <= 1)))(:)), true);
%!test
%! theta = 0.5;
%! df = 2;
%! r = copularnd ('t', theta, df);
%! assert_equal (size (r), [1, 2]);
%! assert_equal (all ((all ((r >= 0) & (r <= 1)))(:)), true);
%!test
%! theta = 0.5;
%! r = copularnd ('Clayton', theta);
%! assert_equal (size (r), [1, 2]);
%! assert_equal (all ((all ((r >= 0) & (r <= 1)))(:)), true);
%!test
%! theta = 0.5;
%! n = 2;
%! r = copularnd ('Clayton', theta, n);
%! assert_equal (size (r), [n, 2]);
%! assert_equal (all ((all ((r >= 0) & (r <= 1)))(:)), true);
%!test
%! theta = [1; 2];
%! n = 2;
%! d = 3;
%! r = copularnd ('Clayton', theta, n, d);
%! assert_equal (size (r), [n, d]);
%! assert_equal (all ((all ((r >= 0) & (r <= 1)))(:)), true);

## The Frank and Gumbel-Hougaard families.  A generator cannot be checked
## against MATLAB value for value, so the sample's rank correlation is checked
## against the rank correlation the family is defined to have.
%!test
%! rand ("seed", 7);
%! for theta = [-5, -2, 2, 5, 10]
%!   r = copularnd ("Frank", theta, 4000);
%!   assert_equal (size (r), [4000, 2]);
%!   assert_equal (all (r(:) >= 0 & r(:) <= 1), true);
%!   rho = corr (tiedrank (r(:,1)), tiedrank (r(:,2)));
%!   assert_equal (rho, copulastat ("Frank", theta, "type", "Spearman"), 0.05);
%! endfor

%!test
%! rand ("seed", 11);
%! for theta = [1.5, 2, 3, 5]
%!   r = copularnd ("Gumbel", theta, 4000);
%!   assert_equal (size (r), [4000, 2]);
%!   assert_equal (all (r(:) >= 0 & r(:) <= 1), true);
%!   rho = corr (tiedrank (r(:,1)), tiedrank (r(:,2)));
%!   assert_equal (rho, copulastat ("Gumbel", theta, "type", "Spearman"), 0.05);
%! endfor

%!test  # each family degenerates to independence at its own boundary
%! rand ("seed", 13);
%! r = copularnd ("Frank", 0, 3000);
%! assert_equal (corr (tiedrank (r(:,1)), tiedrank (r(:,2))), 0, 0.06);
%! r = copularnd ("Gumbel", 1, 3000);
%! assert_equal (corr (tiedrank (r(:,1)), tiedrank (r(:,2))), 0, 0.06);

%!test  # a parameter outside the family's range gives NaN rows
%! assert_equal (copularnd ("Gumbel", 0.5, 2), NaN (2, 2));
%! assert_equal (copularnd ("Frank", Inf, 2), NaN (2, 2));

%!test  # the default is a single bivariate draw
%! rand ("seed", 3);
%! assert_equal (size (copularnd ("Frank", 3)), [1, 2]);
%! assert_equal (size (copularnd ("Gumbel", 2)), [1, 2]);
%! assert_equal (size (copularnd ("Gumbel", 2, 5)), [5, 2]);

%!error<copularnd: the 'Frank' copula is implemented as bivariate only.> ...
%! copularnd ("Frank", 3, 5, 3)
%!error<copularnd: the 'Gumbel' copula is implemented as bivariate only.> ...
%! copularnd ("Gumbel", 2, 5, 3)

## The Ali-Mikhail-Haq and Farlie-Gumbel-Morgenstern families, Octave
## extensions.  As with the other generators the sample's rank correlation is
## checked against the one the family is defined to have.
%!test
%! rand ("seed", 21);
%! for theta = [-0.9, -0.5, 0.5, 0.9]
%!   r = copularnd ("AMH", theta, 4000);
%!   assert_equal (size (r), [4000, 2]);
%!   assert_equal (all (r(:) >= 0 & r(:) <= 1), true);
%!   rho = corr (tiedrank (r(:,1)), tiedrank (r(:,2)));
%!   assert_equal (rho, copulastat ("AMH", theta, "type", "Spearman"), 0.05);
%! endfor

%!test
%! rand ("seed", 23);
%! for theta = [-1, -0.5, 0.5, 1]
%!   r = copularnd ("FGM", theta, 4000);
%!   assert_equal (size (r), [4000, 2]);
%!   assert_equal (all (r(:) >= 0 & r(:) <= 1), true);
%!   rho = corr (tiedrank (r(:,1)), tiedrank (r(:,2)));
%!   assert_equal (rho, copulastat ("FGM", theta, "type", "Spearman"), 0.05);
%! endfor

%!test  # both are the independence copula at a zero parameter
%! rand ("seed", 29);
%! r = copularnd ("AMH", 0, 4000);
%! assert_equal (corr (tiedrank (r(:,1)), tiedrank (r(:,2))), 0, 0.05);
%! r = copularnd ("FGM", 0, 4000);
%! assert_equal (corr (tiedrank (r(:,1)), tiedrank (r(:,2))), 0, 0.05);

%!test  # a parameter outside the family's range gives NaN rows
%! assert_equal (copularnd ("AMH", 1, 2), NaN (2, 2));
%! assert_equal (copularnd ("AMH", -1.5, 2), NaN (2, 2));
%! assert_equal (copularnd ("FGM", 1.5, 2), NaN (2, 2));

%!error<copularnd: the 'AMH' copula is implemented as bivariate only.> ...
%! copularnd ("AMH", 0.5, 5, 3)
%!error<copularnd: the 'FGM' copula is implemented as bivariate only.> ...
%! copularnd ("FGM", 0.5, 5, 3)
