## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {[@var{L}, @var{dL}] =} linearLoss (@var{f}, @var{y}, @var{lossfun}, @var{epsilon})
##
## Per-observation fitted loss of a linear model, and its derivative.
##
## @var{f} is the raw score @code{X * Beta + Bias}, @var{y} the response:
## @math{+1}/@math{-1} for the two classification losses and the observed
## value for the two regression losses.  @var{L} is the loss of each
## observation and @var{dL} its derivative with respect to @var{f}.
##
## The four losses are MATLAB's, and the constants are its own as measured on
## R2024a rather than the textbook ones: @qcode{'mse'} carries the factor of
## one half that makes the reported objective agree, while @qcode{'hinge'},
## @qcode{'logit'} and @qcode{'epsiloninsensitive'} carry none.  Note that
## this is the loss the @emph{fit} minimizes, which is not the loss the
## @code{loss} method reports by default.
##
## At the kink of a non-differentiable loss the derivative is taken to be the
## one-sided value that is zero, which is the subgradient of least norm.
##
## @end deftypefn

function [L, dL] = linearLoss (f, y, lossfun, epsilon)

  switch (lossfun)

    case 'hinge'
      m = 1 - y .* f;
      active = m > 0;
      L = m .* active;
      dL = -y .* active;

    case 'logit'
      ## log (1 + exp (-y * f)) evaluated so that neither tail overflows:
      ## log (1 + exp (u)) is max (u, 0) + log1p (exp (-abs (u))).
      u = -(y .* f);
      L = max (u, 0) + log1p (exp (-abs (u)));
      dL = -y ./ (1 + exp (-u));

    case 'mse'
      r = f - y;
      L = 0.5 * (r .^ 2);
      dL = r;

    case 'epsiloninsensitive'
      r = f - y;
      m = abs (r) - epsilon;
      active = m > 0;
      L = m .* active;
      dL = sign (r) .* active;

    otherwise
      error ("linearLoss: unknown loss function '%s'.", lossfun);

  endswitch

endfunction
