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
## @deftypefn {Private Function} {@var{T} =} variable_level_terms (@var{terms_enc}, @var{enc2raw}, @var{var_idx}, @var{n_vars}, @var{col_pow})
##
## Fold a terms matrix over encoded design columns onto one over the model's
## variables.
##
## @var{terms_enc} is the terms matrix expressed over the encoded design
## columns, @var{enc2raw} maps each of those columns onto a predictor,
## @var{var_idx} maps each predictor onto a variable, and @var{n_vars} is the
## number of variables.  @var{col_pow} gives the exponent each column name
## carries and defaults to all ones; it is what distinguishes a column named
## @qcode{u^2} from one named @qcode{u}.
##
## @var{enc2raw} and @var{col_pow} are indexed by the @emph{columns of
## @var{terms_enc}}, so they must be built from the names of those columns and
## not from the coefficient names, which differ whenever a factor appears only
## inside an interaction or a power.
##
## @var{T} has one column per variable and one row per distinct term.  The
## several indicator columns of a categorical predictor collapse onto the single
## variable they came from, so the duplicate rows this produces are dropped.
##
## This helper is shared by the @code{LinearModel} and
## @code{GeneralizedLinearModel} classes.
##
## @end deftypefn

function T = variable_level_terms (terms_enc, enc2raw, var_idx, n_vars, col_pow)
  if (nargin < 5)
    col_pow = ones (1, columns (terms_enc));
  endif
  n_terms = rows (terms_enc);
  T = zeros (n_terms, n_vars);
  for t = 1:n_terms
    for j = 1:columns (terms_enc)
      if (terms_enc(t, j) != 0 && enc2raw(j) > 0 && var_idx(enc2raw(j)) > 0)
        c = var_idx(enc2raw(j));
        T(t, c) = max (T(t, c), terms_enc(t, j) * col_pow(j));
      endif
    endfor
  endfor
  keep = true (n_terms, 1);
  for t = 2:n_terms
    if (any (all (T(1:t-1, :) == T(t, :), 2)))
      keep(t) = false;
    endif
  endfor
  T = T(keep, :);
endfunction
