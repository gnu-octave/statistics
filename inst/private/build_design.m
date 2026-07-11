## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2026 Avanish Salunke <avanishsalunke16@gmail.com>
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
## @deftypefn {Private Function} {@var{X_design} =} build_design (@var{terms}, @var{X_enc})
##
## Build a model design matrix from a terms matrix and an encoded predictor
## matrix.
##
## @var{X_enc} is the @math{n}-by-@math{p} encoded predictor matrix (categorical
## columns already expanded to indicator variables) and @var{terms} is an
## @math{m}-by-(@math{p}+1) matrix whose @math{t}-th row gives the exponents of
## each encoded predictor in the @math{t}-th model term (the trailing column,
## reserved for the response, is ignored).  A row of all zeros denotes the
## intercept.
##
## The returned @var{X_design} is @math{n}-by-@math{m}; column @math{t} is the
## element-wise product of the encoded predictors raised to their exponents in
## term @math{t}.
##
## This helper is shared by the @code{LinearModel} and
## @code{GeneralizedLinearModel} classes.
##
## @end deftypefn

function X_design = build_design (terms, X_enc)
  n_obs    = rows (X_enc);
  n_coef   = rows (terms);
  p_enc    = columns (X_enc);
  X_design = zeros (n_obs, n_coef);
  for t = 1:n_coef
    term_row = terms(t, 1:p_enc);
    col_t    = ones (n_obs, 1);
    for j = find (term_row != 0)
      col_t = col_t .* (X_enc(:, j) .^ term_row(j));
    endfor
    X_design(:, t) = col_t;
  endfor
endfunction
