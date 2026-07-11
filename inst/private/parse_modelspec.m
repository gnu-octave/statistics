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
## @deftypefn {Private Function} {[@var{terms}, @var{has_intercept}, @var{coef_names}, @var{errmsg}] =} parse_modelspec (@var{modelspec}, @var{pred_names}, @var{n_preds}, @var{intercept_nv})
##
## Turn a model specification into a terms matrix and coefficient names.
##
## @var{modelspec} is either a keyword (@qcode{'constant'}, @qcode{'linear'},
## @qcode{'interactions'}, @qcode{'purequadratic'}, @qcode{'quadratic'}, or
## @qcode{'full'}) or a numeric terms matrix; @var{pred_names} lists the
## @var{n_preds} encoded predictor names, and @var{intercept_nv} is a logical
## flag requesting an intercept term.
##
## @var{terms} is the resulting @math{m}-by-(@var{n_preds}+1) terms matrix,
## @var{has_intercept} indicates whether an intercept row is present, and
## @var{coef_names} are the term (coefficient) names.  On invalid input the
## outputs are empty and @var{errmsg} holds a diagnostic message body; the
## caller should emit it under its own function name.  @var{errmsg} is empty on
## success.
##
## This helper is shared by the @code{LinearModel} and
## @code{GeneralizedLinearModel} classes.
##
## @end deftypefn

function [terms, has_intercept, coef_names, errmsg] = parse_modelspec ( ...
    modelspec, pred_names, n_preds, intercept_nv)

  terms = [];  has_intercept = [];  coef_names = {};  errmsg = '';
  p = n_preds;

  if (isempty (modelspec) ...
      || (ischar (modelspec) && strcmpi (modelspec, 'linear')))
    terms = [zeros(1, p+1); [eye(p), zeros(p, 1)]];

  elseif (ischar (modelspec) && strcmpi (modelspec, 'constant'))
    terms = zeros (1, p+1);

  elseif (ischar (modelspec) && strcmpi (modelspec, 'interactions'))
    linear_part = [zeros(1, p+1); [eye(p), zeros(p, 1)]];
    inter_part  = zeros (0, p+1);
    for i = 1:p
      for j = i+1:p
        row = zeros (1, p+1); row(i) = 1; row(j) = 1;
        inter_part = [inter_part; row];
      endfor
    endfor
    terms = [linear_part; inter_part];

  elseif (ischar (modelspec) && strcmpi (modelspec, 'purequadratic'))
    linear_part = [zeros(1, p+1); [eye(p), zeros(p, 1)]];
    quad_part   = zeros (p, p+1);
    for j = 1:p
      quad_part(j, j) = 2;
    endfor
    terms = [linear_part; quad_part];

  elseif (ischar (modelspec) && strcmpi (modelspec, 'quadratic'))
    linear_part = [zeros(1, p+1); [eye(p), zeros(p, 1)]];
    quad_part   = zeros (p, p+1);
    for j = 1:p
      quad_part(j, j) = 2;
    endfor
    inter_part = zeros (0, p+1);
    for i = 1:p
      for j = i+1:p
        row = zeros (1, p+1); row(i) = 1; row(j) = 1;
        inter_part = [inter_part; row];
      endfor
    endfor
    terms = [linear_part; inter_part; quad_part];

  elseif (ischar (modelspec) && strcmpi (modelspec, 'full'))
    terms = zeros (1, p+1);
    for k = 1:p
      idx_mat = nchoosek (1:p, k);
      for j = 1:rows (idx_mat)
        row               = zeros (1, p+1);
        row(idx_mat(j,:)) = 1;
        terms             = [terms; row];
      endfor
    endfor

  elseif (isnumeric (modelspec))
    terms = double (modelspec);
    if (size (terms, 2) == p)
      terms = [terms, zeros(rows (terms), 1)];
    elseif (size (terms, 2) == p + 1)
      if (! all (terms(:, end) == 0))
        terms = [];
        errmsg = "Last column of terms matrix must be all zeros.";
        return;
      endif
    else
      terms = [];
      errmsg = sprintf ("Terms matrix must have %d or %d columns.", p, p+1);
      return;
    endif

  else
    errmsg = "Unknown model specification.";
    return;
  endif

  if (! intercept_nv)
    int_rows = all (terms(:, 1:end-1) == 0, 2);
    terms    = terms(! int_rows, :);
  endif

  has_intercept = any (all (terms(:, 1:end-1) == 0, 2));

  n_terms    = rows (terms);
  coef_names = cell (1, n_terms);
  for t = 1:n_terms
    term_row = terms(t, 1:end-1);
    if (all (term_row == 0))
      coef_names{t} = '(Intercept)';
    else
      parts_t = {};
      for j = 1:numel (term_row)
        if (term_row(j) != 0)
          if (term_row(j) == 1)
            parts_t{end+1} = pred_names{j};
          else
            parts_t{end+1} = sprintf ("%s^%d", pred_names{j}, term_row(j));
          endif
        endif
      endfor
      coef_names{t} = strjoin (parts_t, ':');
    endif
  endfor
endfunction
