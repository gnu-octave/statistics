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
## @deftypefn {Private Function} {[@var{map}, @var{pow}] =} encodednames_to_row (@var{col_names}, @var{pred_names}, @var{cat_info})
##
## Map each encoded design column back onto the predictor it came from.
##
## @var{col_names} lists the names of the columns to be mapped, @var{pred_names}
## the model's predictors, and @var{cat_info} the structure of categorical names
## and levels used to encode them.
##
## @var{map} is a row vector holding, for each encoded column, the index of its
## predictor in @var{pred_names}, or zero when no predictor claims it.  A
## categorical predictor's indicator columns are named @qcode{@var{name}_@var{level}},
## so a column is matched by the longest predictor name it starts with, which
## keeps @qcode{ab_x} with @qcode{ab} rather than with @qcode{a}.
##
## @var{pow} is a row vector holding the exponent each column name carries, so
## that a column named @qcode{u^2} maps onto predictor @qcode{u} with an
## exponent of 2.  It is 1 for a name carrying no exponent.
##
## This helper is shared by the @code{LinearModel} and
## @code{GeneralizedLinearModel} classes.
##
## @end deftypefn

function [m, pow] = encodednames_to_row (col_names, pred_names, cat_info)
  m   = zeros (1, numel (col_names));
  pow = ones (1, numel (col_names));
  for k = 1:numel (col_names)
    nm = col_names{k};

    ## A column named 'u^2' belongs to predictor 'u' and carries the exponent.
    tok = regexp (nm, '^(.*)\^(\d+)$', 'tokens');
    if (! isempty (tok))
      nm     = tok{1}{1};
      pow(k) = str2double (tok{1}{2});
    endif

    idx = find (strcmp (pred_names, nm), 1);
    if (isempty (idx))
      idx = 0;
      best_len = 0;
      for j = 1:numel (pred_names)
        pj = pred_names{j};
        if (numel (nm) > numel (pj) + 1 ...
            && strncmp (nm, [pj, '_'], numel (pj) + 1) ...
            && numel (pj) > best_len)
          idx = j;
          best_len = numel (pj);
        endif
      endfor
    endif
    m(k) = idx;
  endfor
endfunction
