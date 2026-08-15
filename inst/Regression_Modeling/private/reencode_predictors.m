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
## @deftypefn {Private Function} {@var{X_enc} =} reencode_predictors (@var{X_raw}, @var{pred_names}, @var{cat_info}, @var{enc_names})
##
## Re-encode raw predictor data to match a fitted model's encoded columns.
##
## Given new raw predictor data @var{X_raw} (numeric, one column per predictor
## in @var{pred_names}), the categorical level information @var{cat_info}
## recorded at fit time (fields @qcode{names} and @qcode{levels}), and the
## target encoded column names @var{enc_names}, return the
## @math{n}-by-@code{numel (enc_names)} matrix @var{X_enc} whose columns
## reproduce, in order, the encoding used when the model was fitted.
##
## Each encoded name is matched as a plain predictor, a power term
## @qcode{name^k}, or a categorical indicator @qcode{name_level}, so the result
## aligns with the fitted design.  When @var{enc_names} is omitted it defaults
## to @var{pred_names} (plain pass-through).
##
## This helper is shared by the @code{LinearModel} and
## @code{GeneralizedLinearModel} classes.
##
## @end deftypefn

function X_enc = reencode_predictors (X_raw, pred_names, cat_info, enc_names)
  if (nargin < 4)
    enc_names = pred_names;
  endif
  n     = rows (X_raw);
  X_enc = zeros (n, numel (enc_names));
  for c = 1:numel (enc_names)
    name = enc_names{c};

    j = find (strcmp (pred_names, name), 1);
    if (! isempty (j))
      X_enc(:, c) = X_raw(:, j);
      continue;
    endif

    tok = regexp (name, '^(.+)\^(\d+)$', 'tokens');
    if (! isempty (tok))
      j = find (strcmp (pred_names, tok{1}{1}), 1);
      k = str2double (tok{1}{2});
      X_enc(:, c) = X_raw(:, j) .^ k;
      continue;
    endif

    found = false;
    for j = 1:numel (pred_names)
      ci = [];
      if (! isempty (cat_info.names))
        ci = find (strcmp (cat_info.names, pred_names{j}));
      endif
      if (isempty (ci))
        continue;
      endif
      levels_j = cat_info.levels{ci};
      for L = 2:numel (levels_j)
        if (strcmp (name, sprintf ("%s_%s", pred_names{j}, char (levels_j{L}))))
          X_enc(:, c) = double (X_raw(:, j) == L);
          found = true;
          break;
        endif
      endfor
      if (found)
        break;
      endif
    endfor
  endfor
endfunction
