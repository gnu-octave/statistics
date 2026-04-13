## Copyright (C) 2026 Jayant Chauhan <0001jayant@gmail.com>
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
## @deftypefn {Private Function} {@var{VariableInfo} =} __fitlm_build_varinfo__ (@
##   @var{var_names}, @var{cat_flags}, @var{raw_cols}, @var{in_model})
##
## Build the VariableInfo struct for a LinearModel.
##
## Confirmed from MATLAB Block 2-P4:
##   VariableInfo is a STRUCT (mirroring table layout) with fields:
##     Class         - cell of char: 'double', 'categorical', etc.
##     Range         - cell of cell: {[min max]} for numeric, {[sorted_uniq]} for categorical
##     InModel       - logical column vector
##     IsCategorical - logical column vector
##
## Row names are stored separately (= var_names).
##
## @var{in_model} is a logical row vector of length numel(var_names): true for
## each predictor that appears in the model (false for response).
## @end deftypefn

function VariableInfo = __fitlm_build_varinfo__ (var_names, cat_flags, raw_cols, in_model)

  n_vars = numel (var_names);

  Class         = cell  (n_vars, 1);
  Range         = cell  (n_vars, 1);
  InModel       = false (n_vars, 1);
  IsCategorical = false (n_vars, 1);

  for j = 1:n_vars
    col = raw_cols{j};

    ## IsCategorical: from cat_flags (set by CategoricalVars or type detection)
    IsCategorical(j) = logical (cat_flags(j));

    ## Class
    if (isa (col, "categorical"))
      Class{j} = "categorical";
    elseif (iscellstr (col) || isstring (col))
      Class{j} = "cell";
    elseif (isnumeric (col))
      Class{j} = "double";
    else
      Class{j} = class (col);
    endif

    ## Range
    if (IsCategorical(j))
      ## Sorted unique levels
      if (isa (col, "categorical"))
        levs = cellstr (categories (col));
      elseif (iscellstr (col))
        levs = unique (col);
      elseif (isstring (col))
        levs = unique (cellstr (col));
      else
        ## Numeric forced categorical — unique numeric values
        uvals = unique (col(! isnan (col)));
        levs = uvals(:)';
      endif
      Range{j} = {levs};
    else
      ## Numeric — [min max]
      if (isnumeric (col))
        valid = col(! isnan (col));
        if (isempty (valid))
          Range{j} = {[NaN NaN]};
        else
          Range{j} = {[min(valid), max(valid)]};
        endif
      else
        Range{j} = {[]};
      endif
    endif

    ## InModel
    InModel(j) = logical (in_model(j));
  endfor

  VariableInfo.Class         = Class;
  VariableInfo.Range         = Range;
  VariableInfo.InModel       = InModel;
  VariableInfo.IsCategorical = IsCategorical;

endfunction
