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
## @deftypefn {Private Function} {[@var{P}, @var{errmsg}] =} ensemblePartition (@var{args}, @var{Y}, @var{n}, @var{stratify})
##
## The partition a cross-validated ensemble is fitted on.
##
## @var{args} holds the cross-validation pairs of the call: one of
## @qcode{'CrossVal'}, @qcode{'KFold'}, @qcode{'Holdout'}, @qcode{'Leaveout'}
## and @qcode{'CVPartition'}, or none, which gives ten folds.  @var{Y} is the
## response and @var{n} the number of observations; @var{stratify} says
## whether the folds keep the class proportions of @var{Y}, as a classifier's
## do.  @qcode{'CrossVal'} set to @qcode{'off'} gives an empty @var{P}.
## @var{errmsg} is the body of the message the caller should raise, or empty.
##
## @end deftypefn

function [P, errmsg] = ensemblePartition (args, Y, n, stratify)

  P = [];
  errmsg = "";
  if (mod (numel (args), 2) != 0)
    errmsg = "name-value arguments must be in pairs.";
    return;
  endif
  if (numel (args) > 2)
    errmsg = strcat ("specify only one of 'CrossVal', 'KFold', 'Holdout',", ...
                     " 'Leaveout' and 'CVPartition'.");
    return;
  endif
  if (stratify)
    data = Y;
  else
    data = n;
  endif
  kfold = min (10, n);
  if (isempty (args))
    P = cvpartition (data, 'KFold', kfold);
    return;
  endif
  name = args{1};
  val = args{2};
  if (! ischar (name))
    errmsg = "invalid parameter name in optional pair arguments.";
    return;
  endif
  switch (tolower (name))
    case 'crossval'
      if (! (ischar (val) && any (strcmpi (val, {'on', 'off'}))))
        errmsg = "'CrossVal' must be 'on' or 'off'.";
      elseif (strcmpi (val, 'on'))
        P = cvpartition (data, 'KFold', kfold);
      endif
    case 'kfold'
      if (! (isnumeric (val) && isscalar (val) && isreal (val)
             && val > 1 && val == fix (val)))
        errmsg = "'KFold' must be an integer greater than 1.";
      else
        P = cvpartition (data, 'KFold', val);
      endif
    case 'holdout'
      if (! (isnumeric (val) && isscalar (val) && isreal (val)
             && val > 0 && val < 1))
        errmsg = "'Holdout' must be a number between 0 and 1.";
      else
        P = cvpartition (data, 'Holdout', val);
      endif
    case 'leaveout'
      if (! (ischar (val) && any (strcmpi (val, {'on', 'off'}))))
        errmsg = "'Leaveout' must be 'on' or 'off'.";
      elseif (strcmpi (val, 'on'))
        P = cvpartition (n, 'LeaveOut');
      endif
    case 'cvpartition'
      if (! isa (val, 'cvpartition'))
        errmsg = "'CVPartition' must be a 'cvpartition' object.";
      elseif (numel (test (val, 1)) != n)
        errmsg = strcat ("'CVPartition' must partition the observations", ...
                         " the ensemble was fitted on.");
      else
        P = val;
      endif
    otherwise
      errmsg = "invalid parameter name in optional pair arguments.";
  endswitch

endfunction
