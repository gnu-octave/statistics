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
## @deftypefn {Private Function} {[@var{part}, @var{args}] =} cvPartitionOf (@var{args}, @var{Y}, @var{n}, @var{classname})
##
## Take the cross-validation options out of an argument list and build the
## partition they ask for.
##
## @var{args} is the caller's @code{varargin}; what comes back is the same
## list with @qcode{'CrossVal'}, @qcode{'KFold'}, @qcode{'Holdout'},
## @qcode{'Leaveout'} and @qcode{'CVPartition'} removed, so the remainder
## can be handed to the learner that fits each fold.  Only one of the last
## four may be given.
##
## @var{Y} is the response, used to stratify a classification partition so
## that every fold holds the classes in the proportions the data does; pass
## an empty @var{Y} for a regression model, which is partitioned by count
## alone.  @var{n} is the number of observations.
##
## The default is ten folds, which is what @qcode{'CrossVal', 'on'} asks
## for and what a caller that named no option at all gets.
##
## @end deftypefn

function [part, args] = cvPartitionOf (args, Y, n, classname)

  numFolds = 10;
  Holdout = [];
  Leaveout = 'off';
  CVPartition = [];
  given = 0;
  keep = {};

  while (numel (args) > 0)
    if (numel (args) < 2)
      error (strcat ("%s: optional arguments must be given in Name-Value", ...
                     " pairs."), classname);
    endif
    switch (lower (args{1}))

      case 'crossval'
        val = args{2};
        if (! (ischar (val) && any (strcmpi (val, {'on', 'off'}))))
          error ("%s: 'CrossVal' must be either 'on' or 'off'.", classname);
        endif

      case 'kfold'
        numFolds = args{2};
        if (! (isnumeric (numFolds) && isscalar (numFolds)
               && isreal (numFolds) && numFolds == fix (numFolds)
               && numFolds > 1))
          error (strcat ("%s: 'KFold' must be an integer value greater", ...
                         " than 1."), classname);
        endif
        given++;

      case 'holdout'
        Holdout = args{2};
        if (! (isnumeric (Holdout) && isscalar (Holdout)
               && isreal (Holdout) && Holdout > 0 && Holdout < 1))
          error (strcat ("%s: 'Holdout' must be a numeric value between", ...
                         " 0 and 1."), classname);
        endif
        given++;

      case 'leaveout'
        Leaveout = args{2};
        if (! (ischar (Leaveout)
               && any (strcmpi (Leaveout, {'on', 'off'}))))
          error ("%s: 'Leaveout' must be either 'on' or 'off'.", classname);
        endif
        given++;

      case 'cvpartition'
        CVPartition = args{2};
        if (! isa (CVPartition, 'cvpartition'))
          error ("%s: 'CVPartition' must be a cvpartition object.", ...
                 classname);
        endif
        given++;

      otherwise
        keep = [keep, args(1:2)];

    endswitch
    args(1:2) = [];
  endwhile

  if (given > 1)
    error (strcat ("%s: you can use only one of 'KFold', 'Holdout',", ...
                   " 'Leaveout', or 'CVPartition' options."), classname);
  endif

  if (! isempty (CVPartition))
    part = CVPartition;
    if (part.NumObservations != n)
      error (strcat ("%s: 'CVPartition' must be built over the same", ...
                     " number of observations as the data."), classname);
    endif
  elseif (! isempty (Holdout))
    part = stratified (Y, n, 'Holdout', Holdout);
  elseif (strcmpi (Leaveout, 'on'))
    part = cvpartition (n, 'LeaveOut');
  else
    part = stratified (Y, n, 'KFold', numFolds);
  endif

  args = keep;

endfunction

## Partition by the response when there is one to stratify by, and by count
## alone otherwise.  Leave-one-out has nothing to stratify: every fold holds
## exactly one observation.
function part = stratified (Y, n, kind, value)

  if (isempty (Y))
    part = cvpartition (n, kind, value);
  else
    part = cvpartition (Y, kind, value);
  endif

endfunction
