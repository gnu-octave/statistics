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
##

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{D} =} nbDistNames (@var{val}, @var{p}, @var{classname})
##
## Resolve the @qcode{DistributionNames} argument of a naive Bayes classifier
## into one distribution name per predictor.
##
## @var{val} is what the caller was given, either empty for the default, a
## character vector naming one distribution for every predictor, or a cell
## array of character vectors naming one per predictor.  @var{p} is the number
## of predictors and @var{classname} prefixes any error raised.
##
## @var{D} is a cell array with one name per predictor, except for
## @qcode{'mn'}, which is returned as the character vector it was given as.
## The multinomial is a single distribution over the whole predictor vector
## rather than one per predictor, so it cannot be named for some predictors
## and not others, and there is no per-predictor list to return.
##
## @seealso{ClassificationNaiveBayes, fitcnb}
## @end deftypefn

function D = nbDistNames (val, p, classname)

  known = {'kernel', 'mvmn', 'normal'};

  if (isempty (val))
    D = repmat ({'normal'}, 1, p);
    return;
  endif

  if (ischar (val) && isrow (val))
    ## The name goes into a variable first: inside braces a space before a
    ## call's paren splits it into two elements and the call loses its argument.
    dname = lower (val);
    if (strcmp (dname, 'mn'))
      D = dname;
      return;
    endif
    D = repmat ({dname}, 1, p);
  elseif (iscellstr (val))
    if (any (strcmpi (val, 'mn')))
      error (strcat (classname, ": the 'mn' distribution applies to every", ...
                     " predictor at once and cannot be named for some", ...
                     " of them."));
    endif
    if (numel (val) != p)
      error (strcat (classname, ": 'DistributionNames' must name one", ...
                     " distribution per predictor."));
    endif
    D = cellfun (@lower, val(:)', 'UniformOutput', false);
  else
    error (strcat (classname, ": 'DistributionNames' must be a character", ...
                   " vector or a cell array of character vectors."));
  endif

  bad = ! ismember (D, known);
  if (any (bad))
    error (strcat (classname, ": unsupported distribution", ...
                   sprintf (" '%s'.", D{find (bad, 1)})));
  endif

endfunction
