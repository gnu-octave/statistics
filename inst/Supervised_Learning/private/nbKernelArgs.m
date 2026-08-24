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
## @deftypefn {Private Function} {[@var{K}, @var{S}, @var{W}] =} nbKernelArgs (@var{Kernel}, @var{Support}, @var{Width}, @var{k}, @var{p}, @var{classname})
##
## Resolve the kernel arguments of a naive Bayes classifier into one kernel and
## one support per predictor, and a bandwidth per class and predictor.
##
## Each of @var{Kernel} and @var{Support} may be given once for every predictor
## or once per predictor.  @var{W} is empty when no bandwidth was given, which
## leaves each density to choose its own.
##
## @seealso{ClassificationNaiveBayes, fitcnb}
## @end deftypefn

function [K, S, W] = nbKernelArgs (Kernel, Support, Width, k, p, classname)

  known = {'box', 'epanechnikov', 'normal', 'triangle'};

  ## Kernel
  if (isempty (Kernel))
    K = repmat ({'normal'}, 1, p);
  elseif (ischar (Kernel) && isrow (Kernel))
    kname = lower (Kernel);
    K = repmat ({kname}, 1, p);
  elseif (iscellstr (Kernel) && numel (Kernel) == p)
    K = cellfun (@lower, Kernel(:)', 'UniformOutput', false);
  else
    error (strcat (classname, ": 'Kernel' must be a character vector or a", ...
                   " cell array of character vectors, one per predictor."));
  endif
  bad = ! ismember (K, known);
  if (any (bad))
    error (strcat (classname, ": unsupported kernel", ...
                   sprintf (" '%s'.", K{find (bad, 1)})));
  endif

  ## Support
  if (isempty (Support))
    S = repmat ({'unbounded'}, 1, p);
  elseif (ischar (Support) && isrow (Support))
    sname = lower (Support);
    S = repmat ({sname}, 1, p);
  elseif (isnumeric (Support) && numel (Support) == 2)
    sbnds = sort (Support(:)');
    S = repmat ({sbnds}, 1, p);
  elseif (iscell (Support) && numel (Support) == p)
    S = Support(:)';
  else
    error (strcat (classname, ": 'Support' must be 'unbounded',", ...
                   " 'positive', a two element numeric vector, or a cell", ...
                   " array holding one of those per predictor."));
  endif
  for j = 1:p
    if (ischar (S{j}) && ! any (strcmpi (S{j}, {'positive', 'unbounded'})))
      error (strcat (classname, ": a character vector 'Support' must be", ...
                     " 'unbounded' or 'positive'."));
    endif
  endfor

  ## Width
  if (isempty (Width))
    W = [];
    return;
  endif
  if (! (isnumeric (Width) && isreal (Width) && all (Width(:) > 0)))
    error (strcat (classname, ": 'Width' must be positive and real."));
  endif
  if (isscalar (Width))
    W = repmat (Width, k, p);
  elseif (isequal (size (Width), [k, p]))
    W = Width;
  elseif (isrow (Width) && numel (Width) == p)
    W = repmat (Width(:)', k, 1);
  elseif (iscolumn (Width) && numel (Width) == k)
    W = repmat (Width(:), 1, p);
  else
    error (strcat (classname, ": 'Width' must be a scalar, one value per", ...
                   " predictor, one per class, or a matrix of one per", ...
                   " class and predictor."));
  endif

endfunction
