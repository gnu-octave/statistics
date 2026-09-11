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
## @deftypefn {Private Function} {[@var{learners}, @var{errmsg}] =} ecocSelectModels (@var{learners}, @var{idx})
##
## Keep a subset of the regularization strengths of every binary learner.
##
## @var{learners} is the cell of binary learners of an error correcting code
## and @var{idx} names the strengths to keep, as indices or as a logical
## vector.  Only a learner fitted over several strengths has any to choose
## between, which today is the linear one alone.
##
## @var{errmsg} is empty when every learner was narrowed, and otherwise the
## body of a message the caller raises under its own name.  The learners are
## left untouched when it is not empty, so a code is never half narrowed.
##
## @seealso{CompactClassificationECOC, ClassificationECOC}
## @end deftypefn

function [learners, errmsg] = ecocSelectModels (learners, idx)

  errmsg = '';
  out = learners;

  for j = 1:numel (learners)
    if (! ismethod (learners{j}, 'selectModels'))
      errmsg = strcat ("the binary learners are", " '", ...
                       class (learners{j}), "' models, which are fitted", ...
                       " over one regularization strength and have none", ...
                       " to select between.");
      return;
    endif
    try
      out{j} = selectModels (learners{j}, idx);
    catch err
      errmsg = err.message;
      return;
    end_try_catch
  endfor

  learners = out;

endfunction
