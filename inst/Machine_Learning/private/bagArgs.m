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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {[@var{o}, @var{errmsg}] =} bagArgs (@var{args}, @var{T}, @var{N}, @var{allowed})
##
## Parse the Name-Value arguments of a bagged ensemble's prediction methods.
##
## @var{args} holds the pairs, @var{T} is the number of trees, @var{N} the
## number of observations, and @var{allowed} names the arguments the calling
## method takes, from @qcode{'Mode'}, @qcode{'Trees'}, @qcode{'TreeWeights'},
## @qcode{'UseInstanceForTree'} and @qcode{'Weights'}.
##
## @var{o} has the fields @qcode{trees}, the indices of the trees to use;
## @qcode{tw}, a weight for each of them; @qcode{use}, an @math{NxT} logical
## matrix over those trees saying which may answer for which observation;
## @qcode{w}, the observation weights or empty; and @qcode{mode}.
## @var{errmsg} is the body of the message the caller should raise, or empty.
##
## @end deftypefn

function [o, errmsg] = bagArgs (args, T, N, allowed)

  errmsg = "";
  o = struct ('trees', 1:T, 'tw', [], 'use', [], 'w', [], ...
              'mode', 'cumulative');
  if (mod (numel (args), 2) != 0)
    errmsg = "name-value arguments must be in pairs.";
    return;
  endif

  tw = [];
  for i = 1:2:numel (args)
    name = args{i};
    val = args{i+1};
    if (! (ischar (name) && any (strcmpi (name, allowed))))
      errmsg = "invalid parameter name in optional pair arguments.";
      return;
    endif
    switch (tolower (name))
      case 'trees'
        if (ischar (val) && strcmpi (val, 'all'))
          o.trees = 1:T;
        elseif (isnumeric (val) && isvector (val) && isreal (val)
                && all (val >= 1) && all (val <= T)
                && all (val == fix (val)))
          o.trees = double (val(:)');
        else
          errmsg = strcat ("'Trees' must be 'all' or a vector of", ...
                           " indices of trees in the ensemble.");
          return;
        endif
      case 'treeweights'
        tw = val;
      case 'useinstancefortree'
        if (! (islogical (val) && isequal (size (val), [N, T])))
          errmsg = strcat ("'UseInstanceForTree' must be a logical", ...
                           " matrix with one row per observation and", ...
                           " one column per tree.");
          return;
        endif
        o.use = val;
      case 'weights'
        if (! (isnumeric (val) && isvector (val) && isreal (val)
               && numel (val) == N && all (val >= 0) && sum (val) > 0))
          errmsg = strcat ("'Weights' must be a nonnegative numeric", ...
                           " vector with one element per observation,", ...
                           " not all zero.");
          return;
        endif
        o.w = double (val(:));
      case 'mode'
        if (! (ischar (val)
               && any (strcmpi (val, {'cumulative', 'individual', ...
                                      'ensemble'}))))
          errmsg = strcat ("'Mode' must be 'cumulative', 'individual'", ...
                           " or 'ensemble'.");
          return;
        endif
        o.mode = tolower (val);
    endswitch
  endfor

  if (isempty (tw))
    o.tw = ones (1, numel (o.trees));
  else
    if (strcmp (o.mode, 'individual'))
      errmsg = "'TreeWeights' cannot be used in 'individual' mode.";
      return;
    endif
    if (! (isnumeric (tw) && isvector (tw) && isreal (tw)
           && numel (tw) == numel (o.trees) && all (tw >= 0)
           && sum (tw) > 0))
      errmsg = strcat ("'TreeWeights' must be a nonnegative numeric", ...
                       " vector with one element per tree used, not all", ...
                       " zero.");
      return;
    endif
    o.tw = double (tw(:)');
  endif

  if (isempty (o.use))
    o.use = true (N, numel (o.trees));
  else
    o.use = o.use(:, o.trees);
  endif

endfunction
