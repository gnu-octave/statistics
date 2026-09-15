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
## @deftypefn {Private Function} {[@var{Zx}, @var{qk}, @var{nlev}, @var{levels}, @var{gidx}] =} expand_random (@var{Z}, @var{G}, @var{n})
##
## Build the expanded random-effects design of a mixed model.
##
## @var{Z} and @var{G} are cell arrays with one entry per grouping term:
## @var{Z}@{k@} is the @var{n}-by-q_k random-effects design and @var{G}@{k@} the
## @var{n}-by-1 grouping variable of term k.  @var{Zx} is the sparse
## @var{n}-by-N matrix, N = sum_k q_k*nlev_k, whose column block for term k and
## level l holds the rows of @var{Z}@{k@} that belong to level l and zeros
## elsewhere.  @var{qk} and @var{nlev} are the per-term column and level counts,
## @var{levels}@{k@} the sorted levels of term k and @var{gidx}@{k@} the level
## index of each observation.
##
## This helper is shared by @code{__lmefit__} and @code{__glmefit__}.
##
## @end deftypefn

function [Zx, qk, nlev, levels, gidx] = expand_random (Z, G, n)
  nt = numel (Z);
  N = 0;
  qk = zeros (1, nt);
  nlev = zeros (1, nt);
  levels = cell (1, nt);
  gidx = cell (1, nt);
  I = J = V = cell (1, nt);
  for k = 1:nt
    qk(k) = columns (Z{k});
    [lev, ~, gi] = unique (G{k}(:));
    nlev(k) = numel (lev);
    levels{k} = lev;
    gidx{k} = gi;
    I{k} = repmat ((1:n)', qk(k), 1);
    J{k} = N + repmat ((gi - 1) * qk(k), qk(k), 1) ...
           + kron ((1:qk(k))', ones (n, 1));
    V{k} = Z{k}(:);
    N += qk(k) * nlev(k);
  endfor
  Zx = sparse (vertcat (I{:}), vertcat (J{:}), vertcat (V{:}), n, N);
endfunction
