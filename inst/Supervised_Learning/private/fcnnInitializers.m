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
## @deftypefn {Private Function} {@var{init} =} fcnnInitializers (@var{acts}, @var{nlayers}, @var{outact})
##
## The weight initializer each layer of a fully connected network is built
## with.
##
## @var{acts} names the hidden layers, either as one character vector applying
## to all @var{nlayers} of them or as a cellstring naming them one by one, and
## @var{outact} names the output layer.  @var{init} is a cellstring of
## @var{nlayers}+1 entries, the output layer last, each one @qcode{'he'} or
## @qcode{'glorot'}.
##
## The scheme is not a setting and cannot be chosen: @code{init_scale} in
## @file{src/fcnn.cpp} picks it from the layer's activation alone, the
## rectifiers taking the wider He range because they pass only half their
## input and the symmetric activations taking Glorot, which accounts for the
## backward pass as well.  This function restates that rule so the model can
## report what it was built with; the two must be changed together.
##
## @end deftypefn

function init = fcnnInitializers (acts, nlayers, outact)

  ## The rectifiers, by the names fcnntrain accepts for them.  Everything
  ## else the engine knows ('linear', 'none', 'sigmoid', 'tanh', 'softmax')
  ## takes Glorot.
  rectifiers = {'relu', 'lrelu', 'prelu', 'elu', 'gelu'};

  if (ischar (acts))
    names = repmat ({acts}, 1, nlayers);
  else
    names = cellstr (acts);
    names = names(:)';
  endif
  names = [names, {outact}];

  init = cell (1, numel (names));
  for i = 1:numel (names)
    if (any (strcmpi (rectifiers, names{i})))
      init{i} = 'he';
    else
      init{i} = 'glorot';
    endif
  endfor

endfunction
