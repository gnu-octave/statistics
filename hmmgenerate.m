## Copyright (C) 2006   Arno Onken   <whyly(at)gmx.net>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{sequence}, @var{states}] =} hmmgenerate (@var{len}, @var{transprob}, @var{outprob})
## Generates an output sequence and hidden states for a Hidden Markov Model.
## The model starts with state @code{1} at step @code{0} but will not
## include step @code{0} in the generated states and sequence
##
## Arguments are
##
## @itemize
## @item
## @var{len} is the number of steps to generate. @var{sequence} and
## @var{states} will have @var{len} entries each
## @item
## @var{transprob} is the matrix with the transition probabilities for the
## states. @code{transprob(i,j)} is the probability for a transition to
## state @code{j} given state @code{i}
## @item
## @var{outprob} is the matrix with the output probabilities.
## @code{outprob(i,j)} is the probability for generating output @code{j}
## given state @code{i}
## @end itemize
##
## Return values are
##
## @itemize
## @item
## @var{sequence} is a vector of length @var{len} with the generated
## outputs. The outputs are integers ranging from @code{1} to
## @code{columns(outprob)}
## @item
## @var{states} is a vector of length @var{len} with the generated hidden
## states. The states are integers ranging from @code{1} to
## @code{columns(transprob)}
## @end itemize
##
## @end deftypefn

## Example:
## transprob = [0.8 0.2; 0.4 0.6];
## outprob = [0.2 0.4 0.4; 0.7 0.2 0.1];
## [sequence, states] = hmmgenerate(25, transprob, outprob)

## References:
## - Matlab 7.0 documentation (pdf)
## - http://en.wikipedia.org/wiki/Hidden_Markov_Model

function [sequence, states] = hmmgenerate(len, transprob, outprob)

  # Check arguments
  if (nargin != 3)
    usage("[sequence, states] = hmmgenerate(len, transprob, outprob)");
  endif

  if (! isscalar(len) || len < 0)
    error("hmmgenerate: len must be a non-negative scalar value")
  endif

  if (! ismatrix(transprob))
    error("hmmgenerate: transprob must be a non-empty numeric matrix");
  endif
  if (! ismatrix(outprob))
    error("hmmgenerate: outprob must be a non-empty numeric matrix");
  endif

  # nstate is the number of states of the Hidden Markow Model
  nstate = rows(transprob);
  # noutput is the number of different outputs that the Hidden Markow Model
  # can generate
  noutput = columns(outprob);

  # Check whether transprob and outprob are feasible for a Hidden Markow Model
  if (columns(transprob) != nstate)
    error("hmmgenerate: transprob must be a square matrix");
  endif
  if (rows(outprob) != nstate)
    error("hmmgenerate: outprob must have the same number of rows as transprob");
  endif

  # Each row in transprob and outprob should contain probabilities
  # => scale so that the sum is 1
  # A zero row remains zero
  # - for transprob
  s = sum(transprob, 2);
  s(s==0) = 1;
  transprob = transprob ./ (s * ones(1, columns(transprob)));
  # - for outprob
  s = sum(outprob, 2);
  s(s==0) = 1;
  outsprob = outprob ./ (s * ones(1, columns(outprob))); 

  # Generate sequences of uniformly distributed random numbers between 0 and 1
  # - for the state transitions
  transdraw = rand(1, len);
  # - and the outputs
  outdraw = rand(1, len);

  # Generate the return vectors
  # They remain unchanged if the according probability row of transprob
  # and outprob contain, respectively, only zeros
  sequence = ones(1, len);
  states = ones(1, len);

  # Calculate cumulated probabilities backwards for easy comparison with the
  # generated random numbers
  # Cumulated probability in first column must always be 1
  # We might have a zero row
  # - for transprob
  transprob(:, end:-1:1) = cumsum(transprob(:, end:-1:1), 2);
  transprob(:, 1) = 1;
  # - for outprob
  outprob(:, end:-1:1) = cumsum(outprob(:, end:-1:1), 2);
  outprob(:, 1) = 1;

  # cstate is the current state
  # Start in state 1 but do not include it in the states vector
  cstate = 1;
  for i = 1:len
    # Compare the randon number i of transdraw to the cumulated probability
    # of the state transition and set the transition accordingly
    states(i) = sum(transdraw(i) <= transprob(cstate, :));
    cstate = states(i);
    # The same for the output of the state
    sequence(i) = sum(outdraw(i) <= outprob(cstate, :));
  endfor

endfunction

%!test
%! len = 25;
%! transprob = [0.8 0.2; 0.4 0.6];
%! outprob = [0.2 0.4 0.4; 0.7 0.2 0.1];
%! [sequence, states] = hmmgenerate(len, transprob, outprob);
%! assert(length(sequence), len);
%! assert(length(states), len);
%! assert(min(sequence) >= 1);
%! assert(max(sequence) <= columns(outprob));
%! assert(min(states) >= 1);
%! assert(max(states) <= rows(transprob));
