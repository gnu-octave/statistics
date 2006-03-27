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
## @deftypefn {Function File} {@var{vpath} =} hmmviterbi (@var{sequence}, @var{transprob}, @var{outprob})
## Uses the Viterbi algorithm to find the Viterbi path for a Hidden Markov
## Model given a sequence of outputs. The model assumes that the generation
## begins with state @code{1} in step @code{0} but does not include step
## @code{0} in the generated states and sequence
##
## Arguments are
##
## @itemize
## @item
## @var{sequence} is the vector of length @var{len} with given outputs. The
## outputs must be integers ranging from @code{1} to @code{columns(outprob)}
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
## @var{vpath} is the vector of the same length as @var{sequence} with the
## estimated hidden states. The states are integers ranging from @code{1} to
## @code{columns(transprob)}
## @end itemize
##
## @end deftypefn

## Example:
## transprob = [0.8 0.2; 0.4 0.6];
## outprob = [0.2 0.4 0.4; 0.7 0.2 0.1];
## [sequence, states] = hmmgenerate(25, transprob, outprob)
## vpath = hmmviterbi(sequence, transprob, outprob)

## References:
## - Matlab 7.0 documentation (pdf)
## - http://en.wikipedia.org/wiki/Viterbi_algorithm

function vpath = hmmviterbi(sequence, transprob, outprob)

  # Check arguments
  if (nargin != 3)
    usage("vpath = hmmviterbi(sequence, transprob, outprob)");
  endif

  if (! isvector(sequence) && ! isempty(sequence))
    error("hmmgenerate: sequence must be a vector")
  endif

  if (! ismatrix(transprob))
    error("hmmgenerate: transprob must be a non-empty numeric matrix");
  endif
  if (! ismatrix(outprob))
    error("hmmgenerate: outprob must be a non-empty numeric matrix");
  endif

  len = length(sequence);
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

  # Store the path starting from i in spath(i, :)
  spath = ones(nstate, len+1);
  # Set the first state for each path
  spath(:, 1) = (1:nstate)';
  # Store the probability for path i in spathprob(i)
  spathprob = transprob(1, :);

  # Find the most likely paths for the given output sequence
  for i = 1:len
    # Calculate the new probabilities for the continuation with each state
    nextpathprob = ((spathprob' .* outprob(:, sequence(i))) * ones(1, nstate)) .* transprob;
    # Find the paths with the highest probabilities
    [spathprob, mindex] = max(nextpathprob);
    # Update spath and spathprob with the new paths
    spath = spath(mindex, :);
    spath(:, i+1) = (1:nstate)';  
  endfor

  # Set vpath to the most likely path
  # We do not want the last state because we do not have an output for it
  [m, mindex] = max(spathprob);
  vpath = spath(mindex, 1:len);

endfunction

%!test
%! sequence = [1 2 1 1 1 2 2 1 2 3 3 3 3 2 3 1 1 1 1 3 3 2 3 1 3];
%! transprob = [0.8 0.2; 0.4 0.6];
%! outprob = [0.2 0.4 0.4; 0.7 0.2 0.1];
%! vpath = hmmviterbi(sequence, transprob, outprob);
%! expected = [1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1];
%! assert(vpath, expected);
