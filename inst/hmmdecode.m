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
## @deftypefn  {statistics} {@var{pstates} =} hmmdecode (@var{sequence}, @var{transprob}, @var{outprob})
## @deftypefnx {statistics} {[@var{pstates}, @var{logpseq}] =} hmmdecode (@dots{})
## @deftypefnx {statistics} {[@var{pstates}, @var{logpseq}, @var{fs}, @var{bs}, @var{s}] =} hmmdecode (@dots{})
## @deftypefnx {statistics} {[@dots{}] =} hmmdecode (@dots{}, @code{"symbols"}, @var{symbols})
##
## Posterior state probabilities of a hidden Markov model.
##
## Calculate the posterior state probabilities of the sequence @var{sequence}
## from a hidden Markov model.  The posterior state probabilities are the
## conditional probabilities of being in each state given the whole observed
## sequence.  The model assumes that the generation starts in state @code{1}
## at step @code{0} but does not include step @code{0} in the sequence.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{sequence} is a vector of length @var{len} of given outputs.  The
## outputs must be integers ranging from @code{1} to @code{columns (outprob)}.
##
## @item
## @var{transprob} is the matrix of transition probabilities of the states.
## @code{transprob(i, j)} is the probability of a transition to state
## @code{j} given state @code{i}.
##
## @item
## @var{outprob} is the matrix of output probabilities.
## @code{outprob(i, j)} is the probability of generating output @code{j}
## given state @code{i}.
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{pstates} is the matrix of posterior state probabilities.  It has one
## row for each state and one column for each element of @var{sequence}.
## @code{pstates(i, j)} is the conditional probability that the model is in
## state @code{i} when it generates the @code{j}-th output of @var{sequence},
## given that @var{sequence} is emitted.
##
## @item
## @var{logpseq} is the logarithm of the probability of the sequence
## @var{sequence}.
##
## @item
## @var{fs} and @var{bs} are the scaled forward and backward probabilities,
## respectively, and @var{s} is the vector of scale factors used to keep the
## computation numerically stable.
## @end itemize
##
## If @code{"symbols"} is specified, then @var{sequence} is expected to be a
## sequence of the elements of @var{symbols} instead of integers ranging from
## @code{1} to @code{columns (outprob)}.  @var{symbols} can be a cell array.
##
## @subheading Examples
##
## @example
## @group
## transprob = [0.8, 0.2; 0.4, 0.6];
## outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
## [sequence, states] = hmmgenerate (25, transprob, outprob);
## pstates = hmmdecode (sequence, transprob, outprob);
## @end group
##
## @group
## symbols = @{"A", "B", "C"@};
## [sequence, states] = hmmgenerate (25, transprob, outprob, ...
##                                   "symbols", symbols);
## pstates = hmmdecode (sequence, transprob, outprob, "symbols", symbols);
## @end group
## @end example
##
## @subheading References
##
## @enumerate
## @item
## Wendy L. Martinez and Angel R. Martinez. @cite{Computational Statistics
## Handbook with MATLAB}. Appendix E, pages 547-557, Chapman & Hall/CRC,
## 2001.
##
## @item
## Lawrence R. Rabiner. A Tutorial on Hidden Markov Models and Selected
## Applications in Speech Recognition. @cite{Proceedings of the IEEE},
## 77(2), pages 257-286, February 1989.
## @end enumerate
## @end deftypefn

function [pstates, logpseq, fs, bs, s] = hmmdecode (sequence, transprob, ...
                                                    outprob, varargin)

  # Check arguments
  if (nargin < 3 || mod (numel (varargin), 2) != 0)
    print_usage ();
  endif

  if (! ismatrix (transprob))
    error ("hmmdecode: transprob must be a non-empty numeric matrix.");
  endif
  if (! ismatrix (outprob))
    error ("hmmdecode: outprob must be a non-empty numeric matrix.");
  endif

  # nstate is the number of states of the hidden Markov model
  nstate = rows (transprob);
  # noutput is the number of different outputs that the hidden Markov model
  # can generate
  noutput = columns (outprob);

  # Check whether transprob and outprob are feasible for a hidden Markov model
  if (columns (transprob) != nstate)
    error ("hmmdecode: transprob must be a square matrix.");
  endif
  if (rows (outprob) != nstate)
    error (strcat ("hmmdecode: outprob must have the same number of", ...
                   " rows as transprob."));
  endif

  # Flag for symbols
  usesym = false;

  # Process varargin
  for i = 1:2:numel (varargin)
    # There must be an identifier: 'symbols'
    if (! ischar (varargin{i}))
      print_usage ();
    endif
    # Upper case is also fine
    lowerarg = lower (varargin{i});
    if (strcmp (lowerarg, 'symbols'))
      if (numel (varargin{i + 1}) != noutput)
        error (strcat ("hmmdecode: number of symbols does not match", ...
                       " number of possible outputs."));
      endif
      usesym = true;
      # Use the following argument as symbols
      symbols = varargin{i + 1};
    else
      error (strcat ("hmmdecode: expected 'symbols'", ...
                     sprintf (" but found '%s'.", varargin{i})));
    endif
  endfor

  len = numel (sequence);

  # Transform sequence from symbols to integers if necessary
  if (usesym)
    # sequenceint is used to build the transformed sequence
    sequenceint = zeros (1, len);
    for i = 1:noutput
      # Search for symbols(i) in the sequence; isequal will have 1 at
      # corresponding indices; i is the right integer for that symbol
      isequal = ismember (sequence, symbols(i));
      # We do not want to change sequenceint if the symbol appears a second
      # time in symbols
      if (any ((sequenceint == 0) & (isequal == 1)))
        isequal *= i;
        sequenceint += isequal;
      endif
    endfor
    if (! all (sequenceint) && len > 0)
      index = max ((sequenceint == 0) .* (1:len));
      error (strcat ("hmmdecode: sequence(", int2str (index), ...
                     ") not in symbols."));
    endif
    sequence = sequenceint;
  else
    if (! isvector (sequence) && ! isempty (sequence))
      error ("hmmdecode: sequence must be a vector.");
    endif
    if (! all (ismember (sequence, 1:noutput)))
      index = max ((ismember (sequence, 1:noutput) == 0) .* (1:len));
      error (strcat ("hmmdecode: sequence(", int2str (index), ...
                     ") out of range."));
    endif
  endif

  # Each row in transprob and outprob should contain probabilities
  # => scale so that the sum is 1.  A zero row remains zero.
  # - for transprob
  ts = sum (transprob, 2);
  ts(ts == 0) = 1;
  transprob = transprob ./ ts;
  # - for outprob
  os = sum (outprob, 2);
  os(os == 0) = 1;
  outprob = outprob ./ os;

  # Prepend a dummy output so that the forward and backward recursions have a
  # clean starting column representing the initial state 1 at step 0.  The
  # dummy column is stripped from PSTATES before returning.
  seq = [noutput + 1, sequence(:)'];
  L = len + 1;

  # Scaled forward probabilities.  Column 1 holds the initial distribution:
  # the model starts in state 1 with probability 1.
  fs = zeros (nstate, L);
  fs(1, 1) = 1;
  s = ones (1, L);
  for count = 2:L
    fs(:, count) = outprob(:, seq(count)) .* (transprob' * fs(:, count - 1));
    # The scale factor normalizes each forward column to sum to 1
    s(count) = sum (fs(:, count));
    fs(:, count) = fs(:, count) ./ s(count);
  endfor

  # Scaled backward probabilities using the same scale factors
  bs = ones (nstate, L);
  for count = L - 1:-1:1
    bs(:, count) = (transprob * (bs(:, count + 1) .* ...
                    outprob(:, seq(count + 1)))) ./ s(count + 1);
  endfor

  # The log probability of the sequence is the sum of the log scale factors
  logpseq = sum (log (s));

  # Posterior state probabilities; strip the dummy starting column
  pstates = fs .* bs;
  pstates(:, 1) = [];

endfunction

%!demo
%! ## Posterior probability of each state at every step of an observed sequence.
%!
%! transprob = [0.95, 0.05; 0.10, 0.90];
%! outprob = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6; 1/10, 1/10, 1/10, 1/10, 1/10, 1/2];
%! sequence = hmmgenerate (10, transprob, outprob);
%! [pstates, logpseq] = hmmdecode (sequence, transprob, outprob)

%!test
%! transprob = [0.8, 0.2; 0.4, 0.6];
%! outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
%! sequence = [1, 2, 1, 1, 1, 2, 2, 1, 2, 3];
%! pstates = hmmdecode (sequence, transprob, outprob);
%! assert_equal (size (pstates), [2, 10]);
%! assert_equal (all (abs (sum (pstates, 1) - 1) < 1e-10), true);
%! assert_equal (all (pstates(:) >= 0 & pstates(:) <= 1), true);

%!test
%! transprob = [0.8, 0.2; 0.4, 0.6];
%! outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
%! sequence = [1, 2, 1, 1, 1, 2, 2, 1, 2, 3];
%! [pstates, logpseq] = hmmdecode (sequence, transprob, outprob);
%! ## Independent brute-force forward algorithm for the log probability
%! nstate = 2;
%! alpha = transprob(1, :) .* outprob(:, sequence(1))';
%! for t = 2:numel (sequence)
%!   alpha = (alpha * transprob) .* outprob(:, sequence(t))';
%! endfor
%! assert_equal (logpseq, log (sum (alpha)), 1e-10);

%!test
%! ## Symbols form must match the integer form
%! transprob = [0.8, 0.2; 0.4, 0.6];
%! outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
%! sequence = [1, 2, 1, 1, 1, 2, 2, 1, 2, 3];
%! symbseq = {'A', 'B', 'A', 'A', 'A', 'B', 'B', 'A', 'B', 'C'};
%! p1 = hmmdecode (sequence, transprob, outprob);
%! p2 = hmmdecode (symbseq, transprob, outprob, 'symbols', {'A', 'B', 'C'});
%! assert_equal (p1, p2, 1e-12);

%!test
%! ## Scaled forward/backward reproduce the posterior gamma
%! transprob = [0.9, 0.1; 0.3, 0.7];
%! outprob = [0.5, 0.5; 0.1, 0.9];
%! sequence = [1, 2, 2, 1, 2];
%! [pstates, logpseq, fs, bs, s] = hmmdecode (sequence, transprob, outprob);
%! recovered = fs .* bs;
%! recovered(:, 1) = [];
%! assert_equal (recovered, pstates, 1e-12);
%! assert_equal (logpseq, sum (log (s)), 1e-12);

%!test
%! ## Empty sequence: no columns in the posterior, unit probability
%! transprob = [0.8, 0.2; 0.4, 0.6];
%! outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
%! [pstates, logpseq] = hmmdecode ([], transprob, outprob);
%! assert_equal (size (pstates), [2, 0]);
%! assert_equal (logpseq, 0, 1e-12);

%!error <Invalid call to hmmdecode> hmmdecode ([1, 2])
%!error <transprob must be a square matrix.> ...
%! hmmdecode ([1, 2], [0.8, 0.2; 0.4, 0.6; 0.1, 0.9], [0.5, 0.5; 0.1, 0.9])
%!error <outprob must have the same number of rows as transprob.> ...
%! hmmdecode ([1, 2], [0.8, 0.2; 0.4, 0.6], [0.5, 0.5])
%!error <sequence\(2\) out of range.> ...
%! hmmdecode ([1, 5], [0.8, 0.2; 0.4, 0.6], [0.2, 0.4, 0.4; 0.7, 0.2, 0.1])
%!error <number of symbols does not match number of possible outputs.> ...
%! hmmdecode ([1, 2], [0.8, 0.2; 0.4, 0.6], [0.2, 0.4, 0.4; 0.7, 0.2, 0.1], ...
%!            'symbols', {'A', 'B'})
