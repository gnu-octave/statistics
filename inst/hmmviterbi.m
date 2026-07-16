## Copyright (C) 2006, 2007 Arno Onken <asnelt@asnelt.org>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{vpath} =} hmmviterbi (@var{sequence}, @var{transprob}, @var{outprob})
## @deftypefnx {statistics} {@var{vpath} =} hmmviterbi (@dots{}, @code{"symbols"}, @var{symbols})
## @deftypefnx {statistics} {@var{vpath} =} hmmviterbi (@dots{}, @code{"statenames"}, @var{statenames})
##
## Viterbi path of a hidden Markov model.
##
## Use the Viterbi algorithm to find the Viterbi path of a hidden Markov
## model given a sequence of outputs. The model assumes that the generation
## starts in state @code{1} at step @code{0} but does not include step
## @code{0} in the generated states and sequence.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{sequence} is the vector of length @var{len} of given outputs. The
## outputs must be integers ranging from @code{1} to
## @code{columns (outprob)}.
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
## @var{vpath} is the vector of the same length as @var{sequence} of the
## estimated hidden states. The states are integers ranging from @code{1} to
## @code{columns (transprob)}.
## @end itemize
##
## If @code{"symbols"} is specified, then @var{sequence} is expected to be a
## sequence of the elements of @var{symbols} instead of integers ranging
## from @code{1} to @code{columns (outprob)}. @var{symbols} can be a cell array.
##
## If @code{"statenames"} is specified, then the elements of
## @var{statenames} are used for the states in @var{vpath} instead of
## integers ranging from @code{1} to @code{columns (transprob)}.
## @var{statenames} can be a cell array.
##
## @subheading Examples
##
## @example
## @group
## transprob = [0.8, 0.2; 0.4, 0.6];
## outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
## [sequence, states] = hmmgenerate (25, transprob, outprob);
## vpath = hmmviterbi (sequence, transprob, outprob);
## @end group
##
## @group
## symbols = @{"A", "B", "C"@};
## statenames = @{"One", "Two"@};
## [sequence, states] = hmmgenerate (25, transprob, outprob, ...
##                      "symbols", symbols, "statenames", statenames);
## vpath = hmmviterbi (sequence, transprob, outprob, ...
##         "symbols", symbols, "statenames", statenames);
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

function vpath = hmmviterbi (sequence, transprob, outprob, varargin)

  # Check arguments
  if (nargin < 3 || mod (length (varargin), 2) != 0)
    print_usage ();
  endif

  if (! ismatrix (transprob))
    error ("hmmviterbi: transprob must be a non-empty numeric matrix");
  endif
  if (! ismatrix (outprob))
    error ("hmmviterbi: outprob must be a non-empty numeric matrix");
  endif

  len = length (sequence);
  # nstate is the number of states of the hidden Markov model
  nstate = rows (transprob);
  # noutput is the number of different outputs that the hidden Markov model
  # can generate
  noutput = columns (outprob);

  # Check whether transprob and outprob are feasible for a hidden Markov model
  if (columns (transprob) != nstate)
    error ("hmmviterbi: transprob must be a square matrix");
  endif
  if (rows (outprob) != nstate)
    error (strcat ("hmmviterbi: outprob must have the same number of", ...
                   " rows as transprob"));
  endif

  # Flag for symbols
  usesym = false;
  # Flag for statenames
  usesn = false;

  # Process varargin
  for i = 1:2:length (varargin)
    # There must be an identifier: 'symbols' or 'statenames'
    if (! ischar (varargin{i}))
      print_usage ();
    endif
    # Upper case is also fine
    lowerarg = lower (varargin{i});
    if (strcmp (lowerarg, 'symbols'))
      if (length (varargin{i + 1}) != noutput)
        error (strcat ("hmmviterbi: number of symbols does not match", ...
                       " number of possible outputs"));
      endif
      usesym = true;
      # Use the following argument as symbols
      symbols = varargin{i + 1};
    # The same for statenames
    elseif (strcmp (lowerarg, 'statenames'))
      if (length (varargin{i + 1}) != nstate)
        error (strcat ("hmmviterbi: number of statenames does not match", ...
                       " number of states"));
      endif
      usesn = true;
      # Use the following argument as statenames
      statenames = varargin{i + 1};
    else
      error (strcat ("hmmviterbi: expected 'symbols' or 'statenames'", ...
                     sprintf (" but found '%s'", varargin{i})));
    endif
  endfor

  # Transform sequence from symbols to integers if necessary
  if (usesym)
    # sequenceint is used to build the transformed sequence
    sequenceint = zeros (1, len);
    for i = 1:noutput
      # Search for symbols(i) in the sequence, isequal will have 1 at
      # corresponding indices; i is the right integer for that symbol
      isequal = ismember (sequence, symbols(i));
      # We do not want to change sequenceint if the symbol appears a second
      # time in symbols
      if (any ((sequenceint == 0) & (isequal == 1)))
        isequal *= i;
        sequenceint += isequal;
      endif
    endfor
    if (! all (sequenceint))
      index = max ((sequenceint == 0) .* (1:len));
      error (strcat ("hmmviterbi: sequence(", int2str (index), ...
                     ") not in symbols"));
    endif
    sequence = sequenceint;
  else
    if (! isvector (sequence) && ! isempty (sequence))
      error ("hmmviterbi: sequence must be a vector");
    endif
    if (! all (ismember (sequence, 1:noutput)))
      index = max ((ismember (sequence, 1:noutput) == 0) .* (1:len));
      error (strcat ("hmmviterbi: sequence(", int2str (index), ...
                     ") out of range"));
    endif
  endif

  # Each row in transprob and outprob should contain log probabilities
  # => scale so that the sum is 1 and convert to log space
  # - for transprob
  s = sum (transprob, 2);
  s(s == 0) = 1;
  transprob = log (transprob ./ s);
  # - for outprob
  s = sum (outprob, 2);
  s(s == 0) = 1;
  outprob = log (outprob ./ s);

  # Viterbi recursion.  The model starts in state 1 at step 0, so the first
  # observation is emitted after one transition from state 1.  delta(k) holds
  # the log probability of the most likely path that ends in state k at the
  # current step; back(:, i) stores, for each state, the predecessor state on
  # that best path, used for the traceback.
  if (len == 0)
    vpath = zeros (1, 0);
  else
    delta = transprob(1, :)' + outprob(:, sequence(1));
    back = zeros (nstate, len);
    for i = 2:len
      # best predecessor for each current state, then add the emission
      [delta, back(:, i)] = max (delta + transprob, [], 1);
      delta = delta' + outprob(:, sequence(i));
    endfor
    # Trace back from the most likely final state
    vpath = zeros (1, len);
    [~, vpath(len)] = max (delta);
    for i = len - 1:-1:1
      vpath(i) = back(vpath(i + 1), i + 1);
    endfor
  endif

  # Transform vpath into statenames if requested
  if (usesn)
    vpath = reshape (statenames(vpath), 1, len);
  endif

endfunction

%!demo
%! ## Most likely (Viterbi) state path for a two-state, three-symbol model.
%!
%! transprob = [0.8, 0.2; 0.4, 0.6];
%! outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
%! sequence = [1, 2, 1, 1, 3, 2, 3, 1, 2, 3];
%! vpath = hmmviterbi (sequence, transprob, outprob)

%!demo
%! ## The state path can also be reported using custom state names.
%!
%! transprob = [0.95, 0.05; 0.10, 0.90];
%! outprob = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6; 1/10, 1/10, 1/10, 1/10, 1/10, 1/2];
%! [sequence, states] = hmmgenerate (12, transprob, outprob, ...
%!                                   "statenames", {"fair", "loaded"});
%! vpath = hmmviterbi (sequence, transprob, outprob, ...
%!                     "statenames", {"fair", "loaded"})

%!test
%! sequence = [1, 2, 1, 1, 1, 2, 2, 1, 2, 3, 3, 3, ...
%!             3, 2, 3, 1, 1, 1, 1, 3, 3, 2, 3, 1, 3];
%! transprob = [0.8, 0.2; 0.4, 0.6];
%! outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
%! vpath = hmmviterbi (sequence, transprob, outprob);
%! expected = [1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, ...
%!             1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1];
%! assert_equal (vpath, expected);

%!test
%! sequence = {'A', 'B', 'A', 'A', 'A', 'B', 'B', 'A', 'B', 'C', 'C', 'C', ...
%!             'C', 'B', 'C', 'A', 'A', 'A', 'A', 'C', 'C', 'B', 'C', 'A', 'C'};
%! transprob = [0.8, 0.2; 0.4, 0.6];
%! outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
%! symbols = {'A', 'B', 'C'};
%! statenames = {'One', 'Two'};
%! vpath = hmmviterbi (sequence, transprob, outprob, 'symbols', symbols, ...
%!                                                   'statenames', statenames);
%! expected = {'One', 'One', 'Two', 'Two', 'Two', 'One', 'One', 'One', ...
%!             'One', 'One', 'One', 'One', 'One', 'One', 'One', 'Two', ...
%!             'Two', 'Two', 'Two', 'One', 'One', 'One', 'One', 'One', 'One'};
%! assert_equal (vpath, expected);

%!test
%! ## The returned path is the true maximum-probability path.  A former bug
%! ## scored a spurious transition out of the last state, biasing the final
%! ## states of the path; these cases guard against a regression.
%! transprob = [0.1854, 0.8146; 0.5948, 0.4052];
%! outprob = [0.5536, 0.4464; 0.2712, 0.7288];
%! assert_equal (hmmviterbi ([1, 2], transprob, outprob), [2, 2]);

%!test
%! transprob = [0.6056, 0.3944; 0.2036, 0.7964];
%! outprob = [0.6525, 0.3475; 0.2878, 0.7122];
%! assert_equal (hmmviterbi ([1, 1, 2, 1], transprob, outprob), [1, 1, 1, 1]);

%!test
%! ## An empty sequence yields an empty path
%! vpath = hmmviterbi ([], [0.8, 0.2; 0.4, 0.6], [0.5, 0.5; 0.3, 0.7]);
%! assert_equal (vpath, zeros (1, 0));
