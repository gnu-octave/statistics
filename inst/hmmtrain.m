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
## @deftypefn  {statistics} {[@var{esttr}, @var{estout}] =} hmmtrain (@var{sequence}, @var{transguess}, @var{outguess})
## @deftypefnx {statistics} {[@dots{}] =} hmmtrain (@dots{}, @code{"algorithm"}, @var{algorithm})
## @deftypefnx {statistics} {[@dots{}] =} hmmtrain (@dots{}, @code{"symbols"}, @var{symbols})
## @deftypefnx {statistics} {[@dots{}] =} hmmtrain (@dots{}, @code{"tolerance"}, @var{tol})
## @deftypefnx {statistics} {[@dots{}] =} hmmtrain (@dots{}, @code{"maxiterations"}, @var{maxiter})
## @deftypefnx {statistics} {[@dots{}] =} hmmtrain (@dots{}, @code{"pseudotransitions"}, @var{pseudotransitions})
## @deftypefnx {statistics} {[@dots{}] =} hmmtrain (@dots{}, @code{"pseudoemissions"}, @var{pseudoemissions})
## @deftypefnx {statistics} {[@dots{}] =} hmmtrain (@dots{}, @code{"verbose"}, @var{vflag})
##
## Estimate the parameters of a hidden Markov model from emitted sequences.
##
## Given one or more observed output sequences and initial guesses for the
## transition and output probability matrices, @code{hmmtrain} finds maximum
## likelihood estimates of the two matrices using the Baum-Welch algorithm
## (the default) or Viterbi training.  The model assumes that the generation
## starts in state @code{1} at step @code{0} but does not include step
## @code{0} in the sequence.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{sequence} is a vector of a sequence of given outputs, or, for training
## from several sequences, a cell array of such vectors or a matrix whose rows
## are individual sequences.  The outputs must be integers ranging from
## @code{1} to @code{columns (outguess)}.
##
## @item
## @var{transguess} is the initial guess for the matrix of transition
## probabilities.  @code{transguess(i, j)} is the probability of a transition
## to state @code{j} given state @code{i}.
##
## @item
## @var{outguess} is the initial guess for the matrix of output probabilities.
## @code{outguess(i, j)} is the probability of generating output @code{j}
## given state @code{i}.
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{esttr} is the estimated matrix of transition probabilities.
##
## @item
## @var{estout} is the estimated matrix of output probabilities.
## @end itemize
##
## @subheading Name-Value pair arguments
##
## @itemize @bullet
## @item
## @code{"algorithm"} selects the training algorithm, either
## @code{"BaumWelch"} (default) or @code{"Viterbi"}.  @code{"BaumWelch"}
## performs the standard forward-backward re-estimation and is recommended for
## most uses.  @code{"Viterbi"} performs segmental (hard) re-estimation from
## the most likely state path of each sequence; it is faster but only
## approximates the maximum-likelihood estimate.
##
## @item
## @code{"symbols"} specifies the possible outputs.  If given, @var{sequence}
## is expected to hold the elements of @var{symbols} instead of integers.
## @var{symbols} can be a cell array.
##
## @item
## @code{"tolerance"} is the convergence tolerance (default @code{1e-6}).  The
## algorithm terminates when the change in the log-likelihood and in both
## estimated matrices falls below @var{tol}.
##
## @item
## @code{"maxiterations"} is the maximum number of iterations (default
## @code{500}).  A warning is issued if the algorithm has not converged within
## this many iterations.
##
## @item
## @code{"pseudotransitions"} and @code{"pseudoemissions"} supply pseudo-count
## matrices for Viterbi training, used to keep transitions or outputs that are
## very unlikely to occur from collapsing to zero probability.
##
## @item
## @code{"verbose"}, when true, prints the log-likelihood and the change in the
## estimates at each iteration.
## @end itemize
##
## @subheading Examples
##
## @example
## @group
## transprob = [0.8, 0.2; 0.4, 0.6];
## outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
## sequence = hmmgenerate (100, transprob, outprob);
## transguess = [0.6, 0.4; 0.5, 0.5];
## outguess = [0.3, 0.3, 0.4; 0.5, 0.3, 0.2];
## [esttr, estout] = hmmtrain (sequence, transguess, outguess);
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

function [esttr, estout] = hmmtrain (sequence, transguess, outguess, varargin)

  # Check arguments
  if (nargin < 3 || mod (numel (varargin), 2) != 0)
    print_usage ();
  endif

  if (! ismatrix (transguess))
    error ("hmmtrain: transguess must be a non-empty numeric matrix.");
  endif
  if (! ismatrix (outguess))
    error ("hmmtrain: outguess must be a non-empty numeric matrix.");
  endif

  # nstate is the number of states of the hidden Markov model
  nstate = rows (transguess);
  # noutput is the number of different outputs that the hidden Markov model
  # can generate
  noutput = columns (outguess);

  # Check whether transguess and outguess are feasible for a hidden Markov model
  if (columns (transguess) != nstate)
    error ("hmmtrain: transguess must be a square matrix.");
  endif
  if (rows (outguess) != nstate)
    error (strcat ("hmmtrain: outguess must have the same number of", ...
                   " rows as transguess."));
  endif

  # Defaults
  algorithm = 'baumwelch';
  usesym = false;
  tol = 1e-6;
  maxiter = 500;
  verbose = false;
  pseudotr = [];
  pseudoout = [];

  # Process varargin
  for i = 1:2:numel (varargin)
    if (! ischar (varargin{i}))
      print_usage ();
    endif
    lowerarg = lower (varargin{i});
    if (strcmp (lowerarg, 'algorithm'))
      if (! ischar (varargin{i + 1}))
        error ("hmmtrain: algorithm must be a string.");
      endif
      algorithm = lower (varargin{i + 1});
      if (! any (strcmp (algorithm, {'baumwelch', 'viterbi'})))
        error (strcat ("hmmtrain: expected 'BaumWelch' or 'Viterbi'", ...
                       sprintf (" but found '%s'.", varargin{i + 1})));
      endif
    elseif (strcmp (lowerarg, 'symbols'))
      if (numel (varargin{i + 1}) != noutput)
        error (strcat ("hmmtrain: number of symbols does not match", ...
                       " number of possible outputs."));
      endif
      usesym = true;
      symbols = varargin{i + 1};
    elseif (strcmp (lowerarg, 'tolerance'))
      tol = varargin{i + 1};
      if (! (isscalar (tol) && isnumeric (tol) && tol > 0))
        error ("hmmtrain: tolerance must be a positive scalar.");
      endif
    elseif (strcmp (lowerarg, 'maxiterations'))
      maxiter = varargin{i + 1};
      if (! (isscalar (maxiter) && maxiter >= 1 && round (maxiter) == maxiter))
        error (strcat ("hmmtrain: maxiterations must be a positive", ...
                       " scalar integer."));
      endif
    elseif (strcmp (lowerarg, 'pseudotransitions'))
      pseudotr = varargin{i + 1};
      if (! ismatrix (pseudotr) || rows (pseudotr) != columns (pseudotr))
        error ("hmmtrain: pseudotransitions must be a square matrix.");
      endif
      if (rows (pseudotr) != nstate)
        error (strcat ("hmmtrain: pseudotransitions must have the same", ...
                       " size as transguess."));
      endif
    elseif (strcmp (lowerarg, 'pseudoemissions'))
      pseudoout = varargin{i + 1};
      if (! ismatrix (pseudoout))
        error ("hmmtrain: pseudoemissions must be a numeric matrix.");
      endif
      if (rows (pseudoout) != nstate || columns (pseudoout) != noutput)
        error (strcat ("hmmtrain: pseudoemissions must have the same", ...
                       " size as outguess."));
      endif
    elseif (strcmp (lowerarg, 'verbose'))
      verbose = logical (varargin{i + 1});
    else
      error (strcat ("hmmtrain: unknown parameter name", ...
                     sprintf (" '%s'.", varargin{i})));
    endif
  endfor

  # Pseudo-count matrices default to zero (no smoothing)
  if (isempty (pseudotr))
    pseudotr = zeros (nstate, nstate);
  endif
  if (isempty (pseudoout))
    pseudoout = zeros (nstate, noutput);
  endif

  if (! usesym)
    symbols = {};
  endif

  # Collect the sequences into a cell array of integer row vectors.  A cell
  # array holds one sequence per element (symbol sequences must use this or a
  # single vector form); a numeric matrix with several rows is one sequence per
  # row; anything else is treated as a single sequence.
  if (iscell (sequence) && ! usesym)
    seqs = sequence;
  elseif (iscell (sequence) && usesym && ! isempty (sequence) ...
          && iscell (sequence{1}))
    seqs = sequence;
  elseif (isnumeric (sequence) && rows (sequence) > 1)
    seqs = num2cell (sequence, 2);
  else
    seqs = {sequence};
  endif
  nseq = numel (seqs);
  for j = 1:nseq
    seqs{j} = checkseq (seqs{j}, usesym, symbols, noutput);
  endfor

  # Normalize the initial guesses so each row sums to 1 (a zero row stays zero)
  guesstr = normalizerows (transguess);
  guessout = normalizerows (outguess);

  usebaum = strcmp (algorithm, 'baumwelch');
  converged = false;
  loglik = 1;

  for iter = 1:maxiter
    oldloglik = loglik;
    loglik = 0;
    oldguesstr = guesstr;
    oldguessout = guessout;

    if (usebaum)
      # Baum-Welch: accumulate expected transition and output counts
      TR = zeros (nstate, nstate);
      OUT = zeros (nstate, noutput);
      for j = 1:nseq
        seqj = seqs{j};
        len = numel (seqj);
        if (len == 0)
          continue;
        endif
        [fs, bs, sc] = fwdback (seqj, guesstr, guessout, noutput);
        loglik += sum (log (sc));
        # Expected transition counts (includes the forced transition out of
        # the initial state 1 via the padding column fs(:,1))
        xi = zeros (nstate, nstate);
        for i = 1:len
          xi += (fs(:, i) * (bs(:, i + 1) .* guessout(:, seqj(i)))') ...
                / sc(i + 1);
        endfor
        TR += guesstr .* xi;
        # Expected output counts: sum of posteriors at positions per symbol
        for l = 1:noutput
          pos = find (seqj == l);
          if (! isempty (pos))
            OUT(:, l) += sum (fs(:, pos + 1) .* bs(:, pos + 1), 2);
          endif
        endfor
      endfor
    else
      # Viterbi training: count the transitions and outputs along the best
      # path of each sequence.  Only observed consecutive-state transitions are
      # counted -- no transition out of the initial state -- matching
      # hmmestimate and MATLAB's first Viterbi iteration.
      TR = pseudotr;
      OUT = pseudoout;
      for j = 1:nseq
        seqj = seqs{j};
        len = numel (seqj);
        if (len == 0)
          continue;
        endif
        vpath = hmmviterbi (seqj, guesstr, guessout);
        # Path log-likelihood (used only for the convergence test); this does
        # include the forced initial transition out of state 1.
        lp = log (guesstr(1, vpath(1))) + log (guessout(vpath(1), seqj(1)));
        OUT(vpath(1), seqj(1)) += 1;
        for i = 2:len
          TR(vpath(i - 1), vpath(i)) += 1;
          OUT(vpath(i), seqj(i)) += 1;
          lp += log (guesstr(vpath(i - 1), vpath(i))) ...
                + log (guessout(vpath(i), seqj(i)));
        endfor
        loglik += lp;
      endfor
    endif

    # Normalize the accumulated counts into probability matrices
    guesstr = normalizerows (TR);
    guessout = normalizerows (OUT);

    # Relative changes used for convergence
    dll = abs (loglik - oldloglik) / (1 + abs (oldloglik));
    dtr = norm (guesstr - oldguesstr, Inf) / nstate;
    dout = norm (guessout - oldguessout, Inf) / noutput;

    if (verbose)
      printf (strcat ("hmmtrain: iteration %d, log-likelihood = %g,", ...
                      " rel. change = %g\n"), iter, loglik, dll);
    endif

    if (dll < tol && dtr < tol && dout < tol)
      converged = true;
      break;
    endif
  endfor

  if (! converged)
    warning (strcat ("hmmtrain: algorithm did not converge to within", ...
                     sprintf (" tolerance %g in %d iterations.", ...
                              tol, maxiter)));
  endif

  esttr = guesstr;
  estout = guessout;

endfunction

## Scale each row of M to sum to 1; a row summing to zero is left unchanged.
function M = normalizerows (M)
  s = sum (M, 2);
  s(s == 0) = 1;
  M = M ./ s;
endfunction

## Padded, scaled forward-backward recursion (see hmmdecode).  Returns FS and
## BS of size nstate-by-(len+1) with a leading column for the initial state 1,
## and the vector SC of scale factors (SC(1) = 1).
function [fs, bs, sc] = fwdback (seqj, transprob, outprob, noutput)
  nstate = rows (transprob);
  len = numel (seqj);
  seq = [noutput + 1, seqj(:)'];
  L = len + 1;
  fs = zeros (nstate, L);
  fs(1, 1) = 1;
  sc = ones (1, L);
  for count = 2:L
    fs(:, count) = outprob(:, seq(count)) .* (transprob' * fs(:, count - 1));
    sc(count) = sum (fs(:, count));
    fs(:, count) = fs(:, count) ./ sc(count);
  endfor
  bs = ones (nstate, L);
  for count = L - 1:-1:1
    bs(:, count) = (transprob * (bs(:, count + 1) .* ...
                    outprob(:, seq(count + 1)))) ./ sc(count + 1);
  endfor
endfunction

## Validate one sequence and, when USESYM, map its symbols to integers.
function seqint = checkseq (seqj, usesym, symbols, noutput)
  len = numel (seqj);
  if (usesym)
    seqint = zeros (1, len);
    for i = 1:noutput
      isequal = ismember (seqj, symbols(i));
      if (any ((seqint == 0) & (isequal == 1)))
        isequal *= i;
        seqint += isequal;
      endif
    endfor
    if (! all (seqint) && len > 0)
      index = max ((seqint == 0) .* (1:len));
      error (strcat ("hmmtrain: sequence(", int2str (index), ...
                     ") not in symbols."));
    endif
  else
    if (! isvector (seqj) && ! isempty (seqj))
      error ("hmmtrain: each sequence must be a vector.");
    endif
    if (! all (ismember (seqj, 1:noutput)))
      index = max ((ismember (seqj, 1:noutput) == 0) .* (1:len));
      error (strcat ("hmmtrain: sequence(", int2str (index), ...
                     ") out of range."));
    endif
    seqint = seqj(:)';
  endif
endfunction

%!demo
%! ## Re-estimate a model with Baum-Welch, starting from rough initial guesses.
%!
%! transprob = [0.95, 0.05; 0.10, 0.90];
%! outprob = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6; 1/10, 1/10, 1/10, 1/10, 1/10, 1/2];
%! sequence = hmmgenerate (1000, transprob, outprob);
%! transguess = [0.8, 0.2; 0.2, 0.8];
%! outguess = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6; 1/8, 1/8, 1/8, 1/8, 1/8, 3/8];
%! [esttr, estout] = hmmtrain (sequence, transguess, outguess)

%!test
%! ## Baum-Welch recovers matrices close to the generating model
%! transprob = [0.9, 0.1; 0.1, 0.9];
%! outprob = [0.9, 0.1; 0.1, 0.9];
%! rand ("seed", 42);
%! sequence = hmmgenerate (500, transprob, outprob);
%! transguess = [0.8, 0.2; 0.2, 0.8];
%! outguess = [0.7, 0.3; 0.3, 0.7];
%! [esttr, estout] = hmmtrain (sequence, transguess, outguess);
%! assert_equal (size (esttr), [2, 2]);
%! assert_equal (size (estout), [2, 2]);
%! assert_equal (all (abs (sum (esttr, 2) - 1) < 1e-10), true);
%! assert_equal (all (abs (sum (estout, 2) - 1) < 1e-10), true);

%!test
%! ## Baum-Welch never decreases the data log-likelihood
%! transprob = [0.8, 0.2; 0.4, 0.6];
%! outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
%! rand ("seed", 7);
%! sequence = hmmgenerate (200, transprob, outprob);
%! transguess = [0.5, 0.5; 0.5, 0.5];
%! outguess = [0.4, 0.3, 0.3; 0.2, 0.4, 0.4];
%! [~, ll0] = hmmdecode (sequence, transguess, outguess);
%! [esttr, estout] = hmmtrain (sequence, transguess, outguess);
%! [~, ll1] = hmmdecode (sequence, esttr, estout);
%! assert_equal (ll1 >= ll0 - 1e-8, true);

%!test
%! ## Multiple sequences supplied as a cell array
%! transprob = [0.8, 0.2; 0.4, 0.6];
%! outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
%! rand ("seed", 11);
%! s1 = hmmgenerate (60, transprob, outprob);
%! s2 = hmmgenerate (80, transprob, outprob);
%! transguess = [0.6, 0.4; 0.5, 0.5];
%! outguess = [0.3, 0.3, 0.4; 0.5, 0.3, 0.2];
%! [esttr, estout] = hmmtrain ({s1, s2}, transguess, outguess);
%! assert_equal (all (abs (sum (esttr, 2) - 1) < 1e-10), true);
%! assert_equal (all (abs (sum (estout, 2) - 1) < 1e-10), true);

%!test
%! ## Viterbi training runs and returns stochastic matrices.  Distinguishable
%! ## states keep both visited; pseudo-counts guard against empty rows.
%! transprob = [0.9, 0.1; 0.2, 0.8];
%! outprob = [0.8, 0.2; 0.2, 0.8];
%! rand ("seed", 3);
%! sequence = hmmgenerate (300, transprob, outprob);
%! transguess = [0.8, 0.2; 0.3, 0.7];
%! outguess = [0.7, 0.3; 0.3, 0.7];
%! [esttr, estout] = hmmtrain (sequence, transguess, outguess, ...
%!                             'algorithm', 'Viterbi', ...
%!                             'pseudotransitions', [1, 1; 1, 1], ...
%!                             'pseudoemissions', [1, 1; 1, 1]);
%! assert_equal (all (abs (sum (esttr, 2) - 1) < 1e-10), true);
%! assert_equal (all (abs (sum (estout, 2) - 1) < 1e-10), true);

%!test
%! ## Symbols form matches the integer form
%! transprob = [0.8, 0.2; 0.4, 0.6];
%! outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
%! intseq = [1, 2, 1, 1, 3, 2, 2, 1, 3, 3, 1, 2, 1, 1, 2];
%! symbseq = {'A', 'B', 'A', 'A', 'C', 'B', 'B', 'A', 'C', 'C', ...
%!            'A', 'B', 'A', 'A', 'B'};
%! transguess = [0.6, 0.4; 0.5, 0.5];
%! outguess = [0.3, 0.3, 0.4; 0.5, 0.3, 0.2];
%! [t1, o1] = hmmtrain (intseq, transguess, outguess, 'maxiterations', 5);
%! [t2, o2] = hmmtrain (symbseq, transguess, outguess, 'maxiterations', 5, ...
%!                      'symbols', {'A', 'B', 'C'});
%! assert_equal (t1, t2, 1e-12);
%! assert_equal (o1, o2, 1e-12);

%!warning <did not converge> ...
%! transprob = [0.8, 0.2; 0.4, 0.6];
%! outprob = [0.2, 0.4, 0.4; 0.7, 0.2, 0.1];
%! rand ("seed", 1);
%! sequence = hmmgenerate (50, transprob, outprob);
%! hmmtrain (sequence, [0.5, 0.5; 0.5, 0.5], [0.4, 0.3, 0.3; 0.2, 0.4, 0.4], ...
%!           'maxiterations', 1);

%!error <Invalid call to hmmtrain> hmmtrain ([1, 2], [0.8, 0.2; 0.4, 0.6])
%!error <transguess must be a square matrix.> ...
%! hmmtrain ([1, 2], [0.8, 0.2; 0.4, 0.6; 0.1, 0.9], [0.5, 0.5; 0.1, 0.9])
%!error <outguess must have the same number of rows as transguess.> ...
%! hmmtrain ([1, 2], [0.8, 0.2; 0.4, 0.6], [0.5, 0.5])
%!error <expected 'BaumWelch' or 'Viterbi'> ...
%! hmmtrain ([1, 2], [0.8, 0.2; 0.4, 0.6], [0.5, 0.5; 0.1, 0.9], ...
%!           'algorithm', 'nope')
%!error <sequence\(2\) out of range.> ...
%! hmmtrain ([1, 5], [0.8, 0.2; 0.4, 0.6], [0.2, 0.4, 0.4; 0.7, 0.2, 0.1])
