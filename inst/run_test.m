## Copyright (C) 1995-2017 Friedrich Leisch
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
## @deftypefn  {statistics} {[@var{pval}, @var{chisq}] =} run_test (@var{x})
##
## Perform a chi-square test with 6 degrees of freedom based on the upward
## runs in the columns of @var{x}.
##
## @code{run_test} can be used to decide whether @var{x} contains independent
## data.
##
## The p-value of the test is returned in @var{pval}.
##
## If no output argument is given, the p-value is displayed.
## @end deftypefn

function [pval, chisq] = run_test (x)

  if (nargin != 1)
    print_usage ();
  endif

  A = [4529.4,  9044.9,  13568,  18091,  22615,  27892;
       9044.4,  18097,   27139,  36187,  45234,  55789;
       13568,   27139,   40721,  54281,  67852,  83685;
       18091,   36187,   54281,  72414,  90470, 111580;
       22615,   45234,   67852,  90470, 113262, 139476;
       27892,   55789,   83685, 111580, 139476, 172860];

  b = [1/6; 5/24; 11/120; 19/720; 29/5040; 1/840];

  n = rows (x);
  r = run_count (x, 6) - n * b * ones (1, columns (x));

  chisq = diag (r' * A * r)' / n;
  pval  = chi2cdf (chisq, 6);

  if (nargout == 0)
    printf ("pval: %g\n", pval);
  endif

endfunction


%!xtest
%! x = [45.0000, -78.0000, 1.2550, 55.4000, -9.0000, 27.0000]
%! [pval, chisq] = run_test (x)
%! assert (pval, [1, 1, 1, 1, 1, 1]) 
%! assert (chisq, [1, 1, 1, 1, 1, 1])

## bounds test
%!test
%!  x1 = randn (10,6);
%!  [pval1, chisq1] = run_test(x1);
%!  assert(pval1 >= 0 && pval1 <= 1);
%!  assert(chisq1 >= 0);
%!test 
%!  x2 = randn (1000,6);
%!  [pval2, chisq2] = run_test(x2);
%!  assert(pval2 >= 0 && pval2 <= 1);
%!  assert(chisq2 >= 0);

## dependent data
%!test  
%!  x3 = [1 2 3 4 5 6; 2 4 6 8 10 12; 3 6 9 12 15 18; ...
%!        4 8 12 16 20 24; 5 10 15 20 25 30; 6 12 18 24 30 36];
%!  [pval3, chisq3] = run_test (x3);
%!  assert (pval3 >= 0 && pval3 <= 1);
%!  assert (chisq3 >= 0)

## vectorized inputs
%!test
%!  x4 = ones (600, 1);
%!  [pval, chisq] = run_test (x4);
%!  assert (pval >= 0 && pval <= 1);
%!  assert (chisq >= 0);
%!test
%!  x5 = [ones(100, 1), 2*ones(100, 1), 3*ones(100, 1), ...
%!        4*ones(100, 1), 5*ones(100, 1), 6*ones(100, 1)];
%!  [pval, chisq] = run_test (x5);
%!  assert (pval >= 0 && pval <= 1);
%!  assert (chisq >= 0);
%!test
%!  x6 = randi([1, 6], 600, 1);  
%!  [pval, chisq] = run_test (x6);
%!  assert (pval >= 0 && pval <= 1);
%!  assert (chisq >= 0);

## inconsistent inputs
%!assert (run_test ([true, false, true, true, false, false]));
%!error <run_count: X must be a numeric vector or matrix> ...
%! run_test ('abcdef');
%!assert (run_test ([1,2;3,4;5,6]));
%!error <Invalid call to run_test.> ...
%! run_test ([1,2,3;4,5,6;7,8,9],[10,11,12;13,14,15;16,17,18]);
%!assert (run_test(42));

## input validation tests
%!error id=Octave:invalid-fun-call run_test();
%!error id=Octave:nonconformant-args run_test ([1, 2, 3.5, 4, 5]);
%!error id=Octave:nonconformant-args run_test ([1, 2, 3, 4, 5, 6, 7]);

