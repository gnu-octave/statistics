# Contribution Guidelines

#### 1. License and Documentation

The **statistics** package is distributed under [GNU General Public License (GPL)](https://www.gnu.org/licenses/gpl-3.0.en.html) with few minor exceptions. Some functions are in the Public Domain. If you are submitting a few function, it should be licensed under GPLv3+ with the following header (use appropriate year, name, etc.) as shown below:

```
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

```

New functions should be properly documented with the embedded help file in [Texinfo](https://www.gnu.org/software/texinfo/) format. This part should be placed outside (before) the function's main body and after the License block. The Texinfo manual can be found [here](https://www.gnu.org/software/texinfo/manual/texinfo/) although reading through existing function files should do the trick.

```
## -*- texinfo -*-
## @deftypefn {Function File} @var{p} = anova1 (@var{x})
##
## Help file goes here.
##
## @end deftypefn
```

**Note:** the texinfo is not printed on the command window screen as it appears in the source file. Type e.g. `help anova1` to display the help document in `anova1` function in order to review it.

#### 2. Demos

Although examples of using a function can go into the help documentation, it is always useful to have example embedded as demos, which the user can invoke with the `demo` command.

```
e.g. >> demo anova1
```

These should go at the end of the file after the function or any other helper functions that might be present in the file. A small demo sample is shown below.

```
%!demo
%! x = meshgrid (1:6);
%! x = x + normrnd (0, 1, 6, 6);
%! anova1 (x, [], 'off');
```

Note: demos can also contain comments which are printed when a demo is ran and may be very useful for the user to understand the function's usage. Embed comments in the usual manner after the `%!` starting characters in each line. For example:

```
%!demo
%! x = [46 64 83 105 123 150 150];
%! c = [0 0 0 0 0 0 1];
%! f = [1 1 1 1 1 1 4];
%! ## Subtract 30.92 from x to simulate a 3 parameter wbl with gamma = 30.92
%! wblplot (x - 30.92, c, f, 0.05);
```

#### 3. BISTs (testing suite)

It is also **very important** that function files contain a testing suite that will test for correct output, properly catching errors etc. BISTs should go at the bottom the file using `%!` starting characters. An example of such a testing suite is shown below.

```
%!error chi2gof ()
%!error chi2gof ([2,3;3,4])
%!error chi2gof ([1,2,3,4], "nbins", 3, "ctrs", [2,3,4])
%!error chi2gof ([1,2,3,4], "expected", [3,2,2], "nbins", 5)

%!test
%! x = [1 2 1 3 2 4 3 2 4 3 2 2];
%! [h, p, stats] = chi2gof (x);
%! assert (h, 0);
%! assert (p, NaN);
%! assert (stats.chi2stat, 0.1205375022748029, 1e-14);
```

It is best practice to append tests during development, since it will save time debugging and it will certainly help catching marginal errors that would otherwise have quickly come back as bug reports. It also **important** appending tests when fixing a bug or implementing some new functionality into an existing function.

Please, add BISTs, they help a lot with the maintenance :wink::innocent: of the statistics package. There is not such as a thing as too many tests in a function. :metal::v:

#### 4. Coding style

The coding style of GNU Octave should be used. In general, limit the lines at 80 characters long.
- Use `LF` (unix) for end of lines, and NOT `CRLF` (windows).
- Use `##` for comments. Don't use `%` or `%%` as in Matlab.
- Use `!` instead of `~` for logical NOT.
```
a != 0;
b(! isnnan (a)) = [];
```
- **Don't use tabs!** Indent the bodies of statement blocks with 2 spaces.
- When calling functions, put spaces after commas and before the calling parentheses, as shown below:
```
x = max (sin (y+3), 2);
```
- An exception are matrix or cell constructors.
```
>> a = [sin(x), cos(x)];
>> b = {sin(x), cos(x)};
```
**Note:** spaces in the above example would result in a parse error!

- For an indexing expression, do not put a space after the identifier (this differentiates indexing and function calls nicely).
- When indexing, the space after a comma is not necessary if index expressions are simple, but add a space for complex expressions as shown below:

```
A(:,i,j)
A([1:i-1; i+1:n]
```

- Always use a specific end-of-block statement (like endif, endswitch) rather than the generic end.
- Enclose the if, while, until, and switch conditions in parentheses.
```
if (isvector (a))
  s = sum (a);
endif
```

Do not do this, however, with the iteration counter portion of a for statement!
```
for i = 1:n
  b(i) = sum (a(:,1));
end
```

**That's about it all I suppose!** Just keep more or less consistent with what's in the other function files of the statistics packages and you should be fine :smile:

