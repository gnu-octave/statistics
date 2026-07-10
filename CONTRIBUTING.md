# Contribution Guidelines

Thank you for considering a contribution to the **statistics** package!  These
guidelines describe how to write code that blends into the package so your
submission can be reviewed and merged quickly.  When in doubt, copy the
conventions of an existing neighbouring file rather than introducing new ones.

## 1. License

Every source file in the **statistics** package is licensed under the
[GNU General Public License, version 3 or later](https://www.gnu.org/licenses/gpl-3.0.en.html)
(GPLv3+) — there are no other licenses in the package.  If you are submitting a
new function, it must carry the following header (use the appropriate year and
your own name and email):

```
## Copyright (C) 2026 Your Name <your.email@example.com>
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
```

Use the current year for a brand-new file (e.g. `2026`); when editing an
existing file, extend its year into a range (e.g. `2022-2026`).


## 2. Documentation

New functions should be properly documented with the embedded help text in
[Texinfo](https://www.gnu.org/software/texinfo/) format.  The documentation
block goes **after** the License header and **before** the function body.  The
Texinfo manual is
[here](https://www.gnu.org/software/texinfo/manual/texinfo/), but reading
through existing function files should do the trick.

```
## -*- texinfo -*-
## @deftypefn  {statistics} {@var{p} =} anova1 (@var{x})
## @deftypefnx {statistics} {@var{p} =} anova1 (@var{x}, @var{group})
##
## Help text goes here.
##
## @seealso{anova2, anovan}
## @end deftypefn
```

A few conventions for the Texinfo block:

- Open every block with `## -*- texinfo -*-`.
- Standalone functions use `@deftypefn {statistics} ... @end deftypefn`;
  additional signatures use `@deftypefnx`.
- Mark arguments and variables with `@var{X}`, inline code with `@code{...}`,
  and quoted literals or option names with `@qcode{'name'}`.
- Add cross-references with `@seealso{a, b, c}` as the last line before
  `@end deftypefn`.
- Put **two spaces after a sentence-ending period** inside the block, matching
  GNU Texinfo style.

For **classdef classes**, the documentation is split across several blocks that
use `@deftp` (define-type) rather than `@deftypefn`:

- The **class** itself is documented with a
  `@deftp {statistics} <name> ... @end deftp` block, placed right after the
  `classdef` line.
- Each public **property** carries its own
  `@deftp {<class>} {property} <Name> ... @end deftp` block, immediately above
  the property in its `properties` section.
- Each public **method**, by contrast, is documented exactly like a standalone
  function — with a `@deftypefn {<class>} ... @end deftypefn` block above the
  method definition.  The scope name in the braces is the class name, and the
  method's error messages are prefixed `<class>.<method>:` (see §4).

In short: **classes and properties use `@deftp`; methods use `@deftypefn`**,
just as free functions do.
- Keep every body line **within 80 columns**.  The only lines allowed to run
  over are the `@deftypefn` / `@deftypefnx` header lines, whose signatures may
  be longer.

**Note:** the Texinfo source is not printed verbatim on the command window.
Type e.g. `help anova1` to render and review the documentation as the user
sees it.

## 3. Demos

Although examples can live in the help documentation, it is always useful to
embed runnable examples as demos, which the user invokes with the `demo`
command:

```
>> demo anova1
```

Demos go at the **end of the file**, after the function and any local helper
functions.  A small demo looks like this:

```
%!demo
%! x = meshgrid (1:6);
%! x = x + normrnd (0, 1, 6, 6);
%! anova1 (x, [], 'off');
```

### Demo comments and online documentation

The online reference pages for the package are generated with
[pkg-octave-doc](https://github.com/gnu-octave/pkg-octave-doc), which renders
each `%!demo` block as an interleaved **notebook**: the block is split into
cells so that every statement's console output and any figures appear
immediately beneath the code that produced them.

To make a demo read well both in the terminal (`demo funcname`) and online, its
comment lines support a **small subset of Markdown** instead of Texinfo.  The
supported constructs are:

- inline `` `code` ``, `**bold**`, `*italic*`, and `[text](url)` links;
- unordered lists (lines starting with `-` or `*`) and ordered lists (lines
  starting with `1.`);
- paragraphs, separated by a blank comment line.

Underscore emphasis (`_text_`) and `#` headings are intentionally **not**
supported — they would clash with identifier names and the `#` comment marker.
Do **not** use Texinfo tags inside demo comments.

```
%!demo
%! ## Subtract 30.92 from x to simulate a 3-parameter Weibull with
%! ## `gamma = 30.92`, then plot the result.
%! x = [46 64 83 105 123 150 150];
%! c = [0 0 0 0 0 0 1];
%! f = [1 1 1 1 1 1 4];
%! wblplot (x - 30.92, c, f, 0.05);
```

Consecutive statements that print nothing are merged into a single input box, so
muted setup code reads as one block; a statement prints (and gets an output box)
as soon as it is left unterminated by a semicolon or calls `disp`, `printf`, and
the like.

## 4. BISTs (testing suite)

It is **very important** that function files ship a built-in self-test (BIST)
suite that checks for correct output and properly catches error conditions.
BISTs go at the bottom of the file, using `%!` line prefixes.

```
%!test
%! x = [1 2 1 3 2 4 3 2 4 3 2 2];
%! [h, p, stats] = chi2gof (x);
%! assert_equal (h, 0);
%! assert_equal (p, NaN);
%! assert_equal (stats.chi2stat, 0.1205375022748029, 1e-14);

%!error<chi2gof: too few input arguments.> chi2gof ()
%!error chi2gof ([2, 3; 3, 4])
```

A few points on layout and content:

- **Put positive tests first, then a blank line, then the `%!error` tests.**
- The package requires **Octave >= 11**, so use **`assert_equal`** for value
  comparisons — it gives clearer failure messages than the older `assert`.
- In a `%!error<regex>` test, match the **full** error message; this line is
  exempt from the 80-column limit.  Escape regex metacharacters such as
  parentheses (`\(`, `\)`).
- **Only test errors that this function itself emits.**  Do not re-test errors
  raised by core Octave functions — those are covered in core Octave.
- Use `%!shared var1, var2` to declare fixtures shared across following tests.

Append tests as you develop; it saves debugging time and catches marginal errors
that would otherwise come back as bug reports.  It is equally important to add
tests when fixing a bug or extending an existing function.  There is no such
thing as too many tests. :metal:

## 5. Coding style

The package has its own coding style.  It is close in spirit to the conventions
used in GNU Octave itself, but it is **not** identical — there are deliberate
differences (Allman-style block terminators, a specific quoting policy, naming
rules, and so on) that are spelled out below.  Follow these rules for any new
code, and when in doubt copy the layout of an existing neighbouring file.

The overriding goal is consistency: a reader should not be able to tell which
contributor wrote a given file.

### Layout and whitespace

- Keep lines **within 80 columns**.  Wrap long calls and strings with `...`
  continuation, aligning the continuation under the first argument.
- Use `LF` (unix) line endings, **never** `CRLF` (windows).
- **Indent with 2 spaces.  Never use tabs** — anywhere, including compiled
  sources.
- Use `##` for comments.  Do **not** use `%` or `%%` as in MATLAB.  (`%!` is
  reserved for BIST/demo blocks; see §3 and §4.)
- When **calling or defining** a function, put a space after each comma and
  **before the opening parenthesis**.  This is the single most visible
  Octave-vs-MATLAB tell — do not omit it:

```
x = max (sin (y + 3), 2);
function out = foo (a, b)
```

- The exception is **matrix and cell constructors**, where a space before the
  paren would be parsed as a separator and split the element in two:

```
a = [sin(x), cos(x)];
b = {sin(x), cos(x)};
```

- For an **indexing** expression, do **not** put a space after the identifier —
  this distinguishes indexing from a function call.  The space after a comma is
  optional for simple indices but recommended for complex ones:

```
A(:,i,j)
A([1:i-1; i+1:n])
```

### Operators and control flow

- Use `!` for logical NOT (not `~`) and `!=` for not-equal (not `~=`):

```
a != 0;
b(! isnan (a)) = [];
```

- Enclose `if`, `while`, `until`, and `switch` conditions in parentheses:

```
if (isvector (a))
  s = sum (a);
endif
```

  Do **not** parenthesise the iteration counter of a `for` statement:

```
for ii = 1:numel (a)
  b(ii) = sum (a(:,ii));
endfor
```

- Always close a block with its **specific typed keyword** — `endif`, `endfor`,
  `endwhile`, `endswitch`, `endfunction`, `endclassdef`, `endproperties`,
  `endmethods` — never the generic `end`.  (`end` still means the last index
  inside an indexing expression, e.g. `x(end)`.)

### Quoting

- **Single quotes are the default** for every character vector and option
  string: `'Format'`, `'off'`, `anova2 (x, reps, 'off')`,
  `set (gca, 'XLim', xl)`.
- **Double quotes** are used **only** for the string literals passed to the
  message functions `error`, `warning`, `printf`, `fprintf`, and `sprintf` —
  always, whether or not the string contains an escape (`\n`) or a format
  specifier (`%d`).  Keeping these double-quoted lets apostrophes read
  naturally:

```
error ("anova1: X must be numeric.");                      % yes
error ("Prediction must be 'curve' or 'observation'.");    % yes
error ('Prediction must be ''curve'' or ''observation''.'); % no
```

- **Cell arrays of character vectors** always use single quotes:
  `{'red', 'green', 'blue'}`, never `{"red", "green", "blue"}`.
- In a `switch`/`case`, a single-char-vector label uses single quotes —
  `case 'first'`, and `case {'first', 'last'}` for several.

### Naming

| Kind                     | Convention   | Examples                          |
|--------------------------|--------------|-----------------------------------|
| Class / type             | lowercase    | `cvpartition`, `LinearModel`\*    |
| Public property          | PascalCase   | `Formula`, `NumObservations`      |
| Local variable           | camelCase    | `optNames`, `inputFormat`         |
| Loop index               | `ii` (or `i`)| `for ii = 1:n`                    |
| Object handle in methods | `this`       | `function disp (this)`            |
| Private helper (`private/`) | `__name__.m` | `__disp__`, `__paramci__`      |
| Compiled oct-file (`src/`)  | `__name__.cc`| `__editdist__.cc`              |

Double-underscore wrapping (`__foo__`) marks a function as **internal / not
user-facing**.  (\*Established statistics classes such as `LinearModel` keep
their historical MATLAB-compatible casing; new lower-level types follow the
lowercase rule.)

### Error messages

Error messages follow the format **`"<scope>: <message>."`**:

- Prefix with the function name (`"anova1: ..."`) — or, inside a class method,
  `"<class>.<method>: ..."` (e.g. `"LinearModel.predict: ..."`).
- Start lower-case and **end with a period**.
- Refer to arguments by their UPPERCASE documentation name
  (`"... X must be numeric."`).
- Build long messages with `strcat (...)` across `...`-continued lines, aligned
  under the first fragment — **never** with `[...]` bracket concatenation:

```
error (strcat ("anova1: GROUP must be a vector of the same", ...
               " length as X."));
```

- When forwarding a computed message, use `error ("scope: %s", errmsg)`.

### Function and file structure

A function file is laid out in this order:

1. GPL header (§1).
2. Texinfo documentation block (§2).
3. The `function` line, then **input validation first**, under an
   `## Input validation` comment: guard `nargin` / `nargout`, then type- and
   shape-check each argument, erroring early.
4. The body, with `##` full-line comments marking the logical steps.
5. `endfunction`.
6. Any local helper functions.
7. The `%!demo` blocks (§3), then the `%!test` / `%!error` blocks (§4).

```
function p = anova1 (x, group, displayopt)

  ## Input validation
  if (nargin < 1)
    error ("anova1: too few input arguments.");
  endif
  ...

  ## Compute the group means
  ...
endfunction
```

Default argument values may be set in the signature, e.g.
`function [err, days] = hms2days (H, MI, S, MS = 0)`.

### classdef classes

Inside a class file, order the sections as: the `classdef` line and class
`@deftp` block; `properties` blocks (public first, then attributed blocks such
as `properties (SetAccess = private, Hidden)`), each documented property
carrying its own `@deftp` block; `methods (Hidden)` for `disp`/`display` and
internals; one or more `methods (Access = public)` blocks grouped by theme;
`methods (Access = private)`; `endclassdef`; and finally any local helper
functions after `endclassdef`.  See §2 for the `@deftp`-vs-`@deftypefn`
documentation split.

### Compiled sources (`src/`)

Oct-files follow the same spirit with a few C++ specifics:

- GPL header in a `/* ... */` block (same text as §1).
- `#include <octave/oct.h>`; define functions with `DEFUN_DLD`.
- **2-space indentation, no tabs**; opening brace of a free function on its own
  line.
- `//` line comments; built via `src/Makefile`.

### Keeping metadata in sync

When you add a user-facing function or class, update **`INDEX`** — add the new
name under the right category heading (this drives the function index and the
online docs).  Do **not** touch `DESCRIPTION`: its `Version` and `Date` are
bumped only at release time, by the maintainer.

The package also ships a set of `doc-cache` files (one per function directory).
These feed `lookfor` and let `pkg install` skip a slow regeneration step, so
keep them current.  Do this whenever you:

- add a new function or class, or
- change the **top-level help docstring** of an existing function or class.

Use the repo-root helper (build the oct-files first with `make -C src` so
`src/doc-cache` is complete):

```
./regen-doc-cache.sh                 # rebuild every directory (slow, ~minutes)
./regen-doc-cache.sh anova1          # update just this function/class entry (fast)
./regen-doc-cache.sh inst/dist_fun   # rebuild one directory's cache
```

The single-name form is the fast path: it edits only that one entry in the
existing cache, so a docstring tweak costs a second instead of several minutes.

Changes to a **method's or property's** docstring do **not** require a
regeneration — those never enter the cache (and `help` reads them straight from
the source file regardless).

---

**That's about it!**  Keep more or less consistent with the other function files
in the package and you should be fine. :smile:  If anything here is unclear,
open an issue or ask in your pull request — we are happy to help.
