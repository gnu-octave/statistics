# statistics

This is the official repository for the Statistics package for GNU Octave.

**Content:**

1. About
2. Documentation
3. Install statistics
4. Provide feedback
5. Contribute

## 1. About

The **statistics** package is a collection of functions for statistical analysis. As with GNU Octave, the **statistics** package aims to be mostly compatible with MATLAB's equivalent Statistics and Machine Learning Toolbox. However, this is not always applicable of even possible. Hence, identical (in name) functions do not necessarily share the same functionality or behavior. Nevertheless, they produce consistent and correct results, unless there is a bug: see [Murphy's Law](https://en.wikipedia.org/wiki/Murphy's_law) :smile:.

As of 10.6.2022, the developemnt of the **statistics** package was moved from [SourceForge](https://octave.sourceforge.io/statistics/) and [Mercurial](https://en.wikipedia.org/wiki/Mercurial) to [GitHub](https://github.com/gnu-octave/statistics) and [Git](https://en.wikipedia.org/wiki/Git). Given the opportunity of this transition, the package has been redesigned, as compared to the its previous point [release 1.4.3](https://octave.sourceforge.io/download.php?package=statistics-1.4.3.tar.gz) at SourceForge, with the aim to keep its structure simplified and easier to maintain. To this end, two major decisions have been made:
- Keep a single dependency to the last two major point releases of GNU Octave.
- Deprecate old functions once their fully Matlab compatible equivalents are implemented.

## 2. Documentation
All functions, class definitions, and their respective methods are documented with [texinfo](https://www.gnu.org/software/texinfo/) format, which can be accessed from the Octave command with the `help` function.  Use dot notation to access the help of a particular method and new classdefs or the relative path for old style classes. For example:
```
help ClassificationKNN.predict
help @cvpartition/test
```
You can also find the entire documentation of the **statistics** package along with its function index at [https://gnu-octave.github.io/statistics/](https://gnu-octave.github.io/statistics/). Alternatively, you can build the online documentation locally using the [`pkg-octave-doc`](https://github.com/gnu-octave/pkg-octave-doc) package. Assuming both packages are installed and loaded, browse to any directory of your choice with *write* permission and run:
```
package_texi2html ("statistics")
```

## 3. Install statistics

To install the latest release (1.7.0) you need Octave (>=7.2.0) installed on your system. Install it by typing:

  `pkg install -forge statistics`

You can automatically download and install the latest development version of the **statistics** package found [here](https://github.com/gnu-octave/statistics/archive/refs/heads/main.zip) by typing:

  `pkg install "https://github.com/gnu-octave/statistics/archive/refs/heads/main.zip"`

If you need to install a specific older release, for example `1.4.2`, type:

  `pkg install "https://github.com/gnu-octave/statistics/archive/refs/tags/release-1.4.2.tar.gz"`

After installation, type:
- `pkg load statistics` to load the **statistics** package.
- `news statistics` to review all the user visible changes since last version.
- `pkg test statistics` to run a test suite for all 452 [^1] functions and class definitions currently available and ensure that they work properly on your system.

[^1]: Several functions are still missing from the statistics package, but you are welcome to [contribute](https://github.com/gnu-octave/statistics/blob/main/CONTRIBUTING.md)!

## 4. Provide feedback

You are encouraged to provide feedback regarding possible bugs, missing features[^2], discrepancies or incompatibilities with Matlab functions. You may open an [issue](https://github.com/gnu-octave/statistics/issues) to open a discussion to your particular case. **Please, do NOT use the issue tracker for requesting help.** Use the [discourse group](https://octave.discourse.group/c/help/6) for requesting help with using functions and programming in Octave.

Please, make sure that when reporting a bug you provide as much information as possible for other users to be able to replicate it. Use [markdown tips](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) to make your post clear and easy to read and understand your issue.

[^2]: Don't open an issue just for requesting a missing function! Implement it yourself and make an invaluable contribution :innocent:

## 5. Contribute

The **statistics** package is **open source**! Everyone is welcome to contribute.

There are currently a few open [issues](https://github.com/gnu-octave/statistics/issues) that you can help fixing.

If you find a bug and fix it, just [clone](https://github.com/gnu-octave/statistics.git) this repo with `git clone https://github.com/gnu-octave/statistics.git`, make your changes and add a [pull](https://github.com/gnu-octave/statistics/pulls) request. Alternatively, you may open an issue and add a git-patch file, which will be patched by the maintainer.

Make sure you follow the coding style already used in the **statistics** package (similar to GNU Octave). For a summary of the coding style rules used in the package see [Contribute](https://github.com/gnu-octave/statistics/blob/main/CONTRIBUTING.md).

Contributing is not only about fixing bugs or implementing new functions. Improving the texinfo of the functions help files or adding BISTs and demos at the end of the function files is also important. Fixing a typo in the help file is still of value though. So don't hesitate to contribute! :+1:

