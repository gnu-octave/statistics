# statistics

This is the official repository for the Statistics package for GNU Octave.

**Content:**

1. [About](#1-about)
2. [Planning](#2-planning)
3. [Documentation](#3-documentation)
4. [Install statistics](#4-install-statistics)
5. [Provide feedback](#5-provide-feedback)
6. [Contribute](#6-contribute)
7. [Google Summer of Code](#7-google-summer-of-code)

## 1. About

The **statistics** package is a collection of functions for statistical analysis. As with GNU Octave, the **statistics** package aims to be mostly compatible with MATLAB's equivalent Statistics and Machine Learning Toolbox. However, this is not always applicable, desirable, or even possible. Hence, identical (in name) functions do not necessarily share the same functionality or behavior. Nevertheless, they produce consistent and correct results, unless there is a bug: see [Murphy's Law](https://en.wikipedia.org/wiki/Murphy's_law) :smile:.

## 2. Planning

As of 10.6.2022, the development of the **statistics** package was moved from [SourceForge](https://octave.sourceforge.io/statistics/) and [Mercurial](https://en.wikipedia.org/wiki/Mercurial) to [GitHub](https://github.com/gnu-octave/statistics) and [Git](https://en.wikipedia.org/wiki/Git). Given the opportunity of this transition, the package has been redesigned, as compared to the its previous point [release 1.4.3](https://octave.sourceforge.io/download.php?package=statistics-1.4.3.tar.gz) at SourceForge, with the aim to keep its structure simplified and easier to maintain. To this end, two major decisions were made:
- Keep a single dependency to the last two major point releases of GNU Octave.
- Deprecate old functions once their fully Matlab compatible equivalents are implemented.

As of 19.10.2025, landmarked by point [release 1.7.6](https://github.com/gnu-octave/statistics/releases/tag/release-1.7.6), a new road map was laid out for future development. During this 3.5 years period, after 2356 individual commits, 197 merged PRs, and 63 closed issues, 19 package releases were made, introducing considerable improvements to the **statistics** package. More importantly, there were no longer remaining old functions to deprecate and the rule of a single dependency to core Octave was prohibiting further improvements with respect to MATLAB compatibility. Hence, it was decided to introduce a package dependency to the [`datatypes`](https://github.com/pr0m1th3as/datatypes) package in order to keep advancing the functionality of existing functions (and enhancing their MATLAB compatibility) as well as introduce new functionality, which requires new datatypes essential to statistical analysis, such as categorical values.

As of 30.11.2025, landmarked by [release 1.8.0](https://github.com/gnu-octave/statistics/releases/tag/release-1.8.0), the `datatypes` package is a required dependency.

As of 18.05.2026, landmarked by [release 1.8.3](https://github.com/gnu-octave/statistics/releases/tag/release-1.8.3), **Octave 11.1** is required and it has been decided that the **statistics** package will update its Octave dependency to the latest major release at least 3 months after its official release to allow sufficient time for users to upgrade their core Octave installation.

As of 26.08.2026, landmarked by [release 1.9.0](https://github.com/gnu-octave/statistics/releases/tag/release-1.9.0), extensive MATLAB compatibility work added 113 new functions and classes across supervised learning, regression, clustering and anomaly detection, dimensionality reduction, distributions, hypothesis testing and plotting. The probability distribution classes moved into the `prob` namespace and several defaults changed to match MATLAB.

As of 30.09.2026, landmarked by [release 2.0.0](https://github.com/gnu-octave/statistics/releases/tag/release-2.0.0) and exactly four years after [release 1.5.0](https://github.com/gnu-octave/statistics/releases/tag/release-1.5.0), the first release made from GitHub, the long-term objective of MATLAB compatibility set out with the move to GitHub is concluded. During the year since release 1.7.6, after 1106 individual commits, 154 merged PRs, and 67 closed issues, 12 package releases were made, bringing the total since the move to GitHub to 31 releases. The **statistics** package is now in maintenance mode: bugs and incompatibilities with MATLAB are fixed as they are reported, and contributions are reviewed as they come.

Some parts of MATLAB's Statistics and Machine Learning Toolbox are left out by decision, as out of scope or too narrow for this package: interactive apps and other mouse-driven tools, automated model search and hyperparameter optimization (`bayesopt`, `fitcauto`, `fitrauto`), automated feature generation (`gencfeatures`, `genrfeatures`), incremental learning, Hamiltonian Monte Carlo sampling, and the `gather` and `update` methods for tall arrays and code generation. Bootstrap sampling and confidence intervals (`bootstrp`, `bootci`) are provided by the [`statistics-resampling`](https://github.com/gnu-octave/statistics-resampling) package. Other functions still missing, such as those for design of experiments and statistical process control, are open to [contributions](#6-contribute).

## 3. Documentation
All functions, class definitions, and their respective methods are documented with [texinfo](https://www.gnu.org/software/texinfo/) format, which can be accessed from the Octave command with the `help` function.  Use dot notation to access the help of methods of classdefs objects. For example:
```
help ClassificationKNN.predict
help prob.BetaDistribution.cdf
```
You can also find the entire documentation of the **statistics** package along with its function index at [https://gnu-octave.github.io/statistics/](https://gnu-octave.github.io/statistics/). Alternatively, you can build the online documentation locally using the [`pkg-octave-doc`](https://github.com/gnu-octave/pkg-octave-doc) package. Assuming both packages are installed and loaded, browse to any directory of your choice with *write* permission and run:
```
package_texi2html ("statistics")
```

Since release 1.9.0, the **statistics** package also ships its documentation as a Qt help file, so it can be browsed and searched in the Documentation tab of Octave's GUI once the package is loaded with `pkg load statistics`.

## 4. Install statistics

To install the latest release, you need Octave (>=11.1.0) installed on your system. Install it by typing:

  `pkg install statistics`

You can automatically download and install the latest development version of the **statistics** package found [here](https://github.com/gnu-octave/statistics/archive/refs/heads/main.zip) by typing:

  `pkg install "https://github.com/gnu-octave/statistics/archive/refs/heads/main.zip"`

If you need to install a specific older release, for example `1.4.2`, type:

  `pkg install "https://github.com/gnu-octave/statistics/archive/refs/tags/release-1.4.2.tar.gz"`

MacOS users, who wish to enable parallel processing for the [`fitcnet`](https://gnu-octave.github.io/statistics/fitcnet.html) function and the [`ClassificationNeuralNetwork`](https://gnu-octave.github.io/statistics/ClassificationNeuralNetwork.html) class, should first link to the OpenMP library before installing the **statistics** package. This will enable the compiler to use multithreading in the [`fcnntrain`](https://gnu-octave.github.io/statistics/fcnntrain.html) and [`fcnnpredict`](https://gnu-octave.github.io/statistics/fcnnpredict.html) dynamically linked functions. To do this, type:

  `setenv ("CPPFLAGS", "-I/opt/homebrew/opt/libomp/include -Xclang -fopenmp")`

After installation, type:
- `pkg load statistics` to load the **statistics** package.
- `news statistics` to review all the user visible changes since last version.
- `pkg test statistics` to run a test suite for all 642 functions and class definitions currently available and ensure that they work properly on your system.

## 5. Provide feedback

The **statistics** package is in maintenance mode, so bug reports are the most valuable feedback. Please, open an [issue](https://github.com/gnu-octave/statistics/issues) if a function returns a wrong result, fails on valid input, or behaves differently from its MATLAB equivalent. Include the code that reproduces the problem, the versions of Octave and of the **statistics** package you use (`ver statistics`), and, for an incompatibility, the output MATLAB gives along with its release. Use [markdown tips](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) to keep your report easy to read.

Requests for missing functions are no longer taken up by the maintainer, and the parts left out by decision (see [Planning](#2-planning)) will not be added. If a missing function matters to you, consider [contributing](#6-contribute) it.

**Please, do NOT use the issue tracker for requesting help.** Use the [discourse group](https://octave.discourse.group/c/help/6) for help with using functions and programming in Octave.

## 6. Contribute

The **statistics** package is **open source** and everyone is welcome to contribute. Pull requests that fix a bug, correct an incompatibility with MATLAB, or add a missing function within the scope of the package are reviewed as they come.

To contribute, [fork](https://github.com/gnu-octave/statistics/fork) this repository, make your changes on a branch, and open a [pull request](https://github.com/gnu-octave/statistics/pulls). Alternatively, you may open an issue and attach a git patch file, which the maintainer will apply.

Please, follow the coding style of the **statistics** package, summarized in [CONTRIBUTING.md](https://github.com/gnu-octave/statistics/blob/main/CONTRIBUTING.md), and include built-in self tests (BISTs) for every change: normal use, edge cases and every error your code raises. Where a function has a MATLAB equivalent, the expected values in the tests should come from MATLAB.

Contributing is not only about fixing bugs or implementing new functions. Improving the help text of a function, or adding tests and demos, is just as welcome. Fixing a typo in the help file is still of value. So don't hesitate to contribute! :+1:

## 7. Google Summer of Code

The **statistics** package took part in [Google Summer of Code](https://summerofcode.withgoogle.com/) under the GNU Octave organization from 2023 to 2026. The following students worked on it:

| Year | Student | Project |
|------|---------|---------|
| 2023 | Azmat Khan | [Addition of regression GAM and kNN classification functionality in statistics package](https://azm999.hashnode.dev/) |
| 2024 | Pallav Purbia | [Adding `ClassificationSVM`, `ClassificationNeuralNetwork` to the Statistics Package of Octave](https://pallavpurbia.blogspot.com/) |
| 2024 | Ruchika Sonagote | [Adding GAM and Discriminant Classification Classes and Implementing Missing Methods](https://ruchikasonagote.blogspot.com) |
| 2025 | Swayam Shah | [Adding clustering \*Searcher classes in statistics package](https://swayams.medium.com/list/google-summer-of-code-25-gnu-octave-100ac1340426) |
| 2026 | Avanish Salunke | [Add `LinearModel` and `CompactLinearModel` objects to the statistics package](https://avanishsalunke.github.io/AvanishSalunke-gsoc-2026/) |
| 2026 | Aman Behera | [Add `anova` Classdef Object to the GNU Octave Statistics Package](https://beingamanforever.hashnode.dev/) |

Thank you all for your work. With release 2.0.0 there is little left to build for the **statistics** package, so it will not take part in Google Summer of Code for the foreseeable future.
