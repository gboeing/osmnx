# Contributing guidelines

Thanks for using OSMnx and for considering contributing to it by opening an issue or pull request. Every piece of software is a work in progress. This project is the result of many hours of work contributed freely its contributors and the many people that build the projects it depends on. Thank you for contributing!

#### If you have a "how-to" or usage question:

- please ask your question on [StackOverflow](https://stackoverflow.com/search?q=osmnx), as we reserve the issue tracker for bug reports and new feature development. Any such questions asked in the issue tracker will be automatically closed.

#### If you're having an installation problem:

- make sure you've followed the installation instructions in the [documentation](https://osmnx.readthedocs.io/)
- if you installed OSMnx via conda-forge, please open an issue at its [feedstock](https://github.com/conda-forge/osmnx-feedstock/issues)

#### If you've found a bug:

- read the error message, then review the [documentation](https://osmnx.readthedocs.io/) and [OSMnx examples](https://github.com/gboeing/osmnx-examples) gallery, which cover key concepts, installation, and package usage
- search through the [open issues](https://github.com/gboeing/osmnx/issues?q=is%3Aopen+is%3Aissue) and [closed issues](https://github.com/gboeing/osmnx/issues?q=is%3Aissue+is%3Aclosed) to see if the problem has already been reported
- if the problem is with a dependency of OSMnx, open an issue in the dependency's repo
- if the problem is with OSMnx itself and you can fix it simply, please open a pull request
- if the problem persists, please open an issue in the [issue tracker](https://github.com/gboeing/osmnx/issues), and _provide all the information requested in the template_, including a minimal working example so others can independently and completely reproduce the problem

#### If you have a feature proposal or want to contribute:

- post your proposal on the [issue tracker](https://github.com/gboeing/osmnx/issues), and _provide all the information requested in the template_, so we can review it together (some proposals may not be a good fit for the project)
- fork the repo, make your change, update the [changelog](./CHANGELOG.md), run the [tests](./tests), and submit a PR
- adhere to the project's code and docstring standards by running its [pre-commit](.pre-commit-config.yaml) hooks
- respond to code review

This project requires minimum Python and NumPy versions in accordance with [NEP 29](https://numpy.org/neps/nep-0029-deprecation_policy.html).
