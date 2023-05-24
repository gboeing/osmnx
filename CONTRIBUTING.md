# OSMnx contributing guidelines

Thanks for using OSMnx and for considering contributing to it by opening an issue or pull request. If you open an issue, fill out every section of the new issue template. Your issue will be automatically closed if you don't provide all of the information requested in the template, because we need this information to address your issue. Please note that we reserve the issue tracker for reporting bugs and proposing new features, not for asking usage questions.

#### If you have a "how-to" or usage question rather than a bug report or new feature proposal:

  - please ask your question on [StackOverflow](https://stackoverflow.com/search?q=osmnx)

#### If you're having an installation problem:

  - make sure you've followed the installation instructions in the [documentation](https://osmnx.readthedocs.io/)
  - if you installed OSMnx via conda-forge, please open an issue at its [feedstock](https://github.com/conda-forge/osmnx-feedstock/issues)

#### If you've found a bug:

  - read the error message and [documentation](https://osmnx.readthedocs.io/)
  - search through the [open issues](https://github.com/gboeing/osmnx/issues?q=is%3Aopen+is%3Aissue) and [closed issues](https://github.com/gboeing/osmnx/issues?q=is%3Aissue+is%3Aclosed) to see if it has already been reported
  - if the problem is with a dependency of this project, open an issue in the dependency's repo
  - if the problem is with OSMnx and you can fix it simply, please open a pull request
  - if the problem persists, please open an issue in the [issue tracker](https://github.com/gboeing/osmnx/issues), and *provide all the information requested in the new issue template*, including a minimal working example so others can independently and completely reproduce the problem

#### If you have a feature proposal or want to contribute:

  - post your proposal on the [issue tracker](https://github.com/gboeing/osmnx/issues), and *provide all the information requested in the new issue template*, so we can review it together (some proposals may not be a good fit for the project)
  - fork the repo, make your change, [test it](./tests), and submit a PR
  - respond to code review
  - adhere to the following project standards
    - `black` code style with max line length of 100
    - `isort` sorted imports
    - `numpy` style docstrings

This project requires minimum Python and NumPy versions in accordance with [NEP 29](https://numpy.org/neps/nep-0029-deprecation_policy.html).

Every piece of software is a work in progress. This project is the result of many hours of work contributed freely by myself and the many people that build the projects it depends on. Thank you for contributing!
