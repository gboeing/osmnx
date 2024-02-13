# Contributing guidelines

Thanks for using OSMnx and for considering contributing to it by opening an issue or pull request. Every piece of software is a work in progress. This project is the result of many hours of work contributed freely its contributors and the many people that build the projects on which it depends. Thank you for contributing!

## If you have a "how-to" or usage question

Please ask your question on [StackOverflow](https://stackoverflow.com/search?q=osmnx), as we reserve the issue tracker for bug reports and new feature development. Any such questions asked in the issue tracker will be automatically closed.

## If you have an installation problem

Ensure you have followed the installation instructions in the [documentation](https://osmnx.readthedocs.io/). If you installed OSMnx via conda-forge, please open an issue at its [feedstock](https://github.com/conda-forge/osmnx-feedstock/issues).

## If you've found a bug

- Read the error message, then review the [documentation](https://osmnx.readthedocs.io/) and OSMnx [Examples Gallery](https://github.com/gboeing/osmnx-examples), which cover key concepts, installation, and package usage.
- Search through the open and closed [issues](https://github.com/gboeing/osmnx/issues) to see if the problem has already been reported.
- If the problem is with a dependency of OSMnx, open an issue in that dependency's repo.
- If the problem is with OSMnx itself and you can fix it simply, please open a pull request.
- If the problem persists, please open an issue in the [issue tracker](https://github.com/gboeing/osmnx/issues), and _provide all the information requested in the template_, including a minimal working example so others can independently and completely reproduce the bug.

## If you have a feature proposal

- Post your proposal on the [issue tracker](https://github.com/gboeing/osmnx/issues), and _provide all the information requested in the template_, so we can review it together (some proposals may not be a good fit for the project).
- Fork the repo, make your change, update the [changelog](./CHANGELOG.md), run the [tests](./tests), and submit a pull request.
- Adhere to the project's code and docstring standards by running its [pre-commit](.pre-commit-config.yaml) hooks.
- Respond to code review.
