# Contributing guidelines

Thanks for using OSMnx and for considering contributing to it by opening an issue or pull request. Every piece of software is a work in progress. This project is the result of many hours of work contributed freely its contributors and the many people that build the projects on which it depends. Thank you for contributing!

## If you have a "how-to" or usage question

Please ask your question on [StackOverflow](https://stackoverflow.com/search?q=osmnx), as we reserve the issue tracker for bug reports and new feature development. Any such questions asked in the issue tracker will be automatically closed.

## If you have an installation problem

Ensure you have followed the installation instructions in the [documentation](https://osmnx.readthedocs.io/). If you installed OSMnx via conda-forge, please open an issue at its [feedstock](https://github.com/conda-forge/osmnx-feedstock/issues).

## If you have a feature proposal

The OSMnx project follows three principles when considering new functionality: 1) it is useful for a broad set of users, 2) it generalizes well, and 3) it is not trivially easy for users to implement themselves.

- Post your proposal on the [issue tracker](https://github.com/gboeing/osmnx/issues), and _provide all the information requested in the template_, so we can review it together (some proposals may not be a good fit for the project).
- Fork the repo, make your change, update the [changelog](./CHANGELOG.md), run the [tests](./tests), and submit a pull request.
- Adhere to the project's code and docstring standards by running its [pre-commit](.pre-commit-config.yaml) hooks.
- Respond to code review.

## If you found a bug

- Read the error message, then review the [documentation](https://osmnx.readthedocs.io/) and OSMnx [Examples Gallery](https://github.com/gboeing/osmnx-examples), which cover key concepts, installation, and package usage.
- Search through the open and closed [issues](https://github.com/gboeing/osmnx/issues) to see if the problem has already been reported.
- If the problem is with a dependency of OSMnx, open an issue in that dependency's repo.
- If the problem is with OSMnx itself and you can fix it simply, please open a pull request.
- If the problem persists, please open an issue in the [issue tracker](https://github.com/gboeing/osmnx/issues), and _provide all the information requested in the template_, including a minimal standalone example so others can independently and completely reproduce the bug.

## Creating a minimal standalone reproducible example

We need a minimal standalone example code snippet to be able to reproduce and in turn troubleshoot your problem (also provide the resulting complete error traceback if your code generates an error). This code snippet must be:

*Minimal*: the absolute fewest lines of code necessary to reproduce your problem without any extraneous code or data unrelated to your specific problem. This usually requires reworking your code and data down to just a few lines necessary to generate the error.

*Standalone*: all imports, data, and variables must be completely defined within the snippet itself so others can independently run it from top to bottom by copying/pasting it into a Python interpreter. Do not link to or load external files and do not provide screenshots of code or error messages: provide all code, data, and tracebacks inline as text.

If you're unsure how to create a good reproducible example, read
        [this guide](https://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports).
