# OSMnx tests

First, ensure that you have installed the necessary [dependencies](environment-dev.yml). Then use the scripts in this folder to:

  - format the code according to the project's style
  - lint the code
  - lint the docstrings
  - run unit tests and coverage

You can read more about the project's standards and code/docstring style in the [contributing guidelines](../CONTRIBUTING.md).

## Code format

Format the code and sort imports according to the project's style by changing directories to the repository's root and running:

```
bash ./tests/black.sh
```

## Lint and test

Lint and test the code and docstrings by changing directories to the repository's root and running:

```
bash ./tests/lint_test.sh
```

## Continuous integration

All PRs trigger continuous integration tests via GitHub Actions. See the [configuration](../.github/workflows/tests.yml). The following tests are automatically run:

  - build the docs
  - check code formatting
  - docstrings linter
  - code linter
  - unit tests and coverage
