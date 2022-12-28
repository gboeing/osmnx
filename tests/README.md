# OSMnx tests

First, ensure that you have installed the necessary [dependencies](environment-dev.yml). Then use the repo's [pre-commit hooks](../.pre-commit-config.yaml) and/or the scripts in this folder to:

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

## Release a new version

Update `CHANGELOG.md` and edit the version number in `setup.py`, `osmnx/_version.py`, and `docs/source/conf.py`. If needed, update `LICENSE.txt` dates and `requirements.txt` dependency versions. Then change directories to the repository's root and run:

```
bash ./tests/packaging.sh
```

This will tag the repo with the new version number, upload the PyPI distribution, and update the conda-forge feedstock. Then, open a pull request at the [feedstock](https://github.com/conda-forge/osmnx-feedstock/pulls) to release on conda-forge. Finally, when the new version is available on conda-forge, update the [Docker image](../environments/docker) and the [examples repo](https://github.com/gboeing/osmnx-examples) to use the new version.
