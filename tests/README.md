# OSMnx tests

First, ensure that you have installed the necessary [dependencies](../tests/environments/env-ci.yml) for the test suite. Then use the repository's [pre-commit hooks](../.pre-commit-config.yaml) and the scripts in this folder to:

- format the code and docstrings per the project's style
- lint the code and docstrings
- type check the code
- run tests and coverage

Read more about the project's standards in the [contributing guidelines](../CONTRIBUTING.md).

## Code format

Format the code and sort imports per the project's style by running (from the repository root):

```shell
bash ./tests/format.sh
```

## Run tests

Lint, type check, and test the code/docstrings by running (from the repository root):

```shell
pre-commit install
bash ./tests/lint_test.sh
```

## Continuous integration

Pull requests trigger continuous integration tests via GitHub Actions. See the [configuration](../.github/workflows/ci.yml). This includes the following steps:

- build the docs
- check code formatting
- lint the code and docstrings
- type check the code
- run tests and coverage

## Releases

To package and release a new version, update `CHANGELOG.md` and edit the version number in `osmnx/_version.py`. If necessary, update the dates in `LICENSE.txt` and `docs/source/conf.py` and the dependency versions in `pyproject.toml`. Then change directories to the repository's root and run:

```shell
bash -i ./tests/packaging.sh
```

This will tag the repository with the new version number, upload the PyPI distribution, and update the conda-forge feedstock. Then, open a pull request at the [feedstock](https://github.com/conda-forge/osmnx-feedstock) to release on conda-forge. Finally, when the new version is available on conda-forge, update the [Docker image](../environments/docker) and the OSMnx [Examples Gallery](https://github.com/gboeing/osmnx-examples) to use the new version.
