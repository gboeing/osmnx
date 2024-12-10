# OSMnx tests

Read more about the project's standards in the [contributing guidelines](../CONTRIBUTING.md) and ensure that you have installed the necessary [dependencies](../environments/tests/env-ci.yml) for the test suite.

## Code format

Format the code per the project's style by running the pre-commit hooks:

```shell
pre-commit install
pre-commit run -a
```

## Run tests

Run the test suite locally by running (from the repository root):

```shell
bash ./tests/lint_test.sh
```

## Continuous integration

Pull requests trigger continuous integration tests via GitHub Actions (see [workflow](../.github/workflows/ci.yml)), including the following steps:

- build the docs
- check code formatting
- lint the code and docstrings
- type check the code
- run tests and coverage

## Releases

To publish a new version, update `CHANGELOG.md` and edit the version number in `osmnx/_version.py`. If necessary, update the dates in `LICENSE.txt` and `docs/source/conf.py` and the dependency versions in `pyproject.toml`. Then tag the repository with its new semantic version, like `v1.2.3`.

Pushing the tags will trigger Github Actions to publish the distribution to PyPI (see [workflow](../.github/workflows/build-publish-pypi.yml)) and publish a new image to Docker Hub (see [workflow](../.github/workflows/build-publish-docker.yml)). The `regro-cf-autotick-bot` will open a pull request to update the conda-forge [feedstock](https://github.com/conda-forge/osmnx-feedstock): merge that PR to publish the distribution on conda-forge. Finally, update the [Examples Gallery](https://github.com/gboeing/osmnx-examples) to use the new version.
