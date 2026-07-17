# OSMnx tests

Read more about the project's standards in the [contributing guidelines](../CONTRIBUTING.md) and ensure that you have installed the necessary dev [dependencies](../pyproject.toml) for the test suite.

## Code format

Format the code per the project's style by running the pre-commit hooks:

```shell
pre-commit install
pre-commit run -a
```

## Run tests

Run the test suite locally by running (from the repository root):

```shell
uv sync --all-extras --all-groups
bash ./tests/lint_test.sh
```

The first test run populates `tests/.cache/` with one worker then marks it ready. Later runs use four workers and the cache. Delete `tests/.cache/` to force the next run to rebuild the cache. Set `USE_PERSISTENT_CACHE=false` to bypass the persistent cache and test live APIs with one worker.

## Continuous integration

Pull requests trigger continuous integration tests via GitHub Actions (see [workflow](../.github/workflows/ci.yml)), including the following steps:

- build the docs and the package
- check code formatting
- lint the code and docstrings
- type check the code
- run tests and coverage

## Releases

To publish a new version, update `CHANGELOG.md` and edit the version number in `pyproject.toml`. If necessary, update the dates in `LICENSE.txt` and `docs/source/conf.py` and the dependency versions in `pyproject.toml`. Then tag the repository with its new semantic version, like `git tag -a "v1.2.3" -m "v1.2.3"`.

Pushing the tags will trigger Github Actions to publish the distribution to PyPI (see [workflow](../.github/workflows/build-publish-pypi.yml)) and publish a new image to Docker Hub (see [workflow](../.github/workflows/build-publish-docker.yml)). The `regro-cf-autotick-bot` will open a pull request to update the conda-forge [feedstock](https://github.com/conda-forge/osmnx-feedstock): merge that PR to publish the distribution on conda-forge. Finally, update the [Examples Gallery](https://github.com/gboeing/osmnx-examples) to use the new version.
