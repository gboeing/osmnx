# OSMnx tests

First, ensure that you have installed the dependencies in `requirements-dev.txt`. Then use the scripts in this folder to:

  - format the code according to the project's style
  - lint the code
  - lint the docstrings
  - run unit tests and coverage

You can read more about the project's standards and code/docstring style in the [contributing guidelines](../CONTRIBUTING.md).

## Code format

Format the code according to the project's style by changing directories to the repository's root and running: 

```
bash ./tests/run_black.sh
```

## Run tests

Run the tests and linters by changing directories to the repository's root and running:

```
bash ./tests/run_tests.sh
```
