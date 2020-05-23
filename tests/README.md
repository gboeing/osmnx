# OSMnx tests

Use the scripts in this folder to:

  - format the code according to the project's style
  - lint the code
  - lint the docstrings (numpy style)
  - run unit tests and coverage

First, ensure that you have installed the dependencies in `requirements-dev.txt`.

## Code format

Run black to format the code according to the project's style by changing directories to the repository's root and running: 

```
bash ./tests/run_black.sh
```

## Run tests

Run the tests and linters by changing directories to the repository's root and running:

```
bash ./tests/run_tests.sh
```
