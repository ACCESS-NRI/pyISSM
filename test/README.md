# pyISSM testing


This directory is home to a suite of pyISSM tests. The current test suite consists of the following key directories:

- `/assets`: Common assets used across all tests. This includes data files, *.exp files and parameter files.
- `/ci_tests`: Tests automatically executed for GitHub CI purposes, testing core functionality of pyISSM and ISSM.
- `/general`: General tests that demonstrate pyISSM / ISSM functionality.
- `/pytests`: Pytest-based unit tests for pyISSM classes and tools.

Both `ci_tests` and `general` directories have a `dev` subdirectory that provides a "staging" area for in-progress tests that may not be fully functional.

The general directory structure is below:
```
.
├── assets
│   ├── Data
│   │   └── ...
│   ├── Exp
│   │   └── ...
│   └── Par
│       └── ...
├── ci_tests
│   ├── dev
│   │   └── ...
│   └── ...
├── general
│   ├── dev
│   │   └── ...
│   └── ...
├── pytests
│   └── unit
│       ├── pytest.ini
│       └── tests
│           └── ...
└── ...
```
