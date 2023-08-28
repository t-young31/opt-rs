[![Test](https://github.com/t-young31/opt-rs/actions/workflows/test.yml/badge.svg)](https://github.com/t-young31/opt-rs/actions/workflows/test.yml) [![codecov](https://codecov.io/gh/t-young31/mors/branch/main/graph/badge.svg?token=5KTYG2WJ9L)](https://codecov.io/gh/t-young31/mors)

![alt text](src/common/logo.png)

**optrs** is a lightweight molecular mechanics optimisation code written in rust.

***
### Installation
Install the binary with

```bash
pip install optrs
```

or the [Python API](https://github.com/t-young31/opt-rs/tree/main/api).

***
### Usage
Optimise a molecule provided as a [.xyz](https://en.wikipedia.org/wiki/XYZ_file_format) file

```bash
optrs molecule.xyz
```

which will generate _opt.xyz_ in working directory.
[API examples](https://github.com/t-young31/opt-rs/tree/main/api/examples).

***
### Features
- [UFF](https://doi.org/10.1021/ja00051a040) forcefield
- [RB](https://doi.org/10.1002/anie.202011941) forcefield
- 3D structure generation
- Steepest decent optimisation
