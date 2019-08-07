afmize
====
[![Build Status](https://www.travis-ci.com/ToruNiina/afmize.svg?branch=master)](https://www.travis-ci.com/ToruNiina/afmize)
[![version](https://img.shields.io/github/release/ToruNiina/afmize.svg?style=flat)](https://github.com/ToruNiina/afmize/releases)
[![license](https://img.shields.io/github/license/ToruNiina/afmize.svg?style=flat)](https://github.com/ToruNiina/afmize/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2556445.svg)](https://doi.org/10.5281/zenodo.2556445)

making pseudo AFM images from structure file.

## Usage

```sh
$ afmize example.toml
```

### Configurations

It uses [toml format](https://github.com/toml-lang/toml) as a config file format.
For more information about the format, see [the spec of toml-v0.5.0](https://github.com/toml-lang/toml/blob/master/versions/en/toml-v0.5.0.md).

You can change...
- probe size
- radii of atoms
- x, y, and z resolutions

The content would be like this.

```toml
# toml format v0.5.0
# the default unit is: degree for angle, angstrom for length
# you can also specify the unit explicitly
# pm, angstrom, nm, um, mm are allowed
file.input   = "example.pdb"
file.output  = "output"
probe.size   = {radius = "1.0nm", angle = 10.0}
resolution.x = "1.0nm"
resolution.y = "1.0nm"
resolution.z = "0.64Å"
range.x      = ["0.0nm", "100.0nm"]
range.y      = ["0.0nm", "100.0nm"]

scale_bar.length = "5.0nm"

colormap.min =  "0.0nm"
colormap.max = "10.0nm"

stage.position = "1.0Å" # z coordinate of stage position
stage.align    = true   # if true, afmize moves the structure to make the bottom of the bounding box is equal to stage.position.

radii.atom.H = "0.1angstrom"  # radius of all the elements
radii.residue.ARG.CA = "3.0Å" # radius for a specific pair of residue and atom.
# radii of the other atoms in the residue ARG will be the default values.
```

### Supported file formats

- xyz
- pdb
  - note: it may fail reading pdb files that does not conform [wwPDB 3.3](http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM).

## Installation

### Prerequisites

Use your favorite package managers (e.g. `apt`) to install them.

- CMake
- git
  - to download submodules
- C++14 compliant compiler. tested with ...
  - g++-7 or later
  - clang++-6 or later
  - or others that fully support C++14 ISO standard.

All other dependencies are managed by CMake script and git submodule.

### Building

```sh
$ git submodule update --init
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make test # optional
```

After this, you will find the executable at `bin/` directory.

## Licensing Terms

This product is licensed under the terms of the [MIT License](LICENSE).

- Copyright (c) 2018-2019 Toru Niina

All rights reserved.
