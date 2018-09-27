afmize
====

making pseudo AFM images from structure file.

## Usage

```sh
$ afmize example.toml # generates image file
```

### Configurations

It uses [toml format](https://github.com/toml-lang/toml) as a config file format.
For more information about the format, see [the spec of toml-v0.5.0](https://github.com/toml-lang/toml/blob/master/versions/en/toml-v0.5.0.md).

You can change...
- probe size
- radii of atoms
- x, y, and z resolutions

The content would be like this.

```
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

stage.align    = true
stage.position = "1.0Å" # z coordinate of stage

radii.H      = 0.1 # you can set radii of all the elements found in pdb file.
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
- C++14 compliant compiler
  - gcc 5.1 or later
  - clang 3.4 or later
  - or others that fully support C++14 ISO standard.

All other dependencies are managed by CMake script and git submodule.

### Building

```sh
$ git submodule update
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make test # optional
```

After this, you will find the executable at `bin/` directory.

## Licensing Terms

This product is licensed under the terms of the [MIT License](LICENSE).

- Copyright (c) 2018 Toru Niina

All rights reserved.
