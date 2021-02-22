afmize
====
[![Build Status](https://www.travis-ci.com/ToruNiina/afmize.svg?branch=master)](https://www.travis-ci.com/ToruNiina/afmize)
[![version](https://img.shields.io/github/release/ToruNiina/afmize.svg?style=flat)](https://github.com/ToruNiina/afmize/releases)
[![license](https://img.shields.io/github/license/ToruNiina/afmize.svg?style=flat)](https://github.com/ToruNiina/afmize/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3779937.svg)](https://doi.org/10.5281/zenodo.3779937)

Make pseudo AFM images from a structure file.

## Usage

All the configurations are written in `.toml` file.
No command-line option is available.

```sh
$ afmize example.toml
```

See also: [`example/actin/README.md`](example/actin)

## Configuration

It uses [toml format](https://github.com/toml-lang/toml) as a config file format.
For more information about the format, see [the spec of toml-v1.0.0](https://toml.io/en/v1.0.0).

### Example file

```toml
file.input           = "example.pdb"
file.output.basename = "output"
file.output.formats  = ["tsv", "json", "ppm", "svg"]
probe.size           = {radius = "1.0nm", angle = 10.0}
resolution.x         = "1.0nm"
resolution.y         = "1.0nm"
resolution.z         = "0.64angstrom"
range.x              = ["0.0nm", "100.0nm"]
range.y              = ["0.0nm", "100.0nm"]
scale_bar.length     = "5.0nm"
noise                = "0.3nm"
```

### Reference

Length unit is angstrom by default.
You can explicitly specify `"nm"`, `"pm"`, and `"angstrom"`.

#### `file` table

- `file.input`: String
  - An input structure. `.pdb` or `.xyz` are available.
  - note that it may fail reading pdb files that does not conform [wwPDB 3.3](http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM).
- `file.output.basename`: String
  - An output filename. An extension will be added.
- `file.output.formats`: Array of Strings
  - An output format. One or more formats can be selected from `["csv", "json", "ppm", "svg"]`.

#### `method`

- `method`: String, `"rigid"` or `"smooth"`.
  - By default, `"rigid"`.
  - `"rigid"` is a collision-detection based method.
  - `"smooth"` is a method introduced in the [paper](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00991)

#### `probe` table

In the `"rigid"` method, the AFM cantilever tip will be modeled as a hemisphere on top of truncated cone.
If you specify `"smooth"` method, the value will be ignored.

- `probe.size.radius`: String or Floating
  - The radius of the probe. 
- `probe.size.angle`: Floating
  - The angle of the cone. The unit is degree, not radian.

#### `resolution` table

- `resolution.x`: String or Floating
  - The resolution in x direction. The same as the pixel width.
- `resolution.y`: String or Floating
  - The resolution in y direction. The same as the pixel height.
- `resolution.z`: String or Floating
  - The resolution in z direction. The output height will be rounded using this.

#### `range` table

- `range.x`: Array of Strings or Floatings
  - The minimum and maximum coordinate of scanning range in x direction.
  - Atoms that exceed this boundary will not be scanned.
- `range.y`: Array of Strings or Floatings
  - The minimum and maximum coordinate of scanning range in y direction.
  - Atoms that exceed this boundary will not be scanned.

#### `stage` table

- `stage.position`: String or Floating
  - The stage position in Z direction.
  - By default, `0`.
- `stage.align`: Boolean
  - If `true`, the position in Z axis of the structure will be aligned to `0`.
  - Otherwise, the position will be kept intact.
  - By default, `true`.

#### `radii` table

You can change the atom radius in this table.
By defualt, only a few number of atoms are supported.
If you got some error like `unknown atom`, specify its radius in this table.

- `radii.atom.[name-of-atom]`: String or Floating
  - Replace `[name-of-atom]` by your atom name.
  - The atom `[name-of-atom]` will have the radius you specified here.
- `radii.residue.[name-of-residue]`: String or Floating
  - Replace `[name-of-residue]` by your residue name.
  - All the atoms in the residue `[name-of-residue]` will have the same radius you specified here.

#### `scale_bar` table

Scale bar option is only available with `.svg` output format.
If you specify other format, the value will be ignored.

- `length`: String or Floating
  - The length of the scale bar.

#### `sigma` and `gamma`

Parameters `sigma` and `gamma` that are used in the `"smooth"` method.
If you specify `"rigid"` method, the value will be ignored.

- `sigma`: String or Floating
  - The value of `sigma` in the formula described in the paper.
- `gamma`: String or Floating
  - The value of `gamma` in the formula described in the paper.

## Installation

First, clone this repository via `git`. Do not download auto-generated zip or tar
file that does not include submodule information.

### Prerequisites

Use your favorite package managers (e.g. `apt`) to install them.

- CMake (> 3.2)
  - to generate Makefile.
- git
  - to download submodules.
- C++14 compliant compiler. tested with ...
  - g++-7 or later
  - clang++-6 or later
  - or others that fully support C++14 ISO standard.

All other dependencies are managed by CMake script and git submodule.

### Building

```sh
$ git clone https://github.com/ToruNiina/afmize.git
$ cd afmize
$ git submodule update --init
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make test # optional
```

After this, you will find the executable at `bin/` directory.

## Citation

Please cite the following.

- [Toru Niina, Yasuhiro Matsunaga, & Shoji Takada (2021). Rigid-Body Fitting to Atomic Force Microscopy Images for Inferring Probe Shape and Biomolecular Structure. bioRxiv. doi:10.1101/2021.02.21.432132](https://www.biorxiv.org/content/10.1101/2021.02.21.432132v1)
- [Toru Niina, Sotaro Fuchigami, & Shoji Takada (2020). Flexible Fitting of Biomolecular Structures to Atomic Force Microscopy Images via Biased Molecular Simulations. JCTC. doi:10.1021/acs.jctc.9b00991](https://doi.org/10.1021/acs.jctc.9b00991)
- [Toru Niina, & Suguru Kato. (2019, August 7). ToruNiina/afmize: version 1.1.0 (Version v1.1.0). Zenodo. doi:10.5281/zenodo.2556444](https://doi.org/10.5281/zenodo.2556444)

## Contact

If you have any question, please feel free to make an [issue](https://github.com/ToruNiina/afmize/issues) to this repository.

Issues are public. Everyone can share information about problems and save time.

If you want to share a sensitive data with the repository owner to solve the problem,
you can e-mail to [me](https://github.com/ToruNiina).

## Licensing Terms

This product is licensed under the terms of the [MIT License](LICENSE).

- Copyright (c) 2018-2021 Toru Niina

All rights reserved.
