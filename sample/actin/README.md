# Actin-filament fitting

This example file includes the following.

- input files to generate a pseudo-AFM image from a pdb file
- input files to fit molecular structure into a pseudo-AFM image

This example considers you have already built the afmize.
Please make sure that you have `afmize` and `simulator` in the `afmize/bin/` directory.

## pseudo-AFM image generation procedure

Here, you will use `gen_image.toml` and `actin_filament_reference.pdb`.

The `gen_image.toml` file has the following values.

```toml
file.input           = "actin_filament_reference.pdb"   # structure file used to generate an AFM image
file.output.basename = "actin_filament_2nmpx"           # output file name
file.output.formats  = ["tsv", "svg"]                   # output file format
probe.size           = {radius = "2.0nm", angle = 10.0} # probe shape used to generate an AFM image
resolution.x         = "2.0nm"                          # X resolution (width of a pixel in X direction)
resolution.y         = "2.0nm"                          # Y resolution (width of a pixel in Y direction)
resolution.z         = "0.64angstrom"                   # Z resolution (height discritization)
range.x              = ["0.0nm", "80.0nm"]              # region to be used to generate an AFM image
range.y              = ["0.0nm", "80.0nm"]              # (this value corresponds to the coordinate in the pdb file)
scale_bar.length     = "5.0nm"                          # scale bar written in the image file.
stage.align          = true                             # a flag to align the input structure on to Z == 0.
stage.position       = 0.0                              # stage position (if you don't have any special reason, set this to 0)
noise                = "0.3nm"                          # noise intensity.
```

When you pass this file to `bin/afmize`, you will get `actin_filament_2nmpx.tsv` and `actin_filament_2nmpx.svg`.
It would take some seconds.

```console
$ pwd
/path/to/afmize/sample/actin
$ ls
6BNO_A15.pdb actin_filament_reference.pdb  gen_image.toml  README.md  scan.toml
$ ../../bin/afmize gen_image.toml
```

`.tsv` file contains the height at each pixel (from (x0, y0), (x1, y0), ..., (xN, y0), (x0, y1), ... (xN, yN)) in angstrom.
`.svg` file represents the AFM image.
The `.svg` file is for quick look, and using some visualization software is recommended to visualize the result.

## scanning-based structure fitting procedure

Here, you will use `scan.toml`, `6BNO_A15.pdb`, and `actin_filament_2nmpx.tsv` that is generated in the previous step.
You can find the description about input parameters as comments in the `scan.toml`.

When you pass `scan.toml` to `bin/simulator`, you will get `actin_filament_fitting_2.0nm_20.0deg_x36_CosineSimilarity_2nmpx_A15.*` files.
It will take some tens of minutes. It is because we are using large pdb file (~80k atoms) and using normal orientation resolution (Δθ = 10°).
If you use small pdb file or a large Δθ, the time required for scanning decreases.

```console
$ pwd
/path/to/afmize/sample/actin
$ ../../bin/simulator scan.toml
```
