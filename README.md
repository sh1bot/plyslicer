# plyslicer

Split a thin **STL or 3MF solid into curved laminations** for separate printing
and gluing. Two layers by default; `--layers` supports 2–32. Optional bores accept
removable alignment rods. This is an initial implementation for a single connected
thin solid, including curved shells such as a shoehorn.

## Install

Requires Python 3.11 or newer. From the checkout:

```sh
uv venv .venv
uv pip install --python .venv/bin/python -e '.[test]'
.venv/bin/pytest -q
```

If your home directory is read-only, prefix the `uv` commands with
`UV_CACHE_DIR=/tmp/plyslicer-uv`. Use the existing checkout in a cloud task;
creating a worktree is unnecessary.

## Use

```sh
.venv/bin/plyslicer model.stl -o output/two-layers
.venv/bin/plyslicer model.3mf -o output/three-layers --layers 3 --pins 2
.venv/bin/plyslicer model.3mf -o output/with-print-orientations --print-oriented
.venv/bin/plyslicer model.stl -o output/fine --pitch 0.4 --direction 0,0,1
```

The output directory must be new. Files are `layer-01.stl`, `layer-02.stl`, etc.,
plus `report.json`. STL coordinates are assumed **millimeters**; 3MF units and
component/build transforms are applied. Outputs retain assembly coordinates.
The tool does not modify the input.

For the supplied example, replace `model.3mf` with the path to `Schuhlöffel.3mf`.
The uploaded file is not bundled with the source.

### Assembly and printing

The report records an insertion direction and a local coordinate frame. Assemble
in increasing layer order, lowering the next piece along the **negative** reported
direction. Each glue interface is a continuous piecewise-linear height field:
each line parallel to that direction meets the interface at most once. Translating
an upper part outward along the direction therefore separates it from lower
parts. A peninsula in the silhouette does not in itself prevent assembly.

`--pins 2 --pin-diameter 2 --clearance 0.15` cuts two straight through-bores for
2 mm rods, with 0.15 mm **radial** clearance (2.30 mm bore diameter). Bore centers
are in local XY coordinates in the report. Insert removable rods to register the
stack while the glue sets, then remove them. Bores intentionally remove exterior
material. This version does not generate captive notches, blind pin sockets, or
printed dowels. Trial-print a fit coupon for your printer and adhesive.

Import the layers separately into your slicer. Orient successive pieces so the
**build directions** differ, for example rotate one piece 90° about its long axis
relative to its neighbor. Merely rotating it on the bed changes infill direction,
not the layer plane. Curved pieces may need supports or brims; inspect the preview.
`--print-oriented` also exports a `print/` subdirectory with alternating build
directions 90° apart and each mesh resting at Z=0. Reversible assembly-to-print
transforms are saved in the report. These are suggestions, not support-minimized
orientations: inspect them before printing. The CLI does not generate G-code.

### Geometry and controls

The automatic direction is chosen from PCA directions and world axes using
projected surface area: it favors the broad thin orientation. Ray sampling finds
the outer lower and upper shell envelopes. Equal fractions between these envelopes
define the curved cuts; boolean intersections preserve the original exterior.
The interfaces share their triangulation and do not cross. Source overhangs are
reported; they do not imply that the glue interfaces interlock.

If the first direction produces disconnected layers, automatic selection tries
other axes and tilted directions. It keeps the first split that passes all output
checks, recording rejected attempts in `report.json`. It never drops disconnected
fragments to force a successful export. An explicit `--direction` is respected
without fallback. This is a finite search, so failure does not prove that every
possible cutting direction is unsuitable.

- `--pitch 0.8`: maximum grid spacing in mm. Smaller values resolve more curvature,
  at greater time and memory cost. Compare outputs at half the pitch for convergence.
- `--direction X,Y,Z`: explicitly choose an insertion direction in input coordinates.
- `--min-thickness 0.4`: reject when the sampled fifth-percentile **axial** thickness
  per layer is below this value. It is a bulk guard, not a minimum wall guarantee.
  Rounded edges can taper to zero, and surface-normal thickness can be smaller.
- `--pins 0`: optional alignment bores; their size and clearance are in mm.

Input must be closed, consistently wound, positive-volume, and a single solid.
Outputs must each remain closed and connected; disconnected laminations are rejected.
Volume conservation is checked before optional bores are cut. Interface assembly
follows from the height-field construction; input projection diagnostics use
floating-point polygon unions with the tolerance included in the report.

This is **not general local-region segmentation**: it laminates the whole solid
and estimates its broad thin orientation. It does not independently rank and split
individual thin regions in a complex bulky object. Cavities and overhangs can make
envelope-based thickness estimates misleading; inspect report warnings and outputs.
No glue gap, adhesive shrinkage, flex, load strength, or printer accuracy is modeled.
Do not treat a successful export as a mechanical fit or strength certification.

## Development

```sh
.venv/bin/pytest -q
```

Tests exercise flat and curved splits, multiple layer counts, volume conservation,
straight separation, source overhangs, alignment bores, invalid parameters, invalid
meshes, CLI exports, and overwrite protection. Heavy boolean work uses Manifold.
