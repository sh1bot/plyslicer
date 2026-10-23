# plyslicer requirements

## Goal

A command line tool that accepts an STL model and divides it into thinner,
printable laminations. Print the pieces with layer lines in different axes,
then glue them together like plywood.

## User requirements

- Identify the largest, flattest, thinnest region of the model automatically.
- Slice along a surface following that region's shape; support curved geometry,
  rather than assuming every input can use a flat plane.
- Produce two pieces by default and allow more than two laminations.
- Preserve the original assembled shape, subject to explicitly chosen manufacturing
  clearances and alignment features.
- Ensure that the pieces can actually be pushed together. Avoid mating geometry
  that traps or interlocks the pieces during assembly.
- Optionally add alignment notches or fiducials to hold alignment while glue sets.
- Export separate printable parts.

## Geometry interpretation to confirm

A peninsula in a surface outline does not necessarily prevent assembly. Undercuts
relative to the insertion direction do. A mating surface that is single-valued
along a common insertion direction permits straight insertion for that interface;
the entire pieces and any alignment features must also be checked for collisions.

The intended meaning of following a curve needs confirmation: a curved mid-surface
that follows the model's thin region, or planar slices selected automatically to
follow its broad orientation. Curved mid-surfaces are the more general target.

## Decisions still open

- How to rank regions when area, flatness, and thickness favor different regions.
- Whether an unsuitable model should be rejected or require user selection of a region.
- Equal thickness versus user-specified thickness for individual laminations.
- Alignment geometry, clearances, minimum wall thickness, and mesh units (STL has
  no unit metadata).
- Export orientation guidance so different print-layer axes are actually achieved.

## Acceptance criteria

- A thin flat solid can be split into the requested number of closed printable meshes.
- A suitable curved thin solid can be laminated while retaining an achievable
  assembly path.
- Unsupported geometry or uncertain assembly feasibility produces a clear diagnostic,
  rather than a claim that the pieces will fit.
- Alignment features preserve sufficient wall thickness and permit assembly.
- Report chosen region, insertion direction, thickness estimates, output units,
  and any manufacturing assumptions.

## Supplied example: shoehorn

The user supplied `Schuhlöffel.3mf` as a representative input. Support STL as
originally requested; also consider 3MF input so this example can be used directly.
The archive contains a millimeter model with one mesh of 57,577 vertices and
115,154 triangles, plus a component object. Raw mesh bounds span approximately
19.22 × 80.00 × 10.61 mm before component/build transforms; these bounds are not
a measurement of wall thickness or final assembled dimensions.

The embedded thumbnail shows a curved spoon-like body and a handle. The intended
laminations should follow the broad thin shell through its thickness, rather than
cutting the shoehorn into short segments along its length. This is a working
interpretation to confirm with the user. Neither mesh validity, local thickness,
nor assembly feasibility has yet been established for this example.

## Initial implementation

The CLI implements STL/3MF loading, whole-solid broad-direction selection, curved
height-field interfaces, multiple laminations, removable alignment-rod bores,
closed/connected output validation, and optional alternating print orientations.
See README.md for commands and limitations. General local-region segmentation,
normal-distance thickness control, glue gaps, and captive alignment keys remain
future work. The original exterior may have overhangs: it is the mating interfaces
that must avoid undercuts along the insertion direction.
