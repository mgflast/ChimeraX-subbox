# ChimeraX-Subbox

A ChimeraX bundle that expands a subtomogram-averaging particle STAR file into
**sub-particles**, one per monomer copy you place in the scene.

## What it does

You have a particle STAR file describing a *complex* (e.g. a microtubule
segment: 13/26/39 αβ-tubulin pairs per particle) with a position and
orientation per particle. To average the *monomer* at higher multiplicity:

1. Open the parent/complex map (`#1`).
2. Copy a cropped monomer/reference map and move it onto each sub-position you
   want to average (e.g. one αβ-tubulin on each of the 13 protofilaments,
   `#2`–`#14`).
3. Run Subbox. For each child map it reads:
   * the **rotation** relative to the parent (from the model's scene position),
   * the **offset** relative to the parent (child centre-of-mass − parent
     centre-of-mass, expressed in the parent's local frame).
4. For every particle × every child, it applies that relative transform to the
   parent pose and writes the resulting sub-particle to a new STAR file.
5. Optionally it then deduplicates: within each tomogram, sub-particles that
   land closer together than a distance threshold are collapsed to one.

RELION (`rlnCoordinateX/Y/Z`, `rlnAngleRot/Tilt/Psi`) and Warp
(`wrpCoordinate*`, `wrpAngle*`, with optional `1` suffix) columns are both
auto-detected. Angles use the RELION intrinsic-ZYZ convention. All other
columns and STAR blocks are preserved; an optional `parentID` column records
which parent each sub-particle came from.

No external Python packages are required — rotation math and STAR I/O are
implemented with numpy only, so it runs in a stock ChimeraX.

## Install

Clone the repository, then from the ChimeraX command line:

```
devel install "/path/to/ChimeraX-Subbox"
```

(That builds and installs it in place. Use `devel clean <path>` to remove build
artifacts.) After installing, restart or run `toolshed reload installed`.

## Use — GUI

`Tools ▸ Volume Data ▸ Subbox Particles`

* **Input STAR** / **Output STAR** — browse to the files.
* **Parent map** — model spec of the complex, e.g. `#1`.
* **Child maps** — the placed monomer copies, e.g. `#2-14`
  (the "Current sel" buttons fill these from the current selection).
* **Zero offset component** — X/Y/Z checkboxes drop that component of each
  child's offset. For a microtubule with the reference centred on the tube
  axis, zero X and Y and keep Z (the helical rise); this is the default.
* **Min. particle distance (Å)** — deduplication threshold (0 = off), with an
  optional tomogram-column override next to it. See
  [Deduplication](#deduplication).
* **Preview transforms** logs the parsed offset + angles per child without
  writing anything. **Generate subbox STAR** does the run.

## Use — command

```
subbox  "in.star"  output "out.star"  parent #1  children #2-14 \
        forceX true  forceY true  forceZ false  pixelSize 1.0 \
        minDistance 0
```

| argument        | meaning                                              | default   |
|-----------------|------------------------------------------------------|-----------|
| `input`         | input STAR file (required, positional)               | —         |
| `output`        | output STAR file                                     | required  |
| `parent`        | parent/complex map model                             | required  |
| `children`      | child monomer maps                                   | required  |
| `forceX/Y/Z`    | zero that offset component (parent frame)            | true/true/false |
| `pixelSize`     | Å/px of the star coordinates (see Units below)       | 1.0       |
| `addParentId`   | add a `parentID` column                              | true      |
| `parentIdLabel` | name of that column                                  | `parentID`|
| `minDistance`   | deduplication threshold in Å (0 = off)               | 0         |
| `tomoLabel`     | column identifying the tomogram                      | auto      |

## Deduplication

Overlapping parent particles (a microtubule traced twice, sliding-window
segments, …) produce sub-particles at nearly the same place. `minDistance`
removes them: the output is walked top to bottom and a sub-particle is dropped
if a sub-particle that was already kept **from the same tomogram** lies within
`minDistance` of it. It is deliberately greedy and first-come-first-served —
there is no attempt to choose which of a clashing pair to delete so as to
maximise the final particle count.

* The threshold is in **Ångström**. It is divided by `pixelSize` to reach the
  units of the coordinate columns, exactly like the child offsets — so if your
  coordinates are in tomogram pixels, set `pixelSize` correctly (see
  [Units](#units--pixel-size)) or the threshold will be off by that factor.
* Deduplication is always done **per tomogram**, so two particles at the same
  X/Y/Z in different tomograms are both kept. The tomogram is read from the
  first of `rlnTomoName`, `rlnMicrographName`, `wrpSourceName`, `wrpSourceHash`
  present in the file — `tomoLabel` overrides that choice if your file
  identifies tomograms by some other column.
* If none of those columns is present and no `tomoLabel` was given, particles
  are **not** deduplicated across the whole file. The run continues, the
  sub-particles are written as usual, and the log says:
  `could not deduplicate because the star file has no tomogram label
  (options: …)`. A `tomoLabel` naming a column that is not in the file is an
  error and nothing is written.
* The log reports how many sub-particles were removed.

## Units / pixel size

Rotations are unitless. Translations are not: the child offset is measured from
the maps in the **reference map's physical units** (Ångström × voxel size).
Before it is added to `rlnCoordinate*`/`wrpCoordinate*`, it is divided by
`pixelSize` (Å/px):

* Star coordinates already in the map's units (or map shown at voxel size 1) →
  leave `pixelSize 1.0`.
* Star coordinates in **tomogram pixels** while the map is displayed in Å → set
  `pixelSize` to that tomogram pixel size. The GUI's **From map** button fills
  in the parent map's voxel size, which is correct when the coordinates share
  the reference map's sampling.

## Notes

* Offsets are taken from `measure center` (intensity-weighted centroid above the
  lowest surface level), matching the manual workflow. Keep the same contour
  level on parent and children for a consistent centre.
* If you rotate the parent map, transforms are still computed relative to it;
  the common case (parent left at identity) is unchanged.

## License

MIT — see [LICENSE](LICENSE).
