# ChimeraX-Subbox

Subboxing for ChimeraX: turn particles of a complex into particles of its
subunits. Works with RELION and Warp/M star files.

*Video tutorial coming soon.*

## Install

Download the `.whl` from the
[latest release](https://github.com/mgflast/ChimeraX-subbox/releases/latest)
and, in ChimeraX:

```
toolshed install /path/to/ChimeraX_Subbox-1.1-py3-none-any.whl
```

Then `Tools ▸ Volume Data ▸ Subbox Particles`. To install from source instead:
`devel install /path/to/ChimeraX-Subbox`.

## How it works

1. Open your map.
2. Crop out a monomer.
3. Place duplicates of that monomer onto the corresponding positions in the
   parent map.
4. Point at a star file and press go.

Subbox finds the transform of each monomer relative to the parent map, then for
every particle in the input star file writes one output particle per monomer.

Optionally it deduplicates: set a distance threshold D, and of any pair of
particles within D of each other, one is removed. Per tomogram.

## Command

```
subbox "in.star" output "out.star" parent #1 children #2-14
```

| option | | default |
|---|---|---|
| `pixelSize` | Å/px of the star coordinates — set this if they are in tomogram pixels | 1.0 |
| `minDistance` | deduplication threshold in Å, 0 = off | 0 |
| `forceX/Y/Z` | zero that component of each monomer's offset, in the parent's frame | true/true/false |
| `addParentId` | add a `parentID` column recording which particle each subparticle came from | true |
| `tomoLabel` | column identifying the tomogram, if it isn't `rlnTomoName`/`rlnMicrographName`/`wrpSourceName` | auto |

Offsets are measured between the intensity-weighted centres of mass of the
maps, so keep parent and monomer at the same contour level. **Preview
transforms** draws them in the scene without writing anything.

MIT licensed.
