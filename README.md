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

Then open it via `Tools ▸ Volume Data ▸ Subbox Particles`.

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

Offsets are measured between the intensity-weighted centres of mass of the
maps, so keep the parent and the monomer at the same contour level.

MIT licensed.
