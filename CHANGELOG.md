# Changelog

## 1.4

* **The force-zeroing now checks the parent map's orientation.** Child maps
  arranged around a filament or a ring share one rotation axis; the bundle now
  measures it and warns, in the log and from the preview, when it is more than
  5 degrees away from the axis the zeroing keeps. Caught on a Dam1 ring whose
  axis sits 34 degrees off the parent box's Z: zeroing X and Y kept the 90 A
  Z component and threw away 54 A that was almost entirely *along* the ring
  axis, displacing every sub-particle sideways. The 1.2 warning did not fire,
  because it only triggers when the zeroing removes more than it keeps.
* **Fixed: preview arrows are scaled from the child map's box** rather than
  from the spread of the child positions. Rotational copies of one subunit can
  sit within a few Angstrom of each other, which gave sub-Angstrom arrows
  hidden inside their own shafts, and made the glyph size change with the
  arrangement for no visible reason.
* New `core.common_child_axis` and `core.force_axis_warning`.

## 1.3

* **Child anchors are now box centres too** — every offset runs from the
  parent's box centre to the child's box centre, so the transform is purely
  geometric and density/contour levels no longer affect the result at all.
  The child's box centre is also where the sub-particle box is centred, which
  makes the child map directly usable as the reference for the subboxed
  refinement. Centre the monomer in its box when cropping.
* Removed `measure_center` (the intensity-weighted COM measurement); nothing
  uses it any more.

## 1.2

* **Fixed: offsets are now anchored at the parent's box centre** (voxel N//2,
  the point RELION/Warp coordinates refer to) instead of the parent's density
  centre of mass. For a parent map whose density is off-centre the old anchor
  shifted every sub-particle by the COM-to-centre difference (108.6 Å on the
  map this was caught with). Child anchors are unchanged (intensity-weighted
  COM).
* **Fixed: Euler recovery at tilt = 180°.** The gimbal-lock branch of
  `matrix_to_euler_zyz` applied the tilt = 0 folding formula in both
  degenerate cases, returning `rot` off by 180° when tilt = 180°.
* The log and the transform preview now **warn when force-zeroing removes most
  of a child's offset** (more than what remains, and > 1 Å). The X/Y zeroing
  default suits filaments like the microtubule demo, but silently collapsed
  the ~160 Å inter-monomer offset of an arbitrarily oriented ATP synthase
  dimer, producing a reconstruction of the parent map mixed in two
  orientations.
* `compute_transforms` detail tuples gained a fourth element (the Ångström
  magnitude removed by zeroing).

## 1.1

* Added per-tomogram deduplication of the generated sub-particles
  (`minDistance`, in Ångström; 0 = off). Within one tomogram, a sub-particle
  closer than the threshold to a sub-particle that was already kept is
  discarded, greedily and in file order. Particles are never compared across
  tomograms: if the STAR file has no recognised tomogram column and no
  `tomoLabel` is given, deduplication is skipped and a message says so.
* Renamed the tool to **Subbox Particles**, matching the capitalisation of the
  other entries in `Tools ▸ Volume Data`.
* A non-positive `pixelSize` is now rejected instead of being silently treated
  as 1.0.

## 1.0

* Initial release: expand a particle STAR file into sub-particles, one per
  child map placed on a parent complex map. RELION and Warp columns, ZYZ Euler
  handling, offset measurement from the maps' centres of mass, transform
  preview, `subbox` command and GUI panel.
