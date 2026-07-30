# Changelog

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
