"""
Core logic for the Subbox bundle.

Turns each particle in a STAR file into several "sub-particles", one per child
map that the user has placed onto a parent complex map in the ChimeraX scene.

No dependencies beyond numpy (which ChimeraX always ships).  Rotation math and
STAR I/O are implemented here so the bundle does not need scipy / pandas /
starfile, none of which are part of a stock ChimeraX install.
"""

import re
import numpy as np


# ---------------------------------------------------------------------------
# Rotation helpers -- intrinsic ZYZ Euler angles, degrees.
#
# These reproduce scipy's `Rotation.from_euler('ZYZ', ...)` /
# `.as_euler('ZYZ', ...)` (upper-case = intrinsic) exactly, so results match
# the original prototype script.
# ---------------------------------------------------------------------------

def _Rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _Ry(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def euler_zyz_to_matrix(rot, tilt, psi):
    """Intrinsic ZYZ Euler (degrees) -> 3x3 rotation matrix (Rz*Ry*Rz)."""
    a, b, c = np.radians([rot, tilt, psi])
    return _Rz(a) @ _Ry(b) @ _Rz(c)


def matrix_to_euler_zyz(m):
    """3x3 rotation matrix -> intrinsic ZYZ Euler (rot, tilt, psi) in degrees.

    tilt is returned in [0, 180] to match scipy's symmetric-sequence output.
    """
    tilt = np.arccos(np.clip(m[2, 2], -1.0, 1.0))
    if np.sin(tilt) > 1e-6:
        rot = np.arctan2(m[1, 2], m[0, 2])
        psi = np.arctan2(m[2, 1], -m[2, 0])
    else:
        # Gimbal lock (tilt == 0 or 180): fold rot+psi into rot.
        rot = np.arctan2(m[1, 0], m[0, 0])
        psi = 0.0
    return np.degrees([rot, tilt, psi])


def relion_euler_to_matrix(rot, tilt, psi):
    """RELION particle angles -> rotation matrix.

    Matches the prototype's `R.from_euler('ZYZ', ...).inv()`, i.e. the transpose
    of the plain ZYZ matrix.
    """
    return euler_zyz_to_matrix(rot, tilt, psi).T


def matrix_to_relion_euler(m):
    """Inverse of `relion_euler_to_matrix`; returns (rot, tilt, psi) degrees."""
    return matrix_to_euler_zyz(m.T)


# ---------------------------------------------------------------------------
# Per-particle transform application
# ---------------------------------------------------------------------------

def apply_transforms(mt_pos, mt_euler, transforms):
    """Compute child particle poses from a parent particle pose.

    Parameters
    ----------
    mt_pos : (3,) array           parent coordinate (star units, e.g. pixels)
    mt_euler : (rot, tilt, psi)   parent angles in degrees (RELION convention)
    transforms : list of (rel_pos, child_rot)
        rel_pos   : (3,) offset of the child from the parent, in the parent's
                    local frame (star units)
        child_rot : (3,3) rotation of the child relative to the parent

    Returns
    -------
    list of (abs_pos (3,), abs_euler (rot, tilt, psi) degrees)
    """
    mp = relion_euler_to_matrix(*mt_euler)
    out = []
    for rel_pos, child_rot in transforms:
        abs_pos = np.asarray(mt_pos, float) + mp @ np.asarray(rel_pos, float)
        abs_mat = mp @ np.asarray(child_rot, float)
        abs_euler = matrix_to_relion_euler(abs_mat)
        out.append((abs_pos, abs_euler))
    return out


# ---------------------------------------------------------------------------
# Reading transforms out of the live ChimeraX session
# ---------------------------------------------------------------------------

def child_volumes(models, parent):
    """Keep the Volume models from a selection, dropping the parent and any
    non-volume submodels, ordered by model id."""
    from chimerax.map import Volume
    vols = [m for m in models if isinstance(m, Volume) and m is not parent]
    vols.sort(key=lambda m: m.id)
    return vols


def measure_center(volume):
    """Intensity-weighted centre of mass of a Volume, in scene coordinates.

    Reproduces ChimeraX's ``measure center``: grid points at or above the
    lowest displayed contour level are weighted by their value.  Falls back to
    the mean level if the map has no surface.
    """
    m = volume.full_matrix()
    try:
        levels = [s.level for s in volume.surfaces]
    except Exception:
        levels = []
    level = min(levels) if levels else float(m.mean())

    mask = m >= level
    if not mask.any():
        mask = np.ones_like(m, dtype=bool)

    idx = np.argwhere(mask)                 # rows are (k, j, i) == (z, y, x)
    w = m[mask].astype(np.float64)
    w = np.clip(w, 0.0, None)
    if w.sum() <= 0:
        w = np.ones_like(w)
    mean_kji = (idx * w[:, None]).sum(axis=0) / w.sum()
    ijk = mean_kji[::-1]                     # -> (i, j, k) == (x, y, z)

    xyz_data = volume.data.ijk_to_xyz(ijk)
    scene = volume.scene_position.transform_points(np.array([xyz_data], float))[0]
    return scene


def compute_transforms(parent, children, force=(True, True, False)):
    """Derive (rel_pos, child_rot) for each child relative to the parent.

    Rotation comes from the model's ``scene_position``; translation comes from
    the centre-of-mass difference, expressed in the parent's local frame.  The
    ``force`` triple zeroes the X / Y / Z components of the offset.
    """
    p_place = parent.scene_position
    p_rot = np.asarray(p_place.matrix, float)[:, :3]      # 3x3
    p_center = measure_center(parent)
    fx, fy, fz = force

    transforms = []
    details = []
    for child in children:
        c_place = child.scene_position
        rel_place = p_place.inverse() * c_place
        child_rot = np.asarray(rel_place.matrix, float)[:, :3]

        c_center = measure_center(child)
        offset_scene = np.asarray(c_center, float) - np.asarray(p_center, float)
        rel_pos = p_rot.T @ offset_scene                  # into parent frame
        if fx:
            rel_pos[0] = 0.0
        if fy:
            rel_pos[1] = 0.0
        if fz:
            rel_pos[2] = 0.0

        transforms.append((rel_pos, child_rot))
        details.append((child.id_string, rel_pos.copy(),
                        matrix_to_euler_zyz(child_rot)))
    return transforms, details


def make_preview_bild(parent, children, force=(True, True, False)):
    """Return a BILD graphics description of the child transforms.

    A dot marks each child position (reconstructed from the parsed offset,
    in scene coordinates so it overlays the maps) and red/green/blue arrows
    show its X/Y/Z axes.  A grey dot marks the parent centre.
    """
    p_place = parent.scene_position
    p_rot = np.asarray(p_place.matrix, float)[:, :3]
    p_center = measure_center(parent)
    transforms, _ = compute_transforms(parent, children, force)

    pts, rots = [], []
    for rel_pos, child_rot in transforms:
        pts.append(p_center + p_rot @ rel_pos)
        rots.append(p_rot @ child_rot)
    pts = np.array(pts)

    if len(pts) > 1:
        diag = float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)))
    else:
        diag = 0.0
    alen = diag * 0.12 if diag > 0 else 20.0
    rad = max(alen * 0.06, 1.0)

    axis_colors = [(1.0, 0.25, 0.25), (0.35, 0.85, 0.35), (0.3, 0.45, 1.0)]
    out = [".comment subbox preview",
           ".color 0.6 0.6 0.6",
           ".sphere {:.3f} {:.3f} {:.3f} {:.3f}".format(
               p_center[0], p_center[1], p_center[2], rad * 1.6)]
    for pos, rot in zip(pts, rots):
        out.append(".color 1 1 1")
        out.append(".sphere {:.3f} {:.3f} {:.3f} {:.3f}".format(
            pos[0], pos[1], pos[2], rad))
        for a in range(3):
            tip = pos + rot[:, a] * alen
            out.append(".color {:.3f} {:.3f} {:.3f}".format(*axis_colors[a]))
            out.append(".arrow {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} "
                       "{:.3f} {:.3f}".format(
                           pos[0], pos[1], pos[2], tip[0], tip[1], tip[2],
                           rad * 0.4, rad * 1.1))
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# Minimal STAR reader / writer (loop_ blocks, multi-block)
# ---------------------------------------------------------------------------

class StarBlock:
    def __init__(self, name):
        self.name = name
        self.is_loop = False
        self.labels = []        # column names, without leading underscore
        self.rows = []          # list of list-of-str (loop blocks)
        self.pairs = []         # list of (key, value) for non-loop blocks


def read_star(path):
    """Parse a STAR file into a list of StarBlock objects (order preserved)."""
    with open(path, 'r') as fh:
        lines = fh.read().splitlines()

    blocks = []
    block = None
    i = 0
    n = len(lines)
    while i < n:
        raw = lines[i]
        s = raw.strip()
        if s.startswith('data_'):
            block = StarBlock(s[5:].strip())
            blocks.append(block)
            i += 1
            continue
        if block is None or s == '' or s.startswith('#'):
            i += 1
            continue
        if s == 'loop_':
            block.is_loop = True
            i += 1
            # collect column labels
            while i < n:
                t = lines[i].strip()
                if t.startswith('_'):
                    label = t.split()[0][1:]        # drop '_' and trailing '#n'
                    block.labels.append(label)
                    i += 1
                elif t == '' or t.startswith('#'):
                    i += 1
                else:
                    break
            # collect data rows
            while i < n:
                t = lines[i].strip()
                if t == '' or t.startswith('#'):
                    i += 1
                    continue
                if t.startswith('data_') or t == 'loop_':
                    break
                block.rows.append(t.split())
                i += 1
            continue
        if s.startswith('_'):
            parts = s.split(None, 1)
            key = parts[0][1:]
            val = parts[1] if len(parts) > 1 else ''
            block.pairs.append((key, val))
            i += 1
            continue
        i += 1
    return blocks


def _fmt_num(v):
    return "{:.6f}".format(float(v))


def write_star(path, blocks):
    """Serialize StarBlock objects back to a STAR file."""
    out = []
    for b in blocks:
        out.append("data_{}\n".format(b.name))
        out.append("\n")
        if b.is_loop:
            out.append("loop_\n")
            for k, label in enumerate(b.labels, start=1):
                out.append("_{} #{}\n".format(label, k))
            for row in b.rows:
                out.append("  " + "  ".join(str(x) for x in row) + "\n")
        else:
            for key, val in b.pairs:
                out.append("_{}  {}\n".format(key, val))
        out.append("\n")
    with open(path, 'w') as fh:
        fh.write("".join(out))


# ---------------------------------------------------------------------------
# Column detection (RELION and Warp)
# ---------------------------------------------------------------------------

_COORD_SETS = [
    ('rlnCoordinateX', 'rlnCoordinateY', 'rlnCoordinateZ'),
    ('wrpCoordinateX', 'wrpCoordinateY', 'wrpCoordinateZ'),
    ('wrpCoordinateX1', 'wrpCoordinateY1', 'wrpCoordinateZ1'),
]
_ANGLE_SETS = [
    ('rlnAngleRot', 'rlnAngleTilt', 'rlnAnglePsi'),
    ('wrpAngleRot', 'wrpAngleTilt', 'wrpAnglePsi'),
    ('wrpAngleRot1', 'wrpAngleTilt1', 'wrpAnglePsi1'),
]


_TOMO_LABELS = [
    'rlnTomoName',
    'rlnMicrographName',
    'wrpSourceName',
    'wrpSourceHash',
]


def _pick_set(labels, candidates):
    label_set = set(labels)
    for cand in candidates:
        if all(c in label_set for c in cand):
            return cand
    return None


def find_tomo_column(labels, preferred=None):
    """Return the label that identifies the tomogram, or None.

    `preferred` wins if it is present; otherwise the first of `_TOMO_LABELS`
    that occurs in `labels` is used.
    """
    if preferred:
        return preferred if preferred in labels else None
    for c in _TOMO_LABELS:
        if c in labels:
            return c
    return None


def find_particle_block(blocks):
    """Return (block, coord_cols, angle_cols) for the block holding particles."""
    for b in blocks:
        if not b.is_loop:
            continue
        coord = _pick_set(b.labels, _COORD_SETS)
        angle = _pick_set(b.labels, _ANGLE_SETS)
        if coord and angle:
            return b, coord, angle
    return None, None, None


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------

def deduplicate_rows(rows, coord_idx, min_distance, tomo_idx):
    """Greedily drop particles that sit too close to an already-kept particle.

    Rows are visited in order; a row is kept unless a previously kept row from
    the *same* tomogram lies within `min_distance` of it, in which case it is
    discarded (first-come-first-served, no attempt to maximise the final
    count).  Distances are in the units of the coordinate columns.

    Parameters
    ----------
    rows : list of list-of-str      particle rows
    coord_idx : (ix, iy, iz)        column indices of the coordinates
    min_distance : float            threshold; <= 0 disables and returns as-is
    tomo_idx : int                  column index grouping rows per tomogram;
                                    required -- particles are only ever
                                    compared within one tomogram

    Returns
    -------
    (kept_rows, n_removed)
    """
    d = float(min_distance)
    if d <= 0.0:
        return list(rows), 0

    ix, iy, iz = coord_idx
    d2 = d * d
    # One uniform grid per tomogram, cell size = d, so only the 27 neighbouring
    # cells can hold a particle within d.
    grids = {}
    kept = []
    neighbourhood = [(a, b, c)
                     for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]

    for row in rows:
        x = float(row[ix]); y = float(row[iy]); z = float(row[iz])
        grid = grids.setdefault(row[tomo_idx], {})
        ci = int(np.floor(x / d)); cj = int(np.floor(y / d)); ck = int(np.floor(z / d))

        clash = False
        for a, b, c in neighbourhood:
            for px, py, pz in grid.get((ci + a, cj + b, ck + c), ()):
                if (x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2 < d2:
                    clash = True
                    break
            if clash:
                break
        if clash:
            continue

        grid.setdefault((ci, cj, ck), []).append((x, y, z))
        kept.append(row)

    return kept, len(rows) - len(kept)


# ---------------------------------------------------------------------------
# Top-level driver
# ---------------------------------------------------------------------------

def run_subbox(session, star_in, star_out, parent, children,
               force=(True, True, False), coord_pixel_size=1.0,
               add_parent_id=True, parent_id_label='parentID',
               min_distance=0.0, tomo_label=None):
    """Read `star_in`, expand every particle by the child transforms, write
    `star_out`.  Returns a short human-readable summary string.

    `coord_pixel_size` (Angstrom/pixel of the star coordinates) converts the
    map-derived offsets, which are in the reference map's physical units, into
    the units of rlnCoordinate*/wrpCoordinate*.  Use 1.0 when they already
    match; use the star's pixel size when the coordinates are in tomogram
    pixels but the map is displayed in Angstrom.

    `min_distance` (Angstrom, 0 = off) removes duplicate sub-particles: within
    each tomogram, a particle closer than this to an already-kept particle is
    discarded.  It is divided by `coord_pixel_size` to get the threshold in
    coordinate units, exactly like the child offsets.  `tomo_label` names the
    column that identifies the tomogram; by default the first recognised one is
    used.
    """
    from chimerax.core.errors import UserError

    if not children:
        raise UserError("No child maps selected.")
    if float(coord_pixel_size) <= 0.0:
        raise UserError("coordPixelSize must be positive, got {}".format(
            coord_pixel_size))

    transforms, details = compute_transforms(parent, children, force)

    if coord_pixel_size != 1.0:
        transforms = [(rel_pos / float(coord_pixel_size), child_rot)
                      for rel_pos, child_rot in transforms]

    log = session.logger
    log.info("Subbox: parent #{}, {} child map(s), force(X,Y,Z)={}, "
             "coordPixelSize={} A/px".format(
                 parent.id_string, len(children), force, coord_pixel_size))
    for cid, rel_pos, euler in details:
        log.info("  child #{}: offset=({:.2f}, {:.2f}, {:.2f})  "
                 "rot/tilt/psi=({:.2f}, {:.2f}, {:.2f})".format(
                     cid, rel_pos[0], rel_pos[1], rel_pos[2],
                     euler[0], euler[1], euler[2]))

    blocks = read_star(star_in)
    pblock, coord_cols, angle_cols = find_particle_block(blocks)
    if pblock is None:
        raise UserError(
            "Could not find a particle block with recognised coordinate/angle "
            "columns (rln* or wrp*) in {}".format(star_in))

    labels = pblock.labels
    cx, cy, cz = (labels.index(c) for c in coord_cols)
    ar, at, ap = (labels.index(a) for a in angle_cols)

    if add_parent_id and parent_id_label not in labels:
        labels.append(parent_id_label)
        pid_col = len(labels) - 1
    else:
        pid_col = labels.index(parent_id_label) if add_parent_id else None

    n_parents = len(pblock.rows)
    new_rows = []
    for pidx, row in enumerate(pblock.rows):
        mt_pos = np.array([row[cx], row[cy], row[cz]], float)
        mt_euler = (float(row[ar]), float(row[at]), float(row[ap]))
        for abs_pos, abs_euler in apply_transforms(mt_pos, mt_euler, transforms):
            nr = list(row)
            if pid_col is not None and pid_col == len(nr):
                nr.append('')                 # room for the new parentID column
            nr[cx] = _fmt_num(abs_pos[0])
            nr[cy] = _fmt_num(abs_pos[1])
            nr[cz] = _fmt_num(abs_pos[2])
            nr[ar] = _fmt_num(abs_euler[0])
            nr[at] = _fmt_num(abs_euler[1])
            nr[ap] = _fmt_num(abs_euler[2])
            if pid_col is not None:
                nr[pid_col] = str(pidx)
            new_rows.append(nr)

    n_expanded = len(new_rows)
    dedup_note = ""
    if min_distance and float(min_distance) > 0.0:
        # The threshold is given in Angstrom; the coordinates are in star
        # units, so convert with the same pixel size used for the offsets.
        px = float(coord_pixel_size)
        threshold = float(min_distance) / px

        tomo_col = find_tomo_column(labels, tomo_label)
        if tomo_label and tomo_col is None:
            raise UserError(
                "Requested tomogram column '{}' is not in the star file. "
                "Available columns: {}".format(tomo_label, ", ".join(labels)))
        if tomo_col is None:
            # Deduplication is only ever meaningful within one tomogram, so
            # without a tomogram column we skip it rather than compare
            # coordinates from different tomograms.
            log.warning(
                "Subbox: could not deduplicate because the star file has no "
                "tomogram label (options: {}). Use tomoLabel to name the "
                "column that identifies the tomogram.".format(
                    ", ".join(_TOMO_LABELS)))
            dedup_note = ", deduplication skipped (no tomogram label)"
        else:
            log.info("Subbox: deduplicating per '{}' at {} A = {:.4g} "
                     "coordinate units ({} A/px)".format(
                         tomo_col, min_distance, threshold, px))
            new_rows, n_removed = deduplicate_rows(
                new_rows, (cx, cy, cz), threshold, labels.index(tomo_col))
            dedup_note = ", {} removed by {} A deduplication".format(
                n_removed, min_distance)

    pblock.rows = new_rows

    write_star(star_out, blocks)

    summary = ("Subbox: {} parent particles x {} children = {} sub-particles{} "
               "-> {} written to {}".format(
                   n_parents, len(children), n_expanded, dedup_note,
                   len(new_rows), star_out))
    log.status(summary, log=True)
    return summary
