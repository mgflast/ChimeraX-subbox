"""GUI panel for the Subbox bundle."""

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.ui import MainToolWindow


class SubboxTool(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = False
    help = None

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "Subbox Particles"
        self.tool_window = MainToolWindow(self)
        self._build_ui()
        self.tool_window.manage(placement="side")

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        from Qt.QtWidgets import (
            QWidget, QGridLayout, QVBoxLayout, QHBoxLayout, QLabel,
            QLineEdit, QPushButton, QCheckBox, QGroupBox, QFileDialog,
        )

        parent = self.tool_window.ui_area
        outer = QVBoxLayout()
        parent.setLayout(outer)

        grid = QGridLayout()
        outer.addLayout(grid)
        row = 0

        def browse(line, save):
            def _do():
                if save:
                    fn, _ = QFileDialog.getSaveFileName(
                        parent, "Output STAR file", line.text(),
                        "STAR files (*.star);;All files (*)")
                else:
                    fn, _ = QFileDialog.getOpenFileName(
                        parent, "Input STAR file", line.text(),
                        "STAR files (*.star);;All files (*)")
                if fn:
                    line.setText(fn)
            return _do

        # Input star
        grid.addWidget(QLabel("Input STAR:"), row, 0)
        self.in_edit = QLineEdit()
        grid.addWidget(self.in_edit, row, 1)
        b = QPushButton("Browse…")
        b.clicked.connect(browse(self.in_edit, False))
        grid.addWidget(b, row, 2)
        row += 1

        # Output star
        grid.addWidget(QLabel("Output STAR:"), row, 0)
        self.out_edit = QLineEdit()
        grid.addWidget(self.out_edit, row, 1)
        b = QPushButton("Browse…")
        b.clicked.connect(browse(self.out_edit, True))
        grid.addWidget(b, row, 2)
        row += 1

        # Parent map
        grid.addWidget(QLabel("Parent map:"), row, 0)
        self.parent_edit = QLineEdit("#1")
        self.parent_edit.setToolTip("Model spec of the complex/parent map, e.g. #1")
        grid.addWidget(self.parent_edit, row, 1)
        b = QPushButton("Current sel")
        b.clicked.connect(lambda: self._fill_from_selection(self.parent_edit))
        grid.addWidget(b, row, 2)
        row += 1

        # Child maps
        grid.addWidget(QLabel("Child maps:"), row, 0)
        self.children_edit = QLineEdit("#2-14")
        self.children_edit.setToolTip(
            "Model spec of the placed monomer copies, e.g. #2-14")
        grid.addWidget(self.children_edit, row, 1)
        b = QPushButton("Current sel")
        b.clicked.connect(lambda: self._fill_from_selection(self.children_edit))
        grid.addWidget(b, row, 2)
        row += 1

        # Coordinate pixel size
        grid.addWidget(QLabel("Coord pixel size (Å/px):"), row, 0)
        self.px_edit = QLineEdit("1.0")
        self.px_edit.setToolTip(
            "Ångström/pixel of the star's rlnCoordinate values. The map-derived "
            "offset (in Å) is divided by this before it is added to the star "
            "coordinates. Use 1.0 if the star coordinates are already in the "
            "reference map's units; use the tomogram pixel size if they are in "
            "tomogram pixels.")
        grid.addWidget(self.px_edit, row, 1)
        b = QPushButton("From map")
        b.setToolTip("Fill from the parent map's voxel size")
        b.clicked.connect(self._fill_pixel_size_from_map)
        grid.addWidget(b, row, 2)
        row += 1

        # Deduplication
        grid.addWidget(QLabel("Min. particle distance:"), row, 0)
        self.dedup_edit = QLineEdit("0")
        self.dedup_edit.setToolTip(
            "Deduplication threshold, in the same units as the star's "
            "coordinates (divide an Ångström distance by the coordinate pixel "
            "size above). Within each tomogram, a sub-particle closer than "
            "this to a sub-particle that was already kept is discarded. "
            "0 disables deduplication.")
        grid.addWidget(self.dedup_edit, row, 1)
        self.tomo_edit = QLineEdit()
        self.tomo_edit.setPlaceholderText("tomogram column (auto)")
        self.tomo_edit.setToolTip(
            "Column that identifies the tomogram, e.g. rlnTomoName, "
            "rlnMicrographName or wrpSourceName. Leave empty to auto-detect.")
        grid.addWidget(self.tomo_edit, row, 2)
        row += 1

        # Force-zero options
        force_box = QGroupBox("Zero offset component (parent frame)")
        fl = QHBoxLayout()
        self.fx = QCheckBox("X"); self.fx.setChecked(True)
        self.fy = QCheckBox("Y"); self.fy.setChecked(True)
        self.fz = QCheckBox("Z"); self.fz.setChecked(False)
        fl.addWidget(self.fx); fl.addWidget(self.fy); fl.addWidget(self.fz)
        fl.addStretch(1)
        force_box.setLayout(fl)
        outer.addWidget(force_box)

        # Run / preview
        btn_row = QHBoxLayout()
        preview = QPushButton("Preview transforms")
        preview.clicked.connect(self._preview)
        run_btn = QPushButton("Generate subbox STAR")
        run_btn.clicked.connect(self._run)
        btn_row.addWidget(preview)
        btn_row.addWidget(run_btn)
        outer.addLayout(btn_row)

        self.status = QLabel("")
        self.status.setWordWrap(True)
        outer.addWidget(self.status)

    # -------------------------------------------------------------- helpers
    def _fill_from_selection(self, line):
        sel = [m for m in self.session.selection.models()]
        if not sel:
            self.status.setText("Nothing selected.")
            return
        line.setText(" ".join("#" + m.id_string for m in sel))

    def _build_command(self):
        in_star = self.in_edit.text().strip()
        out_star = self.out_edit.text().strip()
        parent = self.parent_edit.text().strip()
        children = self.children_edit.text().strip()
        if not (in_star and out_star and parent and children):
            self.status.setText(
                "Please fill in input, output, parent and child maps.")
            return None
        try:
            px = float(self.px_edit.text().strip())
        except ValueError:
            self.status.setText("Coordinate pixel size must be a number.")
            return None
        try:
            dedup = float(self.dedup_edit.text().strip() or "0")
        except ValueError:
            self.status.setText("Minimum particle distance must be a number.")
            return None
        cmd = (
            'subbox "{inp}" output "{out}" parent {p} children {c} '
            'forceX {fx} forceY {fy} forceZ {fz} pixelSize {px} '
            'minDistance {dd}'
        ).format(
            inp=in_star, out=out_star, p=parent, c=children,
            fx=str(self.fx.isChecked()).lower(),
            fy=str(self.fy.isChecked()).lower(),
            fz=str(self.fz.isChecked()).lower(),
            px=px, dd=dedup,
        )
        tomo = self.tomo_edit.text().strip()
        if tomo:
            cmd += ' tomoLabel "{}"'.format(tomo)
        return cmd

    def _fill_pixel_size_from_map(self):
        from chimerax.map import Volume
        from chimerax.core.commands import ModelArg
        try:
            m = ModelArg.parse(self.parent_edit.text().strip(), self.session)[0]
        except Exception as e:
            self.status.setText("Could not parse parent map: {}".format(e))
            return
        if not isinstance(m, Volume):
            self.status.setText("Parent is not a volume/map.")
            return
        step = m.data.step  # (x, y, z) voxel size in physical units
        self.px_edit.setText("{:g}".format(step[0]))
        self.status.setText(
            "Set coord pixel size to the parent map voxel size ({:g}).".format(
                step[0]))

    def _preview(self):
        # Report the transforms and draw them in the scene (no files touched).
        import os
        import tempfile
        from chimerax.map import Volume
        from chimerax.core.commands import ModelArg, ModelsArg
        from . import core
        try:
            parent = ModelArg.parse(self.parent_edit.text().strip(), self.session)[0]
            children = ModelsArg.parse(self.children_edit.text().strip(), self.session)[0]
        except Exception as e:
            self.status.setText("Could not parse model specs: {}".format(e))
            return
        if not isinstance(parent, Volume):
            self.status.setText("Parent is not a volume/map.")
            return
        child_vols = core.child_volumes(children, parent)
        if not child_vols:
            self.status.setText("No child volumes found.")
            return

        force = (self.fx.isChecked(), self.fy.isChecked(), self.fz.isChecked())
        _, details = core.compute_transforms(parent, child_vols, force)
        lines = ["Preview of {} child transform(s):".format(len(details))]
        for cid, rel_pos, euler in details:
            lines.append(
                "  #{}: offset=({:.2f}, {:.2f}, {:.2f})  "
                "rot/tilt/psi=({:.2f}, {:.2f}, {:.2f})".format(
                    cid, rel_pos[0], rel_pos[1], rel_pos[2],
                    euler[0], euler[1], euler[2]))
        self.session.logger.info("\n".join(lines))

        # Draw dots + XYZ axes as a BILD model, replacing any previous preview.
        try:
            prev = getattr(self, "_preview_model", None)
            if prev is not None and not prev.deleted:
                run(self.session, "close #!{}".format(prev.id_string))
        except Exception:
            pass
        try:
            bild = core.make_preview_bild(parent, child_vols, force)
            path = os.path.join(tempfile.gettempdir(), "subbox_preview.bild")
            with open(path, "w") as fh:
                fh.write(bild)
            opened = run(self.session, 'open "{}"'.format(path))
            model = opened[0] if isinstance(opened, (list, tuple)) else opened
            if model is not None:
                model.name = "subbox preview"
                self._preview_model = model
        except Exception as e:
            self.session.logger.warning("Could not draw preview: {}".format(e))

        self.status.setText(
            "Previewed {} transforms (axes drawn; details in log).".format(
                len(details)))

    def _run(self):
        cmd = self._build_command()
        if cmd is None:
            return
        try:
            result = run(self.session, cmd)
            if result:
                self.status.setText(str(result))
        except Exception as e:
            self.status.setText("Error: {}".format(e))
            raise
