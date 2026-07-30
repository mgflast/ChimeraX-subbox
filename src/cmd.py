"""The ``subbox`` command.

``minDistance`` is in Angstrom (0 = off) and is converted to coordinate units
with ``pixelSize``, like the child offsets.  ``tomoLabel`` overrides which
column identifies the tomogram for deduplication.
"""

from chimerax.core.commands import (
    CmdDesc, register, OpenFileNameArg, SaveFileNameArg,
    ModelArg, ModelsArg, BoolArg, StringArg, FloatArg,
)
from chimerax.core.errors import UserError


def subbox(session, input, output=None, parent=None, children=None,
           forceX=True, forceY=True, forceZ=False, pixelSize=1.0,
           addParentId=True, parentIdLabel='parentID',
           minDistance=0.0, tomoLabel=None):
    from chimerax.map import Volume
    from . import core

    if output is None:
        raise UserError("Missing required argument: output")
    if parent is None:
        raise UserError("Missing required argument: parent (e.g. parent #1)")
    if children is None:
        raise UserError("Missing required argument: children (e.g. children #2-14)")

    if not isinstance(parent, Volume):
        raise UserError("parent #{} is not a volume/map".format(
            getattr(parent, 'id_string', '?')))

    child_vols = core.child_volumes(children, parent)
    if not child_vols:
        raise UserError("No child volumes found in the 'children' specifier.")

    return core.run_subbox(
        session, input, output, parent, child_vols,
        force=(forceX, forceY, forceZ), coord_pixel_size=pixelSize,
        add_parent_id=addParentId, parent_id_label=parentIdLabel,
        min_distance=minDistance, tomo_label=tomoLabel)


def register_command(command_info, logger):
    desc = CmdDesc(
        required=[('input', OpenFileNameArg)],
        keyword=[
            ('output', SaveFileNameArg),
            ('parent', ModelArg),
            ('children', ModelsArg),
            ('forceX', BoolArg),
            ('forceY', BoolArg),
            ('forceZ', BoolArg),
            ('pixelSize', FloatArg),
            ('addParentId', BoolArg),
            ('parentIdLabel', StringArg),
            ('minDistance', FloatArg),
            ('tomoLabel', StringArg),
        ],
        synopsis="Expand a particle STAR file by placing sub-particles per child map",
    )
    register('subbox', desc, subbox, logger=logger)
