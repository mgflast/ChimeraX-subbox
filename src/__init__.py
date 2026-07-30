from chimerax.core.toolshed import BundleAPI


class _SubboxAPI(BundleAPI):

    api_version = 1

    @staticmethod
    def register_command(bi, ci, logger):
        from . import cmd
        cmd.register_command(ci, logger)

    @staticmethod
    def start_tool(session, bi, ti):
        from .tool import SubboxTool
        return SubboxTool(session, ti.name)


bundle_api = _SubboxAPI()
