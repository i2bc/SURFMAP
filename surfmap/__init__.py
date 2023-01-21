
import importlib.resources as pkg_resources
from pathlib import Path

with pkg_resources.path("surfmap.scripts", 'compute_shell.sh') as resource:
    PATH_TO_SCRIPTS = Path(resource).parent.resolve()
