
import importlib.resources as pkg_resources
from pathlib import Path


PATH_TO_SCRIPTS = Path(pkg_resources.path("surfmap.scripts", 'compute_shell.sh')).parent.resolve()
