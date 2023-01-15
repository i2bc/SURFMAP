from pathlib import Path
import shutil
from typing import Union


class JunkFilePath:
    """Class to handle files or path to remove
    """
    def __init__(self, elements: list=[]) -> None:
        self.to_remove = elements

    def add(self, element: Union[str, list, tuple, Path]):
        if True in [isinstance(element, x) for x in [str, Path]]:
            self.to_remove.append(element)
        elif True in [isinstance(element, x) for x in [list, tuple]]:
            self.to_remove += element

    def empty(self):
        for element in self.to_remove:
            if not element:
                continue
            if Path(element).is_dir():
                try:
                    shutil.rmtree(element)
                except OSError:
                    pass
            elif Path(element).is_file():
                Path(element).unlink()
