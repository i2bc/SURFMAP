import logging
from pathlib import Path
from typing import Union


def set_root_logger(outpath: Union[str, Path]=".", filename: Union[str, Path]="surfmap.log", level=logging.INFO, *args, **kwargs):
    """ Initialize a root_logger from logging

    Args:
        outpath (Union[str, Path]): Path name of the output directory where a log file will be generated
        filename (Union[str, Path]): Name of the log file
    """
    # set output for log filename
    log_filename = Path(outpath) / Path(filename).name

    # create main root_logger
    root_logger = logging.getLogger("surfmap")
    root_logger.setLevel(logging.DEBUG)

    # create a console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)

    # create a file handler
    file_handler = logging.FileHandler(filename=log_filename, mode="w+")
    file_handler.setLevel(level)

    # create formatter
    basic_formatter = logging.Formatter('%(message)s')
    advanced_formatter = logging.Formatter('%(asctime)s — %(name)s — %(levelname)s — %(funcName)-32s %(message)s')

    # add formatter to handlers
    console_handler.setFormatter(basic_formatter) if level >= logging.INFO else console_handler.setFormatter(advanced_formatter)
    file_handler.setFormatter(advanced_formatter)

    # add console and file handlers to root_logger
    root_logger.addHandler(console_handler)
    root_logger.addHandler(file_handler)


def get_logger(name: str):
    logger = logging.getLogger(name=name)
    logger.addHandler(logging.NullHandler())

    return logger
