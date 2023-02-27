import inspect
import logging
from pathlib import Path
from typing import Union


class IndentFormatter(logging.Formatter):
    def __init__( self, fmt=None, datefmt=None ):
        logging.Formatter.__init__(self, fmt, datefmt)
        self.baseline = len(inspect.stack())
    def format( self, rec ):
        stack = inspect.stack()
        rec.indent = '.'*(len(stack)-self.baseline)*3
        rec.function = stack[8][3]
        out = logging.Formatter.format(self, rec)
        del rec.indent; del rec.function
        return out
    

def addLoggingLevel(name: str, value: int, method: str=None):
    """
    Comprehensively adds a new logging level to the `logging` module and the
    currently configured logging class.

    `name` becomes an attribute of the `logging` module with the value
    `value`. `method` becomes a convenience method for both `logging`
    itself and the class returned by `logging.getLoggerClass()` (usually just
    `logging.Logger`). If `method` is not specified, `name.lower()` is
    used.

    To avoid accidental clobberings of existing attributes, this method will
    raise an `AttributeError` if the level name is already an attribute of the
    `logging` module or if the method name is already present

    Example
    -------
    >>> addLoggingLevel('TRACE', logging.DEBUG - 5)
    >>> logging.getLogger(__name__).setLevel("TRACE")
    >>> logging.getLogger(__name__).trace('that worked')
    >>> logging.trace('so did this')
    >>> logging.TRACE
    5

    """
    if not method:
        method = name.lower()

    if hasattr(logging, name):
       raise AttributeError('{} already defined in logging module'.format(name))
    if hasattr(logging, method):
       raise AttributeError('{} already defined in logging module'.format(method))
    if hasattr(logging.getLoggerClass(), method):
       raise AttributeError('{} already defined in logger class'.format(method))

    # This method was inspired by the answers to Stack Overflow post
    # http://stackoverflow.com/q/2183233/2988730, especially
    # http://stackoverflow.com/a/13638084/2988730
    def logForLevel(self, message, *args, **kwargs):
        if self.isEnabledFor(value):
            self._log(value, message, args, **kwargs)
    def logToRoot(message, *args, **kwargs):
        logging.log(value, message, *args, **kwargs)

    logging.addLevelName(value, name)
    setattr(logging, name, value)
    setattr(logging.getLoggerClass(), method, logForLevel)
    setattr(logging, method, logToRoot)


def set_root_logger(outpath: Union[str, Path]=".", filename: Union[str, Path]="surfmap.log", level=logging.INFO, *args, **kwargs):
    """ Initialize a root_logger from logging

    Args:
        outpath (Union[str, Path]): Path name of the output directory where a log file will be generated
        filename (Union[str, Path]): Name of the log file
    """
    # add a new logging level
    addLoggingLevel('TRACE', logging.DEBUG - 5)

    # create main root_logger
    root_logger = logging.getLogger("surfmap")
    root_logger.setLevel(logging.TRACE)

    # create a console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)

    # create a file handler
    log_filename = Path(outpath) / Path(filename).name
    file_handler = logging.FileHandler(filename=log_filename, mode="w+")
    file_handler.setLevel(logging.DEBUG) if level != logging.TRACE else file_handler.setLevel(logging.TRACE)

    # create formatter
    basic_formatter = logging.Formatter('%(message)s')
    advanced_formatter = IndentFormatter('%(asctime)s — %(levelname)s — %(name)s.%(funcName)s %(indent)s %(message)s', "%Y-%m-%d %H:%M:%S")

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
