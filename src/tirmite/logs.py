import logging
import sys


def init_logging(loglevel="DEBUG"):
    numeric_level = getattr(logging, loglevel.upper(), None)

    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % loglevel)

    # fmt = "%(asctime)s | %(levelname)8s | %(module)s:%(lineno)s:%(funcName)20s() | %(message)s"
    fmt = "%(asctime)s | %(levelname)s | %(message)s"
    handler_sh = logging.StreamHandler(sys.stderr)
    handler_sh.setFormatter(CustomFormatter(fmt))
    logging.basicConfig(format=fmt, level=numeric_level, handlers=[handler_sh])


class CustomFormatter(logging.Formatter):
    """Logging colored formatter, adapted from https://alexandra-zaharia.github.io/posts/make-your-own-custom-color-formatter-with-python-logging"""

    grey = "\x1b[38;21m"
    blue = "\x1b[38;5;39m"
    yellow = "\x1b[38;5;226m"
    red = "\x1b[38;5;196m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"

    def __init__(self, fmt):
        super().__init__()
        self.fmt = fmt
        self.FORMATS = {
            logging.DEBUG: self.grey + self.fmt + self.reset,
            logging.INFO: self.blue + self.fmt + self.reset,
            logging.WARNING: self.yellow + self.fmt + self.reset,
            logging.ERROR: self.red + self.fmt + self.reset,
            logging.CRITICAL: self.bold_red + self.fmt + self.reset,
        }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
