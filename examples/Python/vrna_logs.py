import RNA
import logging


# Set the default format and log-level for the Python `logging` module
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)


def reroute_logs(e, data):
    """
    The wrapper callback function that receives RNAlib log messages
    and passes them through to the Python 'logging' module

    Parameters
    ----------

    e: dict
        A dictionary that contains the log message, the log level,
        as well as source code file and line number that issued the
        log message
    data: object
        A data object that is passed through to this callback
    """

    log_reroute = {
      RNA.LOG_LEVEL_ERROR: logging.error,
      RNA.LOG_LEVEL_WARNING: logging.warning,
      RNA.LOG_LEVEL_INFO: logging.info,
      RNA.LOG_LEVEL_DEBUG: logging.debug
    }

    # compose an example log message
    message = f"ViennaRNA: \"{e['message']}\" ({e['file_name']}:{e['line_number']})"

    # send log message to appropriate logging channel
    if e['level'] in log_reroute:
        log_reroute[e['level']](message)


# turn-off RNAlib self logging
RNA.log_level_set(RNA.LOG_LEVEL_SILENT)

# attach RNAlib logging to python logging and capture all
# logs with level at least DEBUG
RNA.log_cb_add(reroute_logs, level = RNA.LOG_LEVEL_DEBUG)


# compose an example call that might issue a few debug- or other
# log messages
md = RNA.md(circ=1, gquad=1)
s = RNA.random_string(10, "ACGU")

fc = RNA.fold_compound(s, md)

fc.mfe()
