"""
Setup for the logging of eos.
"""
import sys
import logging
import logging.handlers

def setup_logging():
    logger = logging.getLogger()  # logging.getLogger('quicknxs')
    logger.setLevel(logging.DEBUG)
    # rename levels to make clear warning is can be a normal message
    logging.addLevelName(logging.INFO, 'VERB')
    logging.addLevelName(logging.WARNING, 'MESG')
    # setting up a logger for console output
    console = logging.StreamHandler(sys.__stdout__)
    console.name = 'console'
    formatter = logging.Formatter('# %(message)s')
    console.setFormatter(formatter)
    console.setLevel(logging.WARNING)
    logger.addHandler(console)

    # if os.path.exists('amor_eos.log'):
    #     rollover = True
    # else:
    #     rollover = False
    logfile = logging.handlers.RotatingFileHandler('amor_eos.log', encoding='utf8', mode='w',
                                                   maxBytes=200*1024**2, backupCount=20)
    # if rollover: logfile.doRollover()
    formatter = logging.Formatter(
        '[%(levelname).4s] - %(asctime)s - %(filename)s:%(lineno)i:%(funcName)s %(message)s',
            '')
    logfile.setFormatter(formatter)
    logfile.setLevel(logging.DEBUG)
    logger.addHandler(logfile)

def update_loglevel(verbose=0):
    if verbose==1:
        logging.getLogger().handlers[0].setLevel(logging.INFO)
    if verbose>1:
        console = logging.getLogger().handlers[0]
        console.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(levelname).1s %(message)s')
        console.setFormatter(formatter)
