"""
Functions used to directly access information from nicos.
"""

import platform
import logging
import subprocess

if  platform.node().startswith('amor'):
    NICOS_CACHE_DIR = '/home/data/nicosdata/cache/'
    GREP = '/usr/bin/grep "value"'
else:
    NICOS_CACHE_DIR = None


def lookup_nicos_value(key, nicos_key, dtype=float, suffix='', year="2026"):
    # TODO: Implement direct communication to nicos to read the cache
    if nicos_key=='ignore':
        return dtype(0)
    if NICOS_CACHE_DIR:
        try:
            logging.info(f"     using parameter {nicos_key} from nicos cache")
            call = f'{GREP} {NICOS_CACHE_DIR}nicos-{nicos_key}/{year}{suffix}'
            value = str(subprocess.getoutput(call)).split('\t')[-1]
            return dtype(value)
        except Exception:
            logging.error(f"Couldn't get value from nicos cache {nicos_key}, {call}")
            return dtype(0)
    else:
        logging.warning(f"     parameter {key} not found, relpace by zero")
        return dtype(0)
