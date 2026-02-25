"""
Functions used to directly access information from nicos.
"""

import socket
import platform
import logging
import subprocess

ON_AMOR = platform.node().startswith('amor')
NICOS_CACHE_DIR = '/home/data/nicosdata/cache/'
GREP = '/usr/bin/grep "value"'


def lookup_nicos_value(key, nicos_key, dtype=float, suffix='', year="2026"):
    # TODO: Implement direct communication to nicos to read the cache
    if nicos_key=='ignore':
        return dtype(0)
    try:
        logging.debug(f'     trying socket request for device {nicos_key}')
        response = nicos_single_request(nicos_key)
        logging.info(f"     using parameter {nicos_key} from nicos cache via socket")
        return dtype(response)
    except Exception as e:
        logging.debug(f'       socket request failed with {e!r}')
    if ON_AMOR:
        logging.debug(f"     trying to extract {nicos_key} from nicos cache files")
        call = f'{GREP} {NICOS_CACHE_DIR}nicos-{nicos_key}/{year}{suffix}'
        try:
            value = str(subprocess.getoutput(call)).split('\t')[-1]
            logging.info(f"     using parameter {nicos_key} from nicos cache file")
            return dtype(value)
        except Exception:
            logging.error(f"       couldn't get value from nicos cache {nicos_key}, {call}")
            return dtype(0)
    else:
        logging.warning(f"     parameter {key} not found, relpace by zero")
        return dtype(0)

def nicos_single_request(device):
    sentinel = b'\n'
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.settimeout(1.0)
        s.connect(('amor', 14869))

        tosend = f'@nicos/{device}/value?\n'

        # write request
        # self.log.debug("get_explicit: sending %r", tosend)
        s.sendall(tosend.encode())

        # read response
        data = b''
        while not data.endswith(sentinel):
            newdata = s.recv(8192)  # blocking read
            if not newdata:
                raise IOError('cache closed connection')
            data += newdata
        s.shutdown(socket.SHUT_RDWR)
    return data.decode('utf-8').split('=')[-1]
