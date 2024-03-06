#!/usr/bin/env python
"""
eos reduces measurements performed on Amor@SINQ, PSI

Author: Jochen Stahn (algorithms, python draft),
        Artur Glavic (structuring and optimisation of code)

conventions (not strictly followed, yet):
- array names end with the suffix '_x[y]' with the meaning
    _e  = events
    _tof
    _l  = lambda
    _t  = theta
    _z  = detector z
    _lz = (lambda, detector z)
    _q  = q_z
"""

import logging

from libeos.command_line import command_line_options
from libeos.logconfig import setup_logging
from libeos.reduction import AmorReduction

#=====================================================================================================
# TODO:
# - calculate resolution using the chopperPhase
# - deal with background correction
# - format of 'call' + add '-Y' if not supplied
#=====================================================================================================

def main():
    setup_logging()
    logging.warning('######## eos - data reduction for Amor ########')

    # read command line arguments and generate classes holding configuration parameters
    config = command_line_options()
    # Create reducer with these arguments
    reducer = AmorReduction(config)
    # Perform actual reduction
    reducer.reduce()

    logging.info('######## eos - finished ########')

if __name__ == '__main__':
    main()
