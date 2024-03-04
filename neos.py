#!/usr/bin/env python
"""
eos reduces measurements performed on Amor@SINQ, PSI

Author: Jochen Stahn

conventions (not strictly followed, yet):
- array names end with the suffix '_x[y]' with the meaning
    _e  = events
    _tof
    _l  = lambda
    _t  = theta
    _z  = detector z
    _lz = (lambda, detector z)
    _q  = q_z
- to come
"""

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

    # read command line arguments and generate classes holding configuration parameters
    reader_config, reduction_config, output_config = command_line_options()
    # Create reducer with these arguments
    reducer = AmorReduction(reader_config, reduction_config, output_config)
    # Perform actual reduction
    reducer.reduce()

if __name__ == '__main__':
    main()
