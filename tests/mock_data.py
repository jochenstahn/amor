"""
Generate a mock dataset in memory for running unit tests.
"""

import h5py
import numpy as np

MOCK_METADATA = {
    'title': 'Testdata',
    'proposal_id': 'none',
    'user/name': 'test user',
    'user/email': 'test@user.de',
    'sample/name': 'test sample',
    'sample/model': 'air | Fe 12 | Si',
    'Amor/source/name': 'SINQ',
    'start_time': '2025-01-01 00:00:01',
    }
MOCK_META_TYPED = {
    'Amor/chopper/pair_separation': (1000.0, np.float32),
    'Amor/detector/transformation/distance': (4000.0, np.float64),
    'Amor/instrument_control_parameters/kappa': (1000.0, np.float64),
    'Amor/instrument_control_parameters/kappa_offset': (1000.0, np.float64),
    'Amor/instrument_control_parameters/div': (1.6, np.float64),
    'Amor/chopper/ch1_trigger_phase': (-9.1, np.float64),
    'Amor/chopper/ch2_trigger_phase': (6.75, np.float64),
    'Amor/chopper/ch2_trigger/event_time_zero': ([0.0]*10, np.uint64),
    'Amor/chopper/ch2_trigger/event_time_offset': ([0.0]*10, np.uint32),
    'Amor/chopper/rotation_speed': (500.0, np.float64),
    'Amor/chopper/phase': (0.0, np.float64),
    'Amor/polarization/configuration/value': (0.0, np.float64),
    }

def mock_data(mu=1.0, nu=2.0):
    hdf = h5py.File.in_memory() # requires h5py >=3.13
    ds = hdf.create_group('entry1')
    for key, value in MOCK_METADATA.items():
        ds.create_dataset(key, data=np.array([value.encode('utf-8')]))
    for key, (value, dtype) in MOCK_META_TYPED.items():
        if type(value) is list:
            ds.create_dataset(key, data=np.array(value), dtype=dtype)
        else:
            ds.create_dataset(key, data=np.array([value]), dtype=dtype)

    ds.create_dataset('Amor/instrument_control_parameters/mu', np.array([mu]), dtype=np.float64)
    ds.create_dataset('Amor/instrument_control_parameters/nu', np.array([nu]), dtype=np.float64)

    return hdf

def compare_with_real_data(fname):
    hdf = h5py.File(fname, 'r')
    ds = hdf['entry1']
    for key, value in MOCK_METADATA.items():
        try:
            ds[key][0].decode('utf-8')
        except KeyError:
            print(f'/entry1/{key} does not exist in file')
    for key, (value, dtype) in MOCK_META_TYPED.items():
        try:
            item = ds[key]
        except KeyError:
            print(f'/entry1/{key} does not exist in file')
        else:
            if item.dtype != dtype:
                print(f'/entry1/{key} does not match {dtype}, dataset is {item.dtype}')
