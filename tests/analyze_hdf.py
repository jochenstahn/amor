"""
Small helper to find information about hdf datafiles for debugging
"""

import h5py

def rec_tree(group, min_size=128):
    if hasattr(group, 'keys'):
        output = {}
        total_size = 0
        for key in group.keys():
            subkeys, size = rec_tree(group[key], min_size)
            total_size += size
            if size>min_size:
                if subkeys:
                    output[key] = subkeys
                else:
                    output[key] = size
        return output, size
    elif hasattr(group, 'size'):
        return None, group.size
    else:
        return None, 0

if __name__ == "__main__":
    import sys
    for fi in sys.argv[1:]:
        print(fi)
        print(rec_tree(sys.argv[1]))

    