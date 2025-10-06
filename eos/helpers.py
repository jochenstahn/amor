"""
Helper functions used during calculations. Uses numba enhanced functions if available, otherwise numpy based
fallback is imported.
"""

try:
    from .helpers_numba import merge_frames, extract_walltime, filter_project_x, calculate_derived_properties_focussing
except ImportError:
    from .helpers_fallback import merge_frames, extract_walltime, filter_project_x, calculate_derived_properties_focussing

