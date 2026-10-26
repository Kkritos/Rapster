"""
Rapster — RAPid cluSTER evolution.

Rapid population synthesis code for compact binary coalescences
and tidal disruption events in dense stellar clusters.

Reference: K. Kritos et al. PRD (2024), arXiv:2210.10055
"""

from importlib.metadata import version as _version, PackageNotFoundError as _PNF
try:
    __version__ = _version("rapster")
except _PNF:
    __version__ = "unknown"

from .cluster_evolution import (
    initialize_cluster,
    compute_cluster_properties,
    compute_timescales,
    form_binaries,
    evolve_interactions,
    compute_external_params,
    evolve_tdes,
    record_evolution,
    flush_records,
    update_cluster,
    print_status,
    write_output,
)
# analysis and plotting helpers are loaded on first use, so that running a
# simulation does not import matplotlib:
_LAZY = {
    'analyze_cluster': ('.analyze_cluster', 'analyze_cluster'),
    'generate_all_plots': ('.plot_cluster', 'generate_all_plots'),
    'load_results': ('.plot_cluster', 'load_results'),
}

def __getattr__(name):
    if name in _LAZY:
        from importlib import import_module
        module, attr = _LAZY[name]
        value = getattr(import_module(module, __name__), attr)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    '__version__',
    # cluster evolution functions:
    'initialize_cluster',
    'compute_cluster_properties',
    'compute_timescales',
    'form_binaries',
    'evolve_interactions',
    'compute_external_params',
    'evolve_tdes',
    'record_evolution',
    'flush_records',
    'update_cluster',
    'print_status',
    'write_output',
    # analysis and plotting:
    'analyze_cluster',
    'generate_all_plots',
    'load_results',
]
