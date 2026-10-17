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
    update_cluster,
    print_status,
    write_output,
)
from .analyze_cluster import analyze_cluster
from .plot_cluster import generate_all_plots, load_results

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
    'update_cluster',
    'print_status',
    'write_output',
    # analysis and plotting:
    'analyze_cluster',
    'generate_all_plots',
    'load_results',
]
