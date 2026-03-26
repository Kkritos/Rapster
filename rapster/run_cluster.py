'''
 Copyright (C) 2026  Konstantinos Kritos <kkritos1@jhu.edu>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''

from .constants import *
from .cluster_evolution import (
    parse_args,
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


def main():
    """Main entry point for running a Rapster cluster simulation."""

    config = parse_args()

    # start initialization clock:
    initialization_time_initial = time.time()
    print('INITIALIZING...')

    state = initialize_cluster(config)

    print('END OF INITIALIZATION. RUNTIME:', "{:.3g}".format(np.abs(time.time() - initialization_time_initial)), 's')
    print('\n')

    print('SIMULATING...')
    print('\n')

    # start simulation clock:
    simulation_time_initial = time.time()

    # Simulation:
    while state['t']<config['t_max'] and state['R_gal']>0 and state['Mcl']>0 and state['mBH'].sum()<fBH_max*state['Mcl']:

        # start local clock:
        local_time_initial = time.time()

        # compute cluster and BH subsystem properties:
        keep_going = compute_cluster_properties(state, config)
        if not keep_going:
            break

        # compute dynamical timescales and adaptive timestep:
        compute_timescales(state, config)

        # binary BH formation (3bb, 2-body capture, exchanges):
        form_binaries(state, config)

        # BBH evolution, pair-pair interactions, triple evolution:
        evolve_interactions(state, config)

        # external parameters (Jacobi radius, escape rate, dynamical friction):
        compute_external_params(state, config)

        # tidal disruption events:
        evolve_tdes(state, config)

        # append evolution record:
        record_evolution(state)

        # cluster structural evolution (mass loss, expansion, time update):
        keep_going = update_cluster(state, config)
        if not keep_going:
            break

        # print runtime information:
        print_status(state, config, local_time_initial, simulation_time_initial)

    write_output(state, config)

    print('END OF SIMULATION. RUNTIME:', "{:.3g}".format(np.abs(time.time() - simulation_time_initial)), 's')
    print('\n')


if __name__ == "__main__":
    main()

# End of file.
