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
from .functions import *
from .stellar_evolution import *
from .three_body_binary import three_body_binary
from .binary_evolution import evolve_BBHs
from .triples import evolve_triples
from .two_body_capture import two_body_capture
from .exchanges import StarStar_to_BHstar, BHstar_to_BBH
from .tidal_disruptions import BH_TidalDisruptions, single_BH_TidalDisruption
from .bhpool import BHPool, FAM
from .compact_accretion import *
from .remnant import configure_kick_model
from .ns_formation import NSFormationTable, sample_ns_formation_times
from .fastrandom import poisson_draw, maxwell_draws
from .snapshots import SnapshotStream, dump_streaming
from .arraycache import ElementwiseCache, RowLog, StreamedRowLog, write_rows
from .populations import StellarPopulations, POPULATION_KEYS
import math
from zipfile import BadZipFile


def _load_initial_bh_masses(input_file):
    """Load and validate initial BH masses from a NumPy .npz archive."""
    try:
        archive = np.load(input_file, allow_pickle=False)
    except (BadZipFile, EOFError, OSError, ValueError) as error:
        raise ValueError(f'could not read initial black-hole input file {input_file!r}') from error

    if not isinstance(archive, np.lib.npyio.NpzFile):
        raise ValueError(f'initial black-hole input file {input_file!r} must be a .npz archive')

    try:
        if 'mBH_ini' not in archive.files:
            raise ValueError(f'initial black-hole input file {input_file!r} must contain an mBH_ini array')
        masses = np.asarray(archive['mBH_ini'])
    finally:
        archive.close()

    if masses.ndim != 1:
        raise ValueError('mBH_ini must be a one-dimensional array of initial BH masses in Msun')
    if not np.issubdtype(masses.dtype, np.number) or np.iscomplexobj(masses):
        raise ValueError('mBH_ini must contain real numeric masses in Msun')

    masses = masses.astype(float, copy=False)
    if not np.all(np.isfinite(masses)) or np.any(masses <= 0):
        raise ValueError('mBH_ini must contain only finite, positive masses in Msun')

    return masses


def _sample_bh_masses_and_kicks(Mcl, Z, config):
    """BHs formed by a stellar population of initial mass Mcl [Msun] and metallicity Z:
    (masses [Msun], natal kicks [km/s]) of all BHs, before retention, drawn as for the
    initial population (-BMD, -RP, -K, Kroupa IMF between -mm and -mM). Used for the
    initial population (initialize_cluster) and for later populations (populations.py);
    the random numbers are drawn in the same order as before this was a function."""
    m_min = config['m_min']; m_max = config['m_max']
    BMD = config['BMD']; RP = config['RP']; NKP = config['NKP']; wSN_kick = config['wSN_kick']
    min_1g_bh_mass = config['min_1g_bh_mass']; max_1g_bh_mass = config['max_1g_bh_mass']

    # number of massive stars:
    N_massive = int(Mcl * integrate.quad(IMF_kroupa, mM_min, m_max)[0] \
                    / integrate.quad(lambda x: x * IMF_kroupa(x), m_min, m_max)[0])

    if BMD==0:
        # default: Kroupa IMF + stellar collapse + SN kicks

        # massive star masses drawn from Kroupa (2002) IMF:
        u_massive = np.random.rand(N_massive)
        m_massive = (u_massive * (m_max**(alphaIMF + 1) - mM_min**(alphaIMF + 1)) \
                     + mM_min**(alphaIMF + 1))**(1 / (alphaIMF + 1))

        # remnant masses:
        if   RP==0:
            m_rem = Mrem_SEVNdelayed(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        elif RP==1:
            m_rem = Mrem_F12d(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        elif RP==2:
            m_rem = Mrem_SEVNrapid(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        elif RP==3:
            m_rem = Mrem_F12r(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)

        # Separate BHs from other remnants:
        mBH = m_rem[m_rem > mBH_min]

        # compute supernova kicks:
        if NKP==0:
            if   RP==0:
                fallback = np.vectorize(f_fb_delayed, otypes=[float])(m_massive[m_rem > mBH_min], np.vectorize(M_CO_SEVN, otypes=[float])(m_massive[m_rem > mBH_min], Z))
                vSN_kick = (1 - fallback) * np.vectorize(get_SN_kick, otypes=[float])(1.4 * np.ones(mBH.size), wSN_kick)
            elif RP==1:
                fallback = np.vectorize(f_fb_delayed, otypes=[float])(m_massive[m_rem > mBH_min], np.vectorize(M_CO_SSE, otypes=[float])(m_massive[m_rem > mBH_min], Z))
                vSN_kick = (1 - fallback) * np.vectorize(get_SN_kick, otypes=[float])(1.4 * np.ones(mBH.size), wSN_kick)
            elif RP==2:
                fallback = np.vectorize(f_fb_rapid, otypes=[float])(m_massive[m_rem > mBH_min], np.vectorize(M_CO_SEVN, otypes=[float])(m_massive[m_rem > mBH_min], Z))
                vSN_kick = (1 - fallback) * np.vectorize(get_SN_kick, otypes=[float])(1.4 * np.ones(mBH.size), wSN_kick)
            elif RP==3:
                fallback = np.vectorize(f_fb_rapid, otypes=[float])(m_massive[m_rem > mBH_min], np.vectorize(M_CO_SSE, otypes=[float])(m_massive[m_rem > mBH_min], Z))
                vSN_kick = (1 - fallback) * np.vectorize(get_SN_kick, otypes=[float])(1.4 * np.ones(mBH.size), wSN_kick)
        elif NKP==1:
            vSN_kick = np.vectorize(get_SN_kick, otypes=[float])(mBH, wSN_kick)


    else:
        # For BMD=1,2: run Kroupa path to determine N_BH, then replace masses.
        u_massive = np.random.rand(N_massive)
        m_massive = (u_massive * (m_max**(alphaIMF + 1) - mM_min**(alphaIMF + 1)) \
                     + mM_min**(alphaIMF + 1))**(1 / (alphaIMF + 1))
        if   RP==0:
            m_rem = Mrem_SEVNdelayed(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        elif RP==1:
            m_rem = Mrem_F12d(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        elif RP==2:
            m_rem = Mrem_SEVNrapid(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        elif RP==3:
            m_rem = Mrem_F12r(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        N_BH_init = (m_rem > mBH_min).sum()

        if BMD==1:
            # uniform BH mass distribution in [min_1g_bh_mass, max_1g_bh_mass]:
            mBH = np.random.uniform(min_1g_bh_mass, max_1g_bh_mass, N_BH_init)

        elif BMD==2:
            # Salpeter power law (m^-2.35) in [min_1g_bh_mass, max_1g_bh_mass]:
            alpha_salpeter = -2.35
            u = np.random.rand(N_BH_init)
            mBH = (u * (max_1g_bh_mass**(alpha_salpeter + 1) - min_1g_bh_mass**(alpha_salpeter + 1)) \
                   + min_1g_bh_mass**(alpha_salpeter + 1))**(1 / (alpha_salpeter + 1))

        elif BMD==3:
            # log-uniform BH mass distribution in [min_1g_bh_mass, max_1g_bh_mass]:
            mBH = 10**(np.random.uniform(np.log10(min_1g_bh_mass), np.log10(max_1g_bh_mass), N_BH_init))

        # Apply momentum-conservation SN kicks (retention is decided by the caller).
        # Fallback kicks (NKP==0) are not available for BMD>0 because there is no
        # stellar progenitor or CO core mass to compute the fallback fraction from.
        # Momentum-conservation kicks only depend on the BH mass, so they apply regardless.
        vSN_kick = np.vectorize(get_SN_kick, otypes=[float])(mBH, wSN_kick)

    return mBH, vSN_kick


def _draw_ns_masses(n, with_NSs, EoS):
    """Draw n NS masses [Msun] from the distribution selected by -NS.

    with_NSs = 1: monochromatic (neutron_star_mass)
    with_NSs = 2: bimodal, Rocha et al. (2023), truncated to the EoS range
    with_NSs = 3: uniform over the EoS range [M_min, M_TOV]
    """
    if with_NSs==1:
        # monochromatic NS mass distribution:
        return neutron_star_mass * np.ones(n)

    # define NS range (min/max masses):
    eos = get_eos(EoS)
    M_NS_min = eos.M_min
    M_NS_max = eos.M_TOV

    if with_NSs==2:
        # bimodal mass distribution from Rocha et al. 2023:
        mNS = []
        while len(mNS) < n:
            # generate a single NS mass:
            m_NS = np.random.normal(loc=1.351, scale=0.084) if np.random.rand() < 0.539 else np.random.normal(loc=1.816, scale=0.260)
            if M_NS_min < m_NS < M_NS_max:
                mNS.append(m_NS)
        return np.array(mNS)

    # uniform mass distribution over the EoS-allowed NS range:
    return np.random.uniform(M_NS_min, M_NS_max, n)


def compact_object_mass(state):
    """Total mass [Msun] of the compact objects (BHs and retained NSs) in the cluster:
    single ones (the pool) and those in binaries, hierarchical triples and BH-star
    pairs (the placeholder first rows of these arrays are zero)."""
    return (float(state['BH'].m.sum())
            + float(state['binaries'][:, 4].sum() + state['binaries'][:, 5].sum())
            + float(state['triples'][:, 4:7].sum())
            + float(state['pairs'][:, 1].sum()))


# rows of the evolution, hardening and TDE records written to their files at a time
# (memory for at most this many rows per record; ~1-16 MB):
STREAM_CHUNK_ROWS = 131072


def flush_records(state):
    """Write the evolution, hardening, TDE and merger rows that can no longer change to
    their output files, once a chunk of them has accumulated (called once per time step)."""
    for key in ('tdes', 'hardening', 'evolution', 'mergers', 'populations'):
        record = state.get(key)
        if isinstance(record, StreamedRowLog):
            record.flush()


def _new_records(config, with_populations):
    """The records of a run, each with one placeholder row: mergers, evolution, hardening,
    TDEs and (with_populations) populations. Their rows are written to their files during
    the run, in chunks (see flush_records and StreamedRowLog), with the same text as
    writing them all at the end."""
    results_dir = os.path.join(os.getcwd(), config.get('results_folder_name', 'Results'))

    def streamed(ncols, name, keys, indicator, keep, prepare=None):
        return StreamedRowLog(ncols, [0.0] * ncols,
                              os.path.join(results_dir, config.get(name + '_file', name) + '.txt'),
                              '# ' + ' '.join(keys) + '\n', _write_rows, keep=keep,
                              chunk=STREAM_CHUNK_ROWS, prepare=prepare,
                              enabled=config.get(indicator, 1) == 1)
    return {
        # (each merger line ends with the final cluster properties, added when the run ends)
        'mergers': streamed(27, 'mergers', merger_keys, 'Mi', keep=0),
        # (the last evolution rows are kept: write_output drops the last one or two rows,
        # depending on how the run ends; the last hardening row can still change)
        'evolution': streamed(71, 'evolution', evolution_keys, 'Ei', keep=2),
        'hardening': streamed(12, 'hardening', hardening_keys, 'Hi', keep=1),
        'tdes': streamed(18, 'tdes', tdes_keys, 'Ti', keep=0, prepare=fill_R_mb),
        'populations': (streamed(len(POPULATION_KEYS), 'populations', POPULATION_KEYS, 'POPi', keep=0)
                        if with_populations else None),
    }


RECORD_KEYS = ('mergers', 'evolution', 'hardening', 'tdes', 'populations')


def initialize_cluster(config):
    """Set up the initial cluster state from the configuration.

    Samples stellar masses from the Kroupa IMF, computes remnant masses using
    the selected prescription, applies supernova natal kicks, retains BHs below
    the escape velocity, initializes BH spins (uniform, monochromatic, or beta),
    optionally includes neutron stars and a massive BH seed, and allocates all
    data arrays (binaries, mergers, evolution, hardening, triples, pairs, tdes).

    Args:
        config (dict): Configuration dictionary from parse_args().

    Returns:
        dict: Simulation state dictionary containing all mutable variables
            (cluster properties, BH arrays, data arrays, counters, time variables).
    """
    # unpack configuration parameters into local variables:
    N = config['N']
    rh = config['rh']
    m_min = config['m_min']
    m_max = config['m_max']
    Z = config['Z']
    seed = config['seed']
    wSN_kick = config['wSN_kick']
    R_gal = config['R_gal']
    v_gal = config['v_gal']
    zCl_form = config['zCl_form']
    n_star = config['n_star']
    NKP = config['NKP']
    fb = config['fb']
    Bi = config['Bi']
    input_BH_file = config['input_BH_file']
    s1g_max = config['s1g_max']
    SD = config['SD']
    RP = config['RP']
    with_NSs = config['with_NSs']
    NSt = config.get('NSt', 0)
    M_BH0 = config['M_BH0']
    s_BH0 = config['s_BH0']
    dt_min = config['dt_min']
    BMD = config['BMD']
    min_1g_bh_mass = config['min_1g_bh_mass']
    max_1g_bh_mass = config['max_1g_bh_mass']
    EoS = config['EoS']
    SFE = config['SFE']
    t_ge = config['t_ge']   # gas expulsion timescale [Myr]
    c_s = config['c_s']
    f_Edd = config['f_Edd']

    # initialize pseudo-random number generator:
    np.random.seed(seed)

    # initial half-mass radius:
    rh0 = rh

    # initial stellar number density:
    n_star0 = n_star

    # average stellar mass:
    Kroupa_norm = 1 / integrate.quad(IMF_kroupa, m_min, m_max)[0]
    m_avg = Kroupa_norm * integrate.quad(lambda x: x * IMF_kroupa(x), m_min, m_max)[0]

    # initial average stellar mass:
    m_avg0 = m_avg

    # initial cluster mass in stars:
    Mcl = N * m_avg
    Mcl0 = Mcl # this includes only the mass in stars (not gas)

    # initial residual gas mass from star formation efficiency:
    M_gas0 = (1 - SFE) / SFE * Mcl0
    M_gas = M_gas0


    # gas remaining at BH formation; reused for NSs as an upper bound
    M_gas_BHform = M_gas0 * np.exp(- tBH_form / t_ge)

    # mass that retains the BHs (and, for -NSt 0, the NSs) of the first population, which
    # form at tBH_form: its stars and remaining gas, and those of the stellar populations
    # born before then (-POP), if any:
    M_retain = Mcl + M_gas_BHform
    M_early = 0.0
    for p in config.get('populations', ()):
        if p['t_birth'] < tBH_form:
            M_k = p['N'] * m_avg
            M_early += M_k + (1 - p['SFE']) / p['SFE'] * M_k * np.exp(- (tBH_form - p['t_birth']) / t_ge)
    if M_early > 0:
        M_retain = M_retain + M_early

    # initial galactocentric radius:
    R_gal0 = R_gal

    # initial multimass relaxation factor:
    psi0 = integrate.quad(lambda x: x**(5/2) * IMF_kroupa(x), m_min, m_max)[0] \
        / integrate.quad(IMF_kroupa, m_min, m_max)[0] / m_avg**(5/2)

    # initial relaxation timescale:
    t_rlx0 = t_relax(Mcl0 + M_gas0, rh0, m_avg, psi0, np.log(np.e + lc * N))

    # core collapse timescale:
    t_cc = k_cc * t_rlx0

    # number of binary stars:
    Nb = int(N * fb / 2)

    # minimum binary-star separation:
    ab_min = 3 * R_sun

    # maximum binary-star separation:
    ab_max = (3 / np.pi / 4 / n_star)**(1/3) / 10

    # initialize binary star population and filter hard binaries:
    if Nb>0:
        # semimajor axes of binaries (loguniform in [ab_min, ab_max]):
        ab = 10**(np.random.uniform(np.log10(ab_min), np.log10(ab_max), Nb))

        # stellar velocity dispersion:
        v_star = np.sqrt(0.4 * G_Newton * Mcl / rh)

        # binary star hard-soft boundary:
        ab_hard = G_Newton * m_avg / v_star

        # filter only hard binary stars:
        ab = ab[ab < ab_hard]
    else:
        ab = np.array([])
        v_star = np.sqrt(0.4 * G_Newton * Mcl / rh)

    # BH masses and natal kicks of the population (see _sample_bh_masses_and_kicks), and
    # the BHs retained (kick below the escape velocity, with the gas left at BH formation):
    mBH, vSN_kick = _sample_bh_masses_and_kicks(Mcl, Z, config)
    mBH = mBH[vSN_kick < v_esc(M_retain, rh)]

    # optionally override BH masses from external file:
    if Bi==1:
        mBH = _load_initial_bh_masses(input_BH_file)
        mBH = mBH + 0.01 * np.random.rand(mBH.size)

    # core collapse of the BH subsystem (onset of balanced evolution): BHs mass-segregate on
    # t_seg ~ (<m>/<m_BH>) t_rh (Spitzer 1969), with t_rh the initial single-mass (psi = 1)
    # relaxation time and k_cc an order-unity constant. If no BHs are retained, keep the
    # IMF-based estimate computed above.
    if mBH.size > 0:
        t_rlx0_single = t_relax(Mcl0 + M_gas0, rh0, m_avg, 1.0, np.log(np.e + lc * N))
        t_cc = k_cc * m_avg / np.mean(mBH) * t_rlx0_single

    # BH spins:
    if SD==0:
        sBH = np.random.uniform(0, s1g_max, mBH.size)
    elif SD==1:
        sBH = s1g_max * np.ones(mBH.size)
    elif SD==2:
        # drawn from the seeded global stream (np.random.seed(seed) above), so that
        # runs with a given -S are reproducible:
        sBH = np.random.beta(1.4, 3.6, mBH.size) * s1g_max

    # NSs that have not formed yet (used only when NSt > 0); filled below:
    mNS_pending = np.array([])
    tNS_pending = np.array([])
    vNS_pending = np.array([])

    # optionally form and retain neutron stars based on natal kicks:
    if with_NSs>0:
        f_NS_form = Kroupa_norm*integrate.quad(IMF_kroupa, 8, 18)[0]
        N_NS_form = int(f_NS_form*N)
        v_NS_natal = maxwell_draws(np.sqrt(3)*wSN_kick, N_NS_form)   # = maxwell.rvs(loc=0, scale=..., size=N_NS_form)

        if NSt==0:
            # legacy: all NSs form at tBH_form; tBH_form is a lower bound on the
            # NS formation time (=> gas upper bound). Random draws are kept in the
            # same order as in earlier versions, so results are bit-identical.
            N_NS_ret = v_NS_natal[v_NS_natal < v_esc(M_retain, rh)].size

            if N_NS_ret > 0:
                mNS = _draw_ns_masses(N_NS_ret, with_NSs, EoS)
                mBH = np.concatenate((mBH, mNS))
                # initialize NS spins to zero:
                sBH = np.concatenate((sBH, np.zeros(N_NS_ret)))

        else:
            # extended NS formation: draw masses and natal kicks for every NS that
            # forms, and formation times from the COMPAS table (marginal p(t) for
            # NSt=1, p(t | m_NS) for NSt=2). Retention is decided at formation, with
            # the escape velocity at that time (form_neutron_stars, main loop).
            if N_NS_form > 0:
                mNS_pending = _draw_ns_masses(N_NS_form, with_NSs, EoS)
                tNS_pending = sample_ns_formation_times(mNS_pending, NSt, tBH_form, table=NSFormationTable())
                vNS_pending = v_NS_natal

                # sort by formation time, so the main loop takes NSs from the front:
                order = np.argsort(tNS_pending, kind='stable')
                mNS_pending = mNS_pending[order]
                tNS_pending = tNS_pending[order]
                vNS_pending = vNS_pending[order]

    # optionally add a massive BH seed to the initial conditions:
    if M_BH0>0:
        mBH = np.concatenate((mBH, [M_BH0]))
        sBH = np.concatenate((sBH, [s_BH0]))

    # BH generations:
    gBH = np.ones(mBH.size)
    hBH = np.zeros(mBH.size)

    # data arrays:
    binaries = np.zeros(shape=(1, 15))
    pairs = np.zeros(shape=(1, 5))
    triples = np.zeros(shape=(1, 24))
    configure_kick_model(config['recoil_kick_model'])

    # per-step BH snapshots (-BOi 1) are streamed to temporary files next to the results:
    results_dir = os.path.join(os.getcwd(), config.get('results_folder_name', 'Results'))

    # records (merger, evolution, hardening and TDE rows), written during the run:
    records = _new_records(config, with_populations=bool(config.get('populations')))
    mergers = records['mergers']; evolution = records['evolution']
    hardening = records['hardening']; tdes = records['tdes']
    snapshot_list = (lambda: SnapshotStream(results_dir)) if config.get('BOi', 0) == 1 else list

    # bundle all mutable simulation variables into the state dictionary:
    state = {
        # cluster properties:
        'Mcl': Mcl, 'Mcl0': Mcl0, 'rh': rh, 'rh0': rh0, 'R_gal': R_gal, 'R_gal0': R_gal0,
        'v_gal': v_gal, 'n_star': n_star, 'n_star0': n_star0,
        'm_avg': m_avg, 'm_avg0': m_avg0, 'v_star': v_star,
        'M_gas': M_gas, 'M_gas0': M_gas0, 't_ge': t_ge,
        'Kroupa_norm': Kroupa_norm, 'Nb': Nb, 'ab': ab,
        't_cc': t_cc, 't_rlx': 0.0,
        'N': N, 'Z': Z,
        # single BHs (and retained NSs): masses, spins, generations, TDE counts:
        'BH': BHPool(mBH, sBH, gBH, hBH),
        'mBH_avg': 0.0, 'vBH': 0.0,
        # NSs not yet formed (NSt > 0), sorted by formation time [Msun, Myr, km/s]:
        'mNS_pending': mNS_pending, 'tNS_pending': tNS_pending, 'vNS_pending': vNS_pending,
        'N_NS_formed': 0, 'N_NS_retained': 0,
        # data arrays:
        'binaries': binaries, 'pairs': pairs, 'triples': triples,
        'mergers': mergers, 'evolution': evolution, 'hardening': hardening, 'tdes': tdes,
        # counters:
        'N_BH': 0, 'N_BBH': 0, 'N_3bb': 0, 'N_me': 0, 'N_meRe': 0, 'N_meEj': 0,
        'N_2cap': 0, 'N_3cap': 0, 'N_dis': 0, 'N_ex': 0, 'N_BHej': 0, 'N_BBHej': 0,
        'N_iter': 0, 'N_bb': 0, 'N_meFi': 0, 'N_me2b': 0,
        'N_ex1': 0, 'N_ex2': 0, 'N_BHstar': 0, 'N_pp': 0,
        'N_Triples': 0, 'N_ZLK': 0, 'N_WD': 0,
        'N_tdeBHWD': 0, 'N_tdeBHstar': 0, 'N_tdeBBHstar': 0, 'N_hardening': 0,
        # time:
        't': 0, 'z': zCl_form, 'dt': dt_min, 'zCl_form': zCl_form, 'seed': seed,
        # aux:
        'i_aux1': 0,
        # tracking lists (snapshots for -BOi 1 are kept on disk, see snapshots.py):
        'simulation_times': [], 'black_hole_masses': snapshot_list(), 'black_hole_spins': snapshot_list(),
        'black_hole_generations': snapshot_list(), 'black_hole_tdes': snapshot_list()
    }

    # stellar populations (populations.py): the initial one, and those born later (-POP):
    state['pops'] = StellarPopulations(Mcl, Z, SFE, config.get('populations', ()))
    set_field_star_mixture(None, None)
    # BHs and NSs of later populations that have not formed yet, sorted by formation time
    # (masses [Msun], spins, formation times [Myr], natal kicks [km/s], 1 for NSs and 0
    # for BHs):
    state['mCO_pending'] = np.array([]); state['sCO_pending'] = np.array([])
    state['tCO_pending'] = np.array([]); state['vCO_pending'] = np.array([])
    state['nCO_pending'] = np.array([])
    state['N_pop_CO_formed'] = 0; state['N_pop_CO_retained'] = 0
    # one row per population and time step, written when populations are configured (-POP):
    state['populations'] = records['populations']
    # (number of evolution rows, number of population rows incl. the placeholder row) after
    # the last few evolution rows:
    state['population_marks'] = []

    return state


def restore_cluster(checkpoint, config, with_populations=False):
    """The simulation state saved in a checkpoint (restart.py), for a restarted run with
    configuration config (results in config['results_folder_name']): the output records
    continue from the checkpoint (StreamedRowLog.restore), and the random number generator
    from its saved state. with_populations: give the run a populations record even if the
    checkpoint has none (a population is added at the restart)."""
    configure_kick_model(config['recoil_kick_model'])
    set_field_star_mixture(None, None)
    state = dict(checkpoint['state'])
    saved = checkpoint['records']
    records = _new_records(config, with_populations=with_populations or saved.get('populations') is not None)
    for key, rec in saved.items():
        if rec is not None and records.get(key) is not None:
            records[key].restore(*rec)
    state.update(records)
    results_dir = os.path.join(os.getcwd(), config.get('results_folder_name', 'Results'))
    snapshot_list = (lambda: SnapshotStream(results_dir)) if config.get('BOi', 0) == 1 else list
    state['simulation_times'] = []
    for key in ('black_hole_masses', 'black_hole_spins', 'black_hole_generations', 'black_hole_tdes'):
        state[key] = snapshot_list()
    state['n_snapshots_before'] = checkpoint['n_snapshots']
    np.random.set_state(checkpoint['rng'])
    return state


def add_black_holes(state, config, masses, spins=None):
    """Add BHs (masses [Msun]; spins, or drawn as for first-generation BHs with -s and -SD)
    to the single compact objects of a restarted run: generation 1, no TDEs."""
    masses = np.asarray(masses, dtype=float)
    n = masses.size
    if spins is None:
        SD = config['SD']; s1g_max = config['s1g_max']
        if SD==0:
            spins = np.random.uniform(0, s1g_max, n)
        elif SD==1:
            spins = s1g_max * np.ones(n)
        else:
            spins = np.random.beta(1.4, 3.6, n) * s1g_max
    spins = np.asarray(spins, dtype=float)
    if n > 0:
        state['BH'].add_many(masses, spins, np.ones(n), np.zeros(n))
        if state['i_aux1'] == 1:
            state['N_BH'] = state['N_BH'] + n


def form_populations(state, config):
    """Give birth to the stellar populations scheduled for the current time (-POP, or a
    population added at a restart). A new population of N stars (Kroupa IMF, as the first
    one) adds its stellar mass to the cluster, its residual gas (star formation efficiency
    SFE) to the gas of the cluster, its hard stellar binaries to those available for
    exchanges, and its BHs and NSs to the compact objects still to form (BHs at an age of
    tBH_form, NSs as set by -NSt; see form_pending_compact_objects). Updates state in-place."""
    pops = state['pops']
    if not pops.scheduled:
        return
    t = state['t']
    for p in pops.due(t):
        if not pops.multi:
            # from now on the initial population's mass is kept by the populations:
            pops.M[0] = state['Mcl']
        m_avg0 = state['m_avg0']
        N_k = p['N']
        M_k = N_k * m_avg0
        Z_k = p['Z']
        SFE_k = p['SFE']
        pops.add(t, N_k, M_k, Z_k, SFE_k, p['x'], p['n_mix'])
        state['Mcl'] = state['Mcl'] + M_k
        state['m_avg'] = state['Mcl'] / sum(pops.numbers(t, m_avg0))

        # residual gas of the new population, expelled with the cluster's gas (-tge):
        state['M_gas'] = state['M_gas'] + (1 - SFE_k) / SFE_k * M_k

        # hard stellar binaries of the new population (as at initialization):
        Nb_k = int(N_k * config['fb'] / 2)
        if Nb_k > 0 and state['n_star'] > 0:
            ab_min = 3 * R_sun
            ab_max = (3 / np.pi / 4 / state['n_star'])**(1/3) / 10
            ab_k = 10**(np.random.uniform(np.log10(ab_min), np.log10(ab_max), Nb_k))
            v_star = np.sqrt(0.4 * G_Newton * state['Mcl'] / state['rh'])
            ab_k = ab_k[ab_k < G_Newton * state['m_avg'] / v_star]
            state['ab'] = np.concatenate((state['ab'], ab_k))
            state['Nb'] = state['Nb'] + Nb_k

        # BHs: masses and natal kicks drawn as for the first population; they form at an
        # age of tBH_form, and are retained if the kick is below the escape velocity then:
        mBH, vBH_kick = _sample_bh_masses_and_kicks(M_k, Z_k, config)
        SD = config['SD']; s1g_max = config['s1g_max']
        if SD==0:
            sBH = np.random.uniform(0, s1g_max, mBH.size)
        elif SD==1:
            sBH = s1g_max * np.ones(mBH.size)
        else:
            sBH = np.random.beta(1.4, 3.6, mBH.size) * s1g_max
        m_new = [mBH]; s_new = [sBH]; t_new = [t + tBH_form * np.ones(mBH.size)]; v_new = [vBH_kick]
        n_new = [np.zeros(mBH.size)]

        # NSs (zero spin), forming at an age of tBH_form (-NSt 0) or at the ages of the
        # COMPAS distribution (-NSt 1, 2):
        with_NSs = config['with_NSs']
        if with_NSs > 0:
            f_NS_form = state['Kroupa_norm'] * integrate.quad(IMF_kroupa, 8, 18)[0]
            N_NS = int(f_NS_form * N_k)
            if N_NS > 0:
                vNS = maxwell_draws(np.sqrt(3) * config['wSN_kick'], N_NS)
                mNS = _draw_ns_masses(N_NS, with_NSs, config['EoS'])
                NSt = config.get('NSt', 0)
                if NSt == 0:
                    tNS = t + tBH_form * np.ones(N_NS)
                else:
                    tNS = t + sample_ns_formation_times(mNS, NSt, tBH_form, table=NSFormationTable())
                m_new.append(mNS); s_new.append(np.zeros(N_NS)); t_new.append(tNS); v_new.append(vNS)
                n_new.append(np.ones(N_NS))

        m_all = np.concatenate([state['mCO_pending']] + m_new)
        s_all = np.concatenate([state['sCO_pending']] + s_new)
        t_all = np.concatenate([state['tCO_pending']] + t_new)
        v_all = np.concatenate([state['vCO_pending']] + v_new)
        n_all = np.concatenate([state['nCO_pending']] + n_new)
        order = np.argsort(t_all, kind='stable')
        state['mCO_pending'] = m_all[order]; state['sCO_pending'] = s_all[order]
        state['tCO_pending'] = t_all[order]; state['vCO_pending'] = v_all[order]
        state['nCO_pending'] = n_all[order]

        print('STELLAR POPULATION {} BORN AT t = {:.6g} Myr: N = {:.6g}, M = {:.6g} Msun, Z = {:.6g}, SFE = {:.6g}, '
              'r_h,k/r_h = {:.6g}{}'.format(pops.size - 1, t, N_k, M_k, Z_k, SFE_k, p['x'],
                                           '' if p['x'] == 1.0 else ', mixing time {:.6g} t_rh'.format(p['n_mix'])))


def form_pending_compact_objects(state, config):
    """Form the BHs and NSs of later stellar populations whose formation time has been
    reached; each is retained if its natal kick is below the current escape velocity, and
    joins the single compact objects (generation 1, no TDEs). Updates state in-place."""
    tq = state['tCO_pending']
    if tq.size == 0:
        return
    k = np.searchsorted(tq, state['t'], side='right')
    if k == 0:
        return
    m = state['mCO_pending'][:k]; s = state['sCO_pending'][:k]; v = state['vCO_pending'][:k]
    retained = v < v_esc(state['Mcl'] + state['M_gas'], state['rh'])
    n_ret = int(retained.sum())
    if n_ret > 0:
        state['BH'].add_many(m[retained], s[retained], np.ones(n_ret), np.zeros(n_ret))
        if state['i_aux1'] == 1:
            state['N_BH'] = state['N_BH'] + n_ret
    is_NS = state['nCO_pending'][:k] == 1
    state['mCO_pending'] = state['mCO_pending'][k:]; state['sCO_pending'] = state['sCO_pending'][k:]
    state['tCO_pending'] = tq[k:]; state['vCO_pending'] = state['vCO_pending'][k:]
    state['nCO_pending'] = state['nCO_pending'][k:]
    state['N_pop_CO_formed'] = state['N_pop_CO_formed'] + k
    state['N_pop_CO_retained'] = state['N_pop_CO_retained'] + n_ret
    # (NSs are also counted with those of the first population, printed at the end for -NSt > 0)
    state['N_NS_formed'] = state['N_NS_formed'] + int(is_NS.sum())
    state['N_NS_retained'] = state['N_NS_retained'] + int((is_NS & retained).sum())


def form_neutron_stars(state, config):
    """Form the NSs whose formation time has been reached (NSt > 0 only).

    NSs with t_form <= t are taken from the front of the time-sorted pending
    arrays. Each is retained if its natal kick is below the escape velocity at
    the current time, and retained NSs join the compact-object arrays with zero
    spin, generation 1 and no TDEs. If the BH subsystem is already active, they
    are also added to N_BH; otherwise they are counted when it activates.
    Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.
    """
    tNS_pending = state['tNS_pending']
    if tNS_pending.size == 0:
        return

    # number of pending NSs that have formed by now (arrays sorted by t_form):
    k = np.searchsorted(tNS_pending, state['t'], side='right')
    if k == 0:
        return

    mNS_new = state['mNS_pending'][:k]
    vNS_new = state['vNS_pending'][:k]

    # retain NSs with natal kick below the current escape velocity:
    retained = vNS_new < v_esc(state['Mcl'] + state['M_gas'], state['rh'])
    mNS_ret = mNS_new[retained]
    N_ret = mNS_ret.size

    if N_ret > 0:
        state['BH'].add_many(mNS_ret, np.zeros(N_ret), np.ones(N_ret), np.zeros(N_ret))
        if state['i_aux1'] == 1:
            state['N_BH'] = state['N_BH'] + N_ret

    # remove the formed NSs from the pending arrays:
    state['mNS_pending'] = state['mNS_pending'][k:]
    state['tNS_pending'] = tNS_pending[k:]
    state['vNS_pending'] = state['vNS_pending'][k:]

    state['N_NS_formed'] = state['N_NS_formed'] + k
    state['N_NS_retained'] = state['N_NS_retained'] + N_ret


def compute_cluster_properties(state, config):
    """Compute BH subsystem and cluster structural properties for the current timestep.

    Saves BH mass/spin/generation snapshots, computes the average and maximum BH
    mass, Spitzer mass-segregation parameter, BH half-mass and core radii, velocity
    dispersions, number densities, relaxation timescales, and the cluster Coulomb
    logarithm. Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.

    Returns:
        bool: True if the simulation should continue, False if the BH subsystem
            has evaporated (N_BH <= 0).
    """

    # unpack current state into local variables:
    t = state['t']
    BH = state['BH']
    binaries = state['binaries']
    pairs = state['pairs']
    i_aux1 = state['i_aux1']
    N_BH = state['N_BH']
    N_BBH = state['N_BBH']
    N_BHstar = state['N_BHstar']
    N_Triples = state['N_Triples']
    Mcl = state['Mcl']
    rh = state['rh']
    m_avg = state['m_avg']
    m_avg0 = state['m_avg0']
    M_gas = state['M_gas']
    Nb = state['Nb']
    ab = state['ab']
    n_star0 = state['n_star0']
    rh0 = state['rh0']
    Mcl0 = state['Mcl0']

    # save current BH properties (only needed for the BH-masses output file, -BOi 1;
    # for clusters with millions of BHs these snapshots otherwise grow to tens of GB):
    # (row 0 of binaries and pairs is an all-zero placeholder, not a BH, and is skipped)
    if config['BOi']==1:
        state['simulation_times'].append(t)
        state['black_hole_masses'].append(np.concatenate((BH.m, binaries[1:, 4], binaries[1:, 5], pairs[1:, 1])))
        state['black_hole_spins'].append(np.concatenate((BH.s, binaries[1:, 6], binaries[1:, 7], pairs[1:, 2])))
        state['black_hole_generations'].append(np.concatenate((BH.g, binaries[1:, 8], binaries[1:, 9], pairs[1:, 3])))
        state['black_hole_tdes'].append(np.concatenate((BH.h, binaries[1:, 13], binaries[1:, 14], pairs[1:, 4])))

    # NSs still to form (only for NSt > 0), and compact objects of later populations (still
    # to form, or of populations to be born before t_max):
    pops = state['pops']
    NSs_pending = (state['tNS_pending'].size > 0 or state['tCO_pending'].size > 0
                   or pops.next_birth_time() < config['t_max'])

    # several stellar populations: draw field stars from them (stellar_evolution.get_star)
    if pops.multi:
        set_field_star_mixture(pops.ages(t), pops.central_weights(t, m_avg0))
    else:
        set_field_star_mixture(None, None)

    # activate BH subsystem once BH formation time is reached (if it would be
    # empty and NSs are still to form, wait for them):
    if t>tBH_form and i_aux1==0 and (BH.size>0 or not NSs_pending):
        N_BH = BH.size
        i_aux1 = 1

    # check if BH subsystem has fully evaporated:
    if i_aux1==1 and N_BH<=0:
        if NSs_pending:
            # deactivate until the next NSs form:
            i_aux1 = 0
        else:
            print('BH SUBSYSTEM EVAPORATED')
            state['i_aux1'] = i_aux1
            state['N_BH'] = N_BH
            return False

    # average BH mass:
    if i_aux1==1:
        mBH_avg = (BH.total(FAM['mass']) + np.sum(binaries[:, 4]) + np.sum(binaries[:, 5]) + np.sum(pairs[:, 1]) ) / N_BH
    else:
        mBH_avg = 0.0

    # maximum BH mass:
    if i_aux1==1:
        mBHs_max = BH.max_mass()
        try:
            mBH1_max = np.max(binaries[:, 4])
        except Exception:
            mBH1_max = 0.0
        try:
            mBH2_max = np.max(binaries[:, 5])
        except Exception:
            mBH2_max = 0.0
        try:
            mBHp_max = np.max(pairs[:, 1])
        except Exception:
            mBHp_max = 0.0
        mBH_max = np.max([mBHs_max, mBH1_max, mBH2_max, mBHp_max])
    else:
        mBH_max = 0.0

    # total BH mass ratio:
    Q_BH = N_BH * mBH_avg / Mcl

    # individual BH mass ratio:
    q_BH = mBH_avg / m_avg

    # Cluster Coulomb logarithm, Lambda = 0.02 N (Breen & Heggie 2013), smoothly floored at 1;
    # N = M_cl / <m> is the current number of stars:
    logLcl = np.log(np.e + lc * Mcl / m_avg)

    # BH Coulomb logarithm, Lambda_BH = 0.02 N_BH (Breen & Heggie 2013), smoothly
    # floored at 1 (their small-N_BH limit): ln(e + 0.02 N_BH) -> 1 as N_BH -> 0
    # and -> ln(0.02 N_BH) for N_BH >> 136:
    if i_aux1==1 and N_BH>0:
        logLBH = np.log(np.e + lc * N_BH)
    else:
        logLBH = 0.0

    # Spitzer parameter (Spitzer 1969); stars and BHs can reach
    # energy equipartition only if S < 0.16:
    S = Q_BH * q_BH**(3/2)

    # BH-to-star temperature ratio (xi = 1: equipartition; Breen & Heggie 2013):
    if S < 0.16 and S > 0:
        xi = 1
    else:
        if i_aux1==1:
            # c_xi: order-unity prefactor (Breen & Heggie 2013 give only the scaling); set by -cxi
            xi = config['c_xi'] * q_BH**(3/5) * Q_BH**(2/5) * (logLBH / logLcl)**(-2/5)
        else:
            xi = 0.0

    # multimass relaxation factor, N-body calibration of Antonini & Gieles (2020, Eq. 12):
    # psi = 1 + a1 f_BH / 0.01, with f_BH = M_BH / M_cl = Q_BH.
    # (Replaces psi = 1 + S, which assumes equipartition and overestimates psi
    # by a factor ~3 once S > 0.16 and the BHs decouple from the stars.)
    psi = 1 + a1_psi * Q_BH / 0.01

    # BH relaxation factor:
    if i_aux1==1 and N_BH>0:
        psi_BH = (BH.total(FAM['p25']) + np.sum(binaries[:, 4]**(5/2)) + \
                  np.sum(binaries[:, 5]**(5/2)) + np.sum(pairs[:, 1]**(5/2)) ) \
                  / N_BH / mBH_avg**(5/2)
    else:
        psi_BH = 0.0

    # BH half-mass radius:
    if i_aux1==1:
        rh_BH = 1/xi * Q_BH * q_BH * rh
    else:
        rh_BH = 0.0

    # BH core radius:
    if i_aux1==1:
        Core_parameter = 1.76 / Hardening_constant / psi_BH * (12 * q_BH / xi - 1) * eta_min**(-11/2) * (1+3*eta_min) * (1 + 6*eta_min) * P_3bb
        rc_BH = rh_BH * N_BH**(-2/3) * (Core_parameter / zeta_BH / logLBH)**(1/3)
    else:
        rc_BH = 0.0

    # BH velocity dispersion:
    if i_aux1==1:
        vBH = np.sqrt(0.4 * G_Newton * mBH_avg * N_BH / rh_BH)
    else:
        vBH = 0.0

    # star velocity dispersion:
    v_star = np.sqrt(0.4 * G_Newton * Mcl / rh)

    # BH half-mass volume:
    Vh_BH = 4 / 3 * np.pi * rh_BH**3

    # BH half-mass density:
    if i_aux1==1:
        nh_BH = N_BH / 2 / Vh_BH
    else:
        nh_BH = 0.0

    # BH core density:
    if i_aux1==1:
        nc_BH = nh_BH * (rc_BH / rh_BH)**(-2)
    else:
        nc_BH = 0.0

    # BH core volume:
    Vc_BH = 4 / 3 * np.pi * rc_BH**3

    # mean BH density:
    if i_aux1==1:
        na_BH = nc_BH * 2/3
    else:
        na_BH = 0.0

    # stellar density updated:
    n_star = n_star0 * (rh0 / rh)**3 * (Mcl / Mcl0) * (m_avg0 / m_avg)
    if pops.multi and pops.segregated:
        # central density of populations with their own half-mass radii (populations.py):
        n_star = n_star * pops.density_factor(t, m_avg0)

    # half-mass volume:
    Vh = 4 * np.pi / 3 * rh**3

    if Nb>0:
        nb = ab.size / Vh
    else:
        nb = 0.0

    # relaxation timescale:
    t_rlx = t_relax(Mcl + M_gas, rh, m_avg, psi, logLcl)

    # BH relaxation timescale:
    if i_aux1==1:
        tBH_rlx = t_relax(N_BH*mBH_avg, rh_BH, mBH_avg, psi_BH, logLBH)
    else:
        tBH_rlx = 1e100

    # BH subsystem core collapse:
    if i_aux1==1:
        tBH_cc = k_cc * tBH_rlx
    else:
        tBH_cc = 1e100

    # Number of BHs in the core:
    Nc_BH = nc_BH * Vc_BH

    # update state:
    state['i_aux1'] = i_aux1
    state['N_BH'] = N_BH
    state['mBH_avg'] = mBH_avg
    state['mBH_max'] = mBH_max
    state['Q_BH'] = Q_BH
    state['q_BH'] = q_BH
    state['S'] = S
    state['xi'] = xi
    state['psi'] = psi
    state['psi_BH'] = psi_BH
    state['rh_BH'] = rh_BH
    state['rc_BH'] = rc_BH
    state['vBH'] = vBH
    state['v_star'] = v_star
    state['nh_BH'] = nh_BH
    state['nc_BH'] = nc_BH
    state['na_BH'] = na_BH
    state['Vc_BH'] = Vc_BH
    state['Nc_BH'] = Nc_BH
    state['n_star'] = n_star
    state['nb'] = nb
    state['t_rlx'] = t_rlx
    state['tBH_rlx'] = tBH_rlx
    state['logLcl'] = logLcl
    state['logLBH'] = logLBH

    return True


# np.mean of the star-star binary separations, recomputed only when `ab` is a new
# array (it is replaced, never modified in place). The array has up to N*fb/2
# entries (~10^6 for heavy clusters), and its mean was needed at every step.
_ab_mean_cache = {'ab': None, 'mean': None}

def _mean_ab(ab):
    if _ab_mean_cache['ab'] is not ab:
        _ab_mean_cache['ab'] = ab
        _ab_mean_cache['mean'] = np.mean(ab)
    return _ab_mean_cache['mean']


def compute_timescales(state, config):
    """Compute dynamical timescales and the adaptive timestep.

    Calculates the three-body binary formation (t_3bb), two-body gravitational
    capture (t_2cap), star-star to BH-star exchange (t_ex1), and BH-star to BBH
    exchange (t_ex2) timescales. Sets the adaptive timestep dt as the minimum of
    these, bounded by [dt_min, dt_max]. Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.
    """

    # unpack current state into local variables:
    t = state['t']
    t_cc = state['t_cc']
    N_BH = state['N_BH']
    N_BBH = state['N_BBH']
    N_BHstar = state['N_BHstar']
    N_Triples = state['N_Triples']
    mBH_avg = state['mBH_avg']
    nc_BH = state['nc_BH']
    vBH = state['vBH']
    Vc_BH = state['Vc_BH']
    Nc_BH = state['Nc_BH']
    nb = state['nb']
    v_star = state['v_star']
    m_avg = state['m_avg']
    binaries = state['binaries']
    pairs = state['pairs']
    ab = state['ab']
    dt_min = config['dt_min']
    dt_max = config['dt_max']

    # number of single (unbound) BHs available for interactions:
    N_BHsin = N_BH - 2*N_BBH - N_BHstar - 3*N_Triples

    # 3bb timescale:
    if t>t_cc and N_BHsin>2:
        try:
            t_3bb = 1 / Rate_3bb(mBH_avg, nc_BH, vBH) / Vc_BH
        except Exception:
            t_3bb = 1e100
    else:
        t_3bb = 1e100

    # 2-body capture timescale:
    if t>t_cc and N_BHsin>1:
        try:
            t_2cap = 1 / Rate_cap(mBH_avg, nc_BH, vBH) / Vc_BH
        except Exception:
            t_2cap = 1e100
    else:
        t_2cap = 1e100

    # star-star -> BH-star timescale:
    if t>t_cc and N_BHsin>0 and nb>0:
        try:
            t_ex1 = 1 / Rate_exc(m_avg, m_avg, mBH_avg, nb, np.sqrt(v_star**2 + vBH**2), _mean_ab(ab)) / Nc_BH / 2
        except Exception:
            t_ex1 = 1e100
    else:
        t_ex1 = 1e100

    # BH-star -> BH-BH timescale:
    if t>t_cc and N_BHstar>0 and N_BHsin>0:
        if N_BHstar>0:
            rate_ex2 = Rate_exc(m_avg, np.mean(pairs[:, 1]), mBH_avg, nc_BH, vBH, np.mean(pairs[:, 0]))
            t_ex2 = 1 / rate_ex2 / N_BHstar if rate_ex2 > 0 else 1e100
        else:
            t_ex2 = 1e100
    else:
        t_ex2 = 1e100

    # adaptive time-step:
    if t<t_cc:
        dt = dt_min
    else:
        dt = np.min([dt_max, np.max([dt_min, np.min([t_3bb, t_2cap, t_ex1, t_ex2])])])

    # at most a fraction f_dt_esc of the cluster mass lost by escape in a step (escape rate
    # of the previous step): otherwise a tidally filling cluster, with a large escape rate
    # xi_e, could lose most of its mass and reach a negative half-mass radius in one step
    xi_e_prev = state.get('xi_e', 0.0)
    t_rlx = state['t_rlx']
    if xi_e_prev > 0 and t_rlx > 0:
        dt_esc = f_dt_esc * t_rlx / xi_e_prev
        if dt > dt_esc:
            dt = max(dt_min, dt_esc)

    if dt > t and t>0:
        dt = t

    # several stellar populations: steps of at most the age of the youngest one (as dt <= t
    # above for one population: stellar-evolution mass loss ~ dt/age), and steps that end at
    # the birth of the next population:
    pops = state['pops']
    if pops.multi:
        youngest = min(t - tk for tk in pops.t_birth)
        dt = min(dt, max(dt_min, youngest))
    t_next = pops.next_birth_time()
    if t < t_next < math.inf:
        dt = min(dt, t_next - t)

    state['t_3bb'] = t_3bb
    state['t_2cap'] = t_2cap
    state['t_ex1'] = t_ex1
    state['t_ex2'] = t_ex2
    state['dt'] = dt


def form_binaries(state, config):
    """Form new binary BHs through three-body, two-body capture, and exchange channels.

    Draws Poisson-sampled event counts for each formation channel and executes:
      - Three-body binary (3bb) formation (channel 3)
      - Two-body gravitational wave capture (channel 2)
      - Star-star -> BH-star exchanges (channel 1, step 1)
      - BH-star -> BBH exchanges (channel 1, step 2)
    Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.
    """
    
    # unpack current state into local variables:
    t = state['t']
    z = state['z']
    dt = state['dt']
    zCl_form = state['zCl_form']
    mBH_avg = state['mBH_avg']
    binaries = state['binaries']
    BH = state['BH']
    vBH = state['vBH']
    v_star = state['v_star']
    m_avg = state['m_avg']
    ab = state['ab']
    pairs = state['pairs']
    mergers = state['mergers']
    seed = state['seed']
    nc_BH = state['nc_BH']
    N_BH = state['N_BH']
    N_BBH = state['N_BBH']
    N_BHstar = state['N_BHstar']
    N_Triples = state['N_Triples']
    N_3bb = state['N_3bb']
    N_2cap = state['N_2cap']
    N_me = state['N_me']
    N_meRe = state['N_meRe']
    N_meEj = state['N_meEj']
    N_ex1 = state['N_ex1']
    N_ex2 = state['N_ex2']
    t_3bb = state['t_3bb']
    t_2cap = state['t_2cap']
    N_tdeBHstar = state['N_tdeBHstar']
    tdes = state['tdes']
    m_min = config['m_min']
    m_max = config['m_max']
    with_tdes = config['with_tdes']
    f_accreted = config['f_accreted']
    EoS = config['EoS']

    # number of single (unbound) BHs available for interactions:
    N_BHsin = N_BH - 2*N_BBH - N_BHstar - 3*N_Triples

    # number of 3bbs:
    k_3bb = np.min([poisson_draw(np.min([dt / t_3bb, 1e8])), int(N_BHsin / 3)])

    # 3bb formation:
    t, z, k_3bb, mBH_avg, binaries, BH, vBH, N_3bb, N_BBH = three_body_binary(t, z, k_3bb, mBH_avg, binaries, BH, vBH, N_3bb, N_BBH, random_pairing=config['random_pairing'])

    # number of 2-body captures:
    N_BHsin = N_BH - 2*N_BBH - N_BHstar - 3*N_Triples
    k_2cap = np.min([poisson_draw(np.min([dt / t_2cap, 1e8])), int(N_BHsin / 2)])

    # 2-body capture(s):
    seed, t, dt, z, zCl_form, k_2cap, mBH_avg, binaries, BH, vBH, v_star, N_2cap, N_BH, N_BBH, N_me, N_meRe, N_meEj, mergers = \
        two_body_capture(seed, t, dt, z, zCl_form, k_2cap, mBH_avg, binaries, BH, vBH, v_star, N_2cap, N_BH, N_BBH, N_me, N_meRe, N_meEj, mergers, random_pairing=config['random_pairing'])

    # number of star-star -> BH-star exchanges:
    N_BHsin = N_BH - 2*N_BBH - N_BHstar - 3*N_Triples
    k_ex1 = np.min([poisson_draw(np.min([dt / state['t_ex1'], 1e8])), int(N_BHsin)])

    # star-star -> BH-star exchange(s):
    if k_ex1 > 0:
        seed, t, z, k_ex1, N_ex1, m_avg, BH, ab, pairs, N_BHstar, N_tdeBHstar, v_star, vBH, tdes, binaries, m_min, m_max, with_tdes, f_accreted, EoS = StarStar_to_BHstar(seed, t, z, k_ex1, N_ex1, m_avg, BH, ab, pairs, N_BHstar, N_tdeBHstar, v_star, vBH, tdes, binaries, m_min, m_max, with_tdes, f_accreted, EoS)
    
    # number of BH-star -> BH-BH exchanges:
    N_BHsin = N_BH - 2*N_BBH - N_BHstar - 3*N_Triples
    k_ex2 = np.min([poisson_draw(np.min([dt / state['t_ex2'], 1e8])), int(N_BHsin), int(N_BHstar)])

    # BH-star -> BBH exchange(s):
    if k_ex2 > 0:
        seed, t, z, k_ex2, N_ex2, m_avg, BH, pairs, binaries, N_BBH, N_BHstar, N_tdeBHstar, v_star, vBH, tdes, m_min, m_max, with_tdes, f_accreted, EoS = BHstar_to_BBH(seed, t, z, k_ex2, N_ex2, m_avg, BH, pairs, binaries, N_BBH, N_BHstar, N_tdeBHstar, v_star, vBH, tdes, m_min, m_max, with_tdes, f_accreted, EoS)
    
    # write back:
    state['t'] = t; state['z'] = z; state['dt'] = dt; state['seed'] = seed
    state['mBH_avg'] = mBH_avg; state['binaries'] = binaries
    state['vBH'] = vBH; state['v_star'] = v_star
    state['m_avg'] = m_avg; state['ab'] = ab; state['pairs'] = pairs; state['mergers'] = mergers
    state['N_BH'] = N_BH; state['N_BBH'] = N_BBH; state['N_BHstar'] = N_BHstar
    state['N_3bb'] = N_3bb; state['N_2cap'] = N_2cap; state['N_me'] = N_me
    state['N_meRe'] = N_meRe; state['N_meEj'] = N_meEj
    state['N_ex1'] = N_ex1; state['N_ex2'] = N_ex2
    state['k_3bb'] = k_3bb; state['k_2cap'] = k_2cap
    state['k_ex1'] = k_ex1; state['k_ex2'] = k_ex2
    state['zCl_form'] = zCl_form
    state['N_tdeBHstar'] = N_tdeBHstar; state['tdes'] = tdes


def evolve_interactions(state, config):
    """Evolve BBH binaries, pair-pair interactions, and hierarchical triples.

    Calls evolve_BBHs to harden/disrupt/merge existing BBHs through single-binary
    and binary-binary encounters, processes BH-star pair-pair interactions that can
    form new BBHs, and evolves hierarchical triples via von Zeipel-Lidov-Kozai
    (ZLK) oscillations. Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.
    """
    
    # unpack current state into local variables:
    seed = state['seed']
    t = state['t']; z = state['z']; dt = state['dt']; zCl_form = state['zCl_form']
    binaries = state['binaries']; hardening = state['hardening']; mergers = state['mergers']
    BH = state['BH']
    n_star = state['n_star']; v_star = state['v_star']; vBH = state['vBH']
    t_rlx = state['t_rlx']; m_avg = state['m_avg']; mBH_avg = state['mBH_avg']
    na_BH = state['na_BH']; nc_BH = state['nc_BH']
    N_BH = state['N_BH']; N_BBH = state['N_BBH']; N_me = state['N_me']
    N_me2b = state['N_me2b']; N_3cap = state['N_3cap']; N_meFi = state['N_meFi']
    N_meRe = state['N_meRe']; N_meEj = state['N_meEj']
    N_dis = state['N_dis']; N_ex = state['N_ex']
    N_BHej = state['N_BHej']; N_BBHej = state['N_BBHej']
    N_hardening = state['N_hardening']; Vc_BH = state['Vc_BH']
    N_bb = state['N_bb']; triples = state['triples']; N_Triples = state['N_Triples']
    pairs = state['pairs']; N_BHstar = state['N_BHstar']
    N_pp = state['N_pp']; N_ZLK = state['N_ZLK']
    i_aux1 = state['i_aux1']
    t_cc = state['t_cc']
    tdes = state['tdes']
    N_tdeBHstar = state['N_tdeBHstar']
    N_tdeBBHstar = state['N_tdeBBHstar']
    m_min = config['m_min']
    m_max = config['m_max']
    with_tdes = config['with_tdes']
    f_accreted = config['f_accreted']
    EoS = config['EoS']

    # BBH evolution:
    seed, t, z, dt, zCl_form, binaries, hardening, mergers, BH, n_star, v_star, vBH, t_rlx, m_avg, mBH_avg, na_BH, nc_BH, N_BH, N_BBH, N_me, N_me2b, N_3cap, N_meFi, N_meRe, N_meEj, N_dis, N_ex, N_BHej, N_BBHej, N_hardening, Vc_BH, N_bb, triples, N_Triples, tdes, N_tdeBBHstar, m_min, m_max, f_accreted, EoS = evolve_BBHs(seed, t, z, dt, zCl_form, binaries, hardening, mergers, BH, n_star, v_star, vBH, t_rlx, m_avg, mBH_avg, na_BH, nc_BH, N_BH, N_BBH, N_me, N_me2b, N_3cap, N_meFi, N_meRe, N_meEj, N_dis, N_ex, N_BHej, N_BBHej, N_hardening, Vc_BH, N_bb, triples, N_Triples, tdes, N_tdeBBHstar, m_min, m_max, f_accreted, EoS)

    # average BBH number density:
    if i_aux1==1:
        n_BBH = N_BBH / Vc_BH
    else:
        n_BBH = 0.0

    # average BH-star number density:
    if i_aux1==1:
        n_BHstar = N_BHstar / Vc_BH
    else:
         n_BHstar = 0.0

    # Binary-binary interaction timescale:
    if t>t_cc and N_BBH>1:
        try:
            t_bb = 1 / Rate_int(np.mean(binaries[:, 4]+binaries[:, 5]), n_BBH, vBH, 2 * np.mean(binaries[:, 2])) / N_BBH
        except Exception:
            t_bb = 1e100
    else:
        t_bb = 1e100

    # Pair-pair interaction timescale:
    if t>t_cc and N_BHstar>1:
        try:
            t_pp = 1 / Rate_int(np.mean(pairs[:, 1]) + m_avg, n_BHstar, vBH, 2 * np.mean(pairs[:, 0])) / N_BHstar
        except Exception:
            t_pp = 1e100
    else:
        t_pp = 1e100

    # number of pair-pair interactions:
    dt = state['dt']
    
    # cap mu to avoid scipy Poisson overflow for extremely high interaction rates:
    k_pp = np.min([poisson_draw(np.min([dt / t_pp, 1e8])), int(N_BHstar / 2)])

    # execute pair-pair interactions: two BH-star pairs form a new BBH:
    if k_pp>0:
        for i in range(k_pp):
            N_pp+=1

            smas = pairs[:, 0]
            masses = pairs[:, 1]

            valid = np.where(smas * masses > 0)[0]
            if len(valid) < 2:
                break
            a1, a2 = np.random.choice(smas[valid], size=2, replace=False, p=(smas*masses)[valid] / np.sum((smas*masses)[valid]))

            k1 = np.squeeze(np.where(smas==a1))+0
            k2 = np.squeeze(np.where(smas==a2))+0

            k1 = int(np.atleast_1d(k1)[0])
            k2 = int(np.atleast_1d(k2)[0])

            if k1 == k2:
                candidates = np.where(smas == a2)[0]
                k2 = int(candidates[1]) if len(candidates) > 1 else None
                if k2 is None:
                    continue

            m1 = pairs[k1][1]; s1 = pairs[k1][2]; g1 = pairs[k1][3]; h1 = pairs[k1][4]
            m2 = pairs[k2][1]; s2 = pairs[k2][2]; g2 = pairs[k2][3]; h2 = pairs[k2][4]

            pairs = np.delete(pairs, [k1, k2], axis=0)
            N_BHstar = N_BHstar - 2

            # probabilities for TDE:
            p1_TDE = min(R_sun*(m1/m_avg)**(1/3)/(a1 + a2), 1.0) if with_tdes==1 else 0.0
            p2_TDE = min(R_sun*(m2/m_avg)**(1/3)/(a1 + a2), 1.0) if with_tdes==1 else 0.0

            if np.random.rand() < p1_TDE:
                # then TDE-1 occurs:

                tde_type = 4 # TDE during BH-star+BH-star encounter

                # Sample star mass from evolving mass function:
                m_star, R_star = get_star(t, tBH_form, m_min, m_max)

                N_tdeBHstar, m1, s1, g1, h1, tdes = single_BH_TidalDisruption(seed, t, z, N_tdeBHstar, tde_type, m_avg, m_star, R_star, m1, s1, g1, h1, v_star, vBH, tdes, f_accreted, EoS)

            if np.random.rand() < p2_TDE:
                # then TDE-2 occurs:

                tde_type = 4 # TDE during BH-star+BH-star encounter

                # Sample star mass from evolving mass function:
                m_star, R_star = get_star(t, tBH_form, m_min, m_max)

                N_tdeBHstar, m2, s2, g2, h2, tdes = single_BH_TidalDisruption(seed, t, z, N_tdeBHstar, tde_type, m_avg, m_star, R_star, m2, s2, g2, h2, v_star, vBH, tdes, f_accreted, EoS)

            # the following exchange happens:
            # BH-star + BH-star -> BH-BH + star' + star'

            if m2/a2 < m1/a1:
                sma = m2 / m_avg * a1
            else:
                sma = m1 / m_avg * a2

            eccen = np.sqrt(np.random.rand()) # thermal eccentricity

            # append new BBH:
            binaries = np.append(binaries, [[np.random.randint(0, 999999999), 1, sma, eccen, m1, m2, s1, s2, g1, g2, t, z, 0, h1, h2]], axis=0)
            N_BBH+=1

    # Triple evolution:
    seed, t, z, zCl_form, triples, binaries, BH, mBH_avg, N_Triples, N_BBH, N_BH, N_me, N_meRe, N_meEj, N_ZLK, v_star, vBH, nc_BH, mergers = \
        evolve_triples(seed, t, z, zCl_form, triples, binaries, BH, mBH_avg, N_Triples, N_BBH, N_BH, N_me, N_meRe, N_meEj, N_ZLK, v_star, vBH, nc_BH, mergers)

    # write back:
    state['seed'] = seed; state['t'] = t; state['z'] = z; state['dt'] = dt; state['zCl_form'] = zCl_form
    state['binaries'] = binaries; state['hardening'] = hardening; state['mergers'] = mergers
    state['n_star'] = n_star; state['v_star'] = v_star; state['vBH'] = vBH
    state['t_rlx'] = t_rlx; state['m_avg'] = m_avg; state['mBH_avg'] = mBH_avg
    state['na_BH'] = na_BH; state['nc_BH'] = nc_BH
    state['N_BH'] = N_BH; state['N_BBH'] = N_BBH; state['N_me'] = N_me
    state['N_me2b'] = N_me2b; state['N_3cap'] = N_3cap; state['N_meFi'] = N_meFi
    state['N_meRe'] = N_meRe; state['N_meEj'] = N_meEj
    state['N_dis'] = N_dis; state['N_ex'] = N_ex
    state['N_BHej'] = N_BHej; state['N_BBHej'] = N_BBHej
    state['N_hardening'] = N_hardening; state['Vc_BH'] = Vc_BH
    state['N_bb'] = N_bb; state['triples'] = triples; state['N_Triples'] = N_Triples
    state['pairs'] = pairs; state['N_BHstar'] = N_BHstar
    state['N_pp'] = N_pp; state['N_ZLK'] = N_ZLK
    state['t_bb'] = t_bb; state['t_pp'] = t_pp; state['k_pp'] = k_pp
    state['tdes'] = tdes
    state['N_tdeBHstar'] = N_tdeBHstar
    state['N_tdeBBHstar'] = N_tdeBBHstar


def evolve_tdes(state, config):
    """Evolve white dwarf population and process tidal disruption events.

    Computes WD formation, evaporation, and TDE rates. Draws Poisson-sampled
    BH-WD and BH-star TDE counts, executes each TDE (updating BH mass and spin
    via accretion), and records event parameters. Skipped if config['with_tdes']
    is False. Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.
    """
    
    # skip TDEs if disabled; zero out all TDE-related state:
    if not config['with_tdes']:
        state['N_WD'] = 0; state['v_WD'] = 0
        state['k_tdeBHWD'] = 0; state['N_tdeBHWD'] = state.get('N_tdeBHWD', 0)
        state['dN_WDformdt'] = 0; state['dN_WDevdt'] = 0; state['dN_tdeBHWDdt'] = 0
        state['k_tdeBHstar'] = 0; state['dN_tdeBHstardt'] = 0
        state['N_tdeBHstar'] = state.get('N_tdeBHstar', 0)
        return

    # unpack current state into local variables:
    seed = state['seed']; t = state['t']; z = state['z']; dt = state['dt']
    BH = state['BH']
    vBH = state['vBH']; v_star = state['v_star']
    mBH_avg = state['mBH_avg']; m_avg = state['m_avg']
    N_BH = state['N_BH']; N_Triples = state['N_Triples']
    N_WD = state['N_WD']; N_tdeBHWD = state['N_tdeBHWD']; N_tdeBHstar = state['N_tdeBHstar']
    n_star = state['n_star']; t_rlx = state['t_rlx']
    Mcl = state['Mcl']; rh = state['rh']
    binaries = state['binaries']; pairs = state['pairs']; tdes = state['tdes']
    Kroupa_norm = state['Kroupa_norm']
    xi_e = state['xi_e']
    m_min = config['m_min']
    m_max = config['m_max']
    f_accreted = config['f_accreted']
    EoS = config['EoS']

    # white dwarf formation parameters (solar lifetime and max WD progenitor mass):
    solar_life = 1.0e4
    m_WN = 8.0
    t_WN = solar_life/m_WN**(2.5)
    N_strs = Mcl/m_avg

    # WD formation rate (/Myr):
    if not state['pops'].multi:
        dN_WDformdt = N_strs/t/2.5 * Kroupa_norm*IMF_kroupa((solar_life/t)**(1/2.5)) * (solar_life/t)**(1/2.5) if t>t_WN else 0.0
    else:
        # sum over the stellar populations, each at its own age:
        dN_WDformdt = state['pops'].white_dwarf_formation_rate(t, state['m_avg0'], Kroupa_norm, IMF_kroupa)

    # WD evaporation rate (/Myr):
    dN_WDevdt = xi_e*N_WD/t_rlx if N_WD>0 else 0.0

    # white dwarf properties and BH-WD TDE rate:
    m_WD = white_dwarf_mass
    R_WD = R_WhiteDwarf()
    v_WD = np.sqrt(m_avg/m_WD)*v_star
    rh_WD = rh
    nc_WD = 3*N_WD/4/np.pi/(rh_WD/1.3)**3
    v_BHWD = np.sqrt(vBH**2 + v_WD**2)
    dN_tdeBHWDdt = 2*np.sqrt(2*(3*np.pi-8))*G_Newton*(mBH_avg + m_WD)*N_BH*nc_WD*R_WD*(mBH_avg/m_WD)**(1/3)/v_BHWD if N_WD>0 else 0.0

    # effective WD evolution equation:
    dN_WDdt = dN_WDformdt - dN_WDevdt - dN_tdeBHWDdt
    N_WD = N_WD + dN_WDdt * dt

    # Poisson number of BH-WD TDEs:
    k_tdeBHWD = np.min([poisson_draw(np.min([dt*dN_tdeBHWDdt, 1e8])), int(N_BH-3*N_Triples), int(N_WD)]) if N_WD>0 else 0

    # execute BH-WD tidal disruption events:
    if k_tdeBHWD > 0:
        tde_type = 11
        seed, t, z, k_tdeBHWD, N_tdeBHWD, tde_type, m_avg, m_WD, R_WD, BH, v_WD, vBH, tdes, binaries, pairs, f_accreted, EoS = BH_TidalDisruptions(seed, t, z, k_tdeBHWD, N_tdeBHWD, tde_type, m_avg, m_WD, R_WD, BH, v_WD, vBH, tdes, binaries, pairs, f_accreted, EoS)

    # micro-TDEs:
    v_BHstar = np.sqrt(v_star**2 + vBH**2)
    dN_tdeBHstardt = 2*np.sqrt(2*(3*np.pi-8))*G_Newton*(mBH_avg + m_avg)*n_star*N_BH*R_sun*(mBH_avg/m_avg)**(1/3)/v_BHstar if (N_BH>0)*(n_star>0) else 0.0

    # Poisson number of BH-star TDEs at this timestep:
    N_strs = Mcl/m_avg
    k_tdeBHstar = np.min([poisson_draw(np.min([dt*dN_tdeBHstardt, 1e8])), int(N_BH-3*N_Triples), int(N_strs)]) if (N_BH>0)*(n_star>0) else 0

    # execute BH-star tidal disruption events (micro-TDEs):
    if k_tdeBHstar > 0:
        tde_type = 1

        # Sample star mass from evolving mass function:
        m_star, R_star = get_star(t, tBH_form, m_min, m_max)

        seed, t, z, k_tdeBHstar, N_tdeBHstar, tde_type, m_avg, m_star, R_star, BH, v_star, vBH, tdes, binaries, pairs, f_accreted, EoS = BH_TidalDisruptions(seed, t, z, k_tdeBHstar, N_tdeBHstar, tde_type, m_avg, m_star, R_star, BH, v_star, vBH, tdes, binaries, pairs, f_accreted, EoS)

    # write back:
    state['seed'] = seed; state['t'] = t; state['z'] = z
    state['v_star'] = v_star; state['vBH'] = vBH
    state['binaries'] = binaries; state['pairs'] = pairs; state['tdes'] = tdes
    state['N_WD'] = N_WD; state['v_WD'] = v_WD
    state['N_tdeBHWD'] = N_tdeBHWD; state['N_tdeBHstar'] = N_tdeBHstar
    state['k_tdeBHWD'] = k_tdeBHWD; state['k_tdeBHstar'] = k_tdeBHstar
    state['dN_WDformdt'] = dN_WDformdt; state['dN_WDevdt'] = dN_WDevdt
    state['dN_tdeBHWDdt'] = dN_tdeBHWDdt; state['dN_tdeBHstardt'] = dN_tdeBHstardt


def record_evolution(state):
    """Append a row of 71 time-dependent quantities to the evolution array.

    Records the current cluster and BH subsystem state (masses, radii, densities,
    velocities, timescales, event counts, TDE statistics) as a single row in
    state['evolution']. Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
    """
    
    # unpack current state into local variables:
    seed = state['seed']; t = state['t']; z = state['z']; dt = state['dt']
    m_avg = state['m_avg']; Mcl = state['Mcl']; rh = state['rh']; M_gas = state['M_gas']
    R_gal = state['R_gal']; v_gal = state['v_gal']
    t_rlx = state['t_rlx']; tBH_rlx = state.get('tBH_rlx', 1e100)
    n_star = state['n_star']; N_BH = state['N_BH']
    mBH_avg = state['mBH_avg']; mBH_max = state.get('mBH_max', 0.0)
    rh_BH = state.get('rh_BH', 0.0); rc_BH = state.get('rc_BH', 0.0)
    S = state.get('S', 0.0); xi = state.get('xi', 0.0)
    psi = state.get('psi', 0.0); psi_BH = state.get('psi_BH', 0.0)
    t_3bb = state.get('t_3bb', 1e100); t_2cap = state.get('t_2cap', 1e100)
    k_3bb = state.get('k_3bb', 0); k_2cap = state.get('k_2cap', 0)
    N_me = state['N_me']; N_BBH = state['N_BBH']
    N_meRe = state['N_meRe']; N_meEj = state['N_meEj']
    v_star = state['v_star']; vBH = state['vBH']
    nh_BH = state.get('nh_BH', 0.0); nc_BH = state.get('nc_BH', 0.0)
    na_BH = state.get('na_BH', 0.0)
    N_3bb = state['N_3bb']; N_2cap = state['N_2cap']
    N_3cap = state['N_3cap']; N_BHej = state['N_BHej']
    N_BBHej = state['N_BBHej']; N_dis = state['N_dis']
    N_ex = state['N_ex']
    t_bb = state.get('t_bb', 1e100); N_bb = state['N_bb']
    N_meFi = state['N_meFi']; N_me2b = state['N_me2b']
    t_ex1 = state.get('t_ex1', 1e100); t_ex2 = state.get('t_ex2', 1e100)
    k_ex1 = state.get('k_ex1', 0); k_ex2 = state.get('k_ex2', 0)
    N_ex1 = state['N_ex1']; N_ex2 = state['N_ex2']
    N_BHstar = state['N_BHstar']
    t_pp = state.get('t_pp', 1e100); k_pp = state.get('k_pp', 0)
    N_pp = state['N_pp']
    N_Triples = state['N_Triples']; N_ZLK = state['N_ZLK']
    N_WD = state['N_WD']; v_WD = state.get('v_WD', 0.0)
    k_tdeBHWD = state.get('k_tdeBHWD', 0); N_tdeBHWD = state['N_tdeBHWD']
    dN_WDformdt = state.get('dN_WDformdt', 0.0); dN_WDevdt = state.get('dN_WDevdt', 0.0)
    dN_tdeBHWDdt = state.get('dN_tdeBHWDdt', 0.0)
    k_tdeBHstar = state.get('k_tdeBHstar', 0)
    dN_tdeBHstardt = state.get('dN_tdeBHstardt', 0.0)
    N_tdeBHstar = state['N_tdeBHstar']
    N_tdeBBHstar = state['N_tdeBBHstar']

    # append a row of 71 time-dependent quantities to the evolution array:
    state['evolution'].append([seed, t, z, dt, m_avg, Mcl, rh, R_gal, v_gal, t_rlx, tBH_rlx, n_star, N_BH, mBH_avg, mBH_max, rh_BH, rc_BH, S,
                                       xi, psi, psi_BH, t_3bb, t_2cap, k_3bb, k_2cap, N_me, N_BBH, N_meRe, N_meEj, v_star, vBH,
                                       nh_BH, nc_BH, na_BH, N_3bb, N_2cap, N_3cap, N_BHej, N_BBHej, N_dis, N_ex, t_bb, N_bb,
                                       N_meFi, N_me2b, t_ex1, t_ex2, k_ex1, k_ex2, N_ex1, N_ex2, N_BHstar, t_pp, k_pp, N_pp, v_esc_cl(v_star, vBH),
                                       2*vBH, N_Triples, N_ZLK, N_WD, v_WD, k_tdeBHWD, N_tdeBHWD, dN_WDformdt, dN_WDevdt, dN_tdeBHWDdt, k_tdeBHstar,
                                       dN_tdeBHstardt, N_tdeBHstar, N_tdeBBHstar, M_gas])

    # rows of populations.txt (when populations are configured, -POP); the number of rows
    # written after each evolution row is kept for the last few, as the last evolution
    # rows are cut at the end of the run (write_output):
    record = state.get('populations')
    if record is not None:
        for row in state['pops'].rows(t, state['m_avg0'], Mcl):
            record.append(row)
        record.keep = 3 * state['pops'].size
        marks = state['population_marks']
        marks.append((len(state['evolution']), len(record)))
        del marks[:-4]

def compute_external_params(state, config):
    """Compute external/environmental parameters for the cluster.

    Calculates the Jacobi (tidal) radius, the dimensionless escape rate xi_e,
    and the dynamical friction timescale t_df. These are needed before evolving
    TDEs (xi_e enters the WD evaporation rate) and before updating the cluster
    mass and radius. Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.
    """
    
    # unpack current state into local variables:
    Mcl = state['Mcl']; rh = state['rh']; R_gal = state['R_gal']
    v_gal = state['v_gal']

    if R_gal > 0:
        # Jacobi radius for an isothermal (flat rotation curve) galaxy:
        # r_J^3 = G M / (Omega^2 - Phi'') = G M / (2 Omega^2), with Omega = v_gal / R_gal
        rJ = (G_Newton * Mcl * R_gal**2 / 2 / v_gal**2)**(1/3)

        # dimensionless escape rate:
        xi_e = xi_e0 * np.exp(10 * rh / rJ)

        # dynamical friction timescale
        # (Gnedin et al. 2014, Eq. 8; lnL = 5.8, eccentricity factor f_eps = 0.5):
        t_df = 0.45e3 * (R_gal / 1e3)**2 * v_gal / (Mcl / 1e5) / 2
    else:
        # at the galactic center (R_gal = 0, -R 0 or -GC 1): no galactic tides (the limit
        # r_J -> infinity, escape rate of an isolated cluster) and no dynamical friction
        xi_e = xi_e0
        t_df = np.inf
    state['xi_e'] = xi_e
    state['t_df'] = t_df


def update_cluster(state, config):
    """Update cluster mass, half-mass radius, galactocentric radius, and time.

    Applies stellar-evolution mass loss, relaxation-driven mass loss, adiabatic
    and relaxation-driven expansion, galactocentric inspiral via dynamical
    friction, and advances the simulation time and redshift. Updates state
    in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.

    Returns:
        bool: True if the simulation should continue, False if the cluster has
            dissolved (Mcl < 0), reached the galaxy center (R_gal < 0), or
            the redshift has reached zero.
    """
    
    # unpack current state into local variables:
    t = state['t']; dt = state['dt']; z = state['z']
    Mcl = state['Mcl']; rh = state['rh']; R_gal = state['R_gal']
    m_avg = state['m_avg']; m_avg0 = state['m_avg0']
    t_rlx = state['t_rlx']; t_cc = state['t_cc']
    zCl_form = state['zCl_form']
    xi_e = state['xi_e']
    t_df = state['t_df']
    M_gas = state['M_gas']; t_ge = state['t_ge']

    pops = state['pops']
    if not pops.multi:
        # average mass evolution:
        if t>t_sev:
            m_avg = m_avg0 * (t / t_sev)**(-nu_sev)
        else:
            m_avg = m_avg0

        # stellar mass loss:
        if t>t_sev:
            dMcl_sev = - nu_sev * Mcl * dt / t
        else:
            dMcl_sev = 0.0

        # relaxation mass loss:
        dMcl_rlx = - xi_e * Mcl * dt / t_rlx
    else:
        # several stellar populations (populations.py): mean mass of all stars, escape
        # (relaxation) mass loss of the cluster shared among the populations, and the
        # stellar-evolution mass loss of each population at its own age:
        m_avg = Mcl / sum(pops.numbers(t, m_avg0))
        dMcl_rlx = - xi_e * Mcl * dt / t_rlx
        dMcl_sev = pops.evolve(t, dt, m_avg0, dMcl_rlx, t_rlx)

    # total mass loss:
    dMcl = dMcl_sev + dMcl_rlx

    # total dynamical mass (stars + gas) — sets potential-dependent responses:
    M_dyn = Mcl + M_gas

    # gas expelled this step (exact exponential):
    M_gas_new = M_gas * np.exp(-dt / t_ge)
    dM_gas = M_gas_new - M_gas

    # stellar-wind adiabatic expansion (total potential):
    drh_sev = - dMcl_sev * rh / M_dyn

    # gas-expulsion adiabatic expansion (Hills, total potential):
    drh_gas = - dM_gas * rh / M_dyn

    # expansion due to relaxation (balanced +2 term uses stellar mass):
    if t>t_cc:
        drh_rlx = zeta * rh * dt / t_rlx + 2 * dMcl_rlx * rh / Mcl
    else:
        drh_rlx = 0.0

    # total expansion:
    drh = drh_sev + drh_rlx + drh_gas

    # galactocentric radius step 
    # (Gnedin et al. 2014, Eq. 7: dR^2/dt = -R^2/t_df  =>  dR/dt = -R/(2 t_df)):
    dR_gal = - R_gal * dt / (2 * t_df) if R_gal > 0 else 0.0

    # cluster mass update:
    Mcl = Mcl + dMcl
    if pops.multi and Mcl >= 0:
        Mcl = sum(pops.M)      # (the same up to rounding; the populations are the reference)

    # gas expulsion (exact exponential drain on t_ge):
    M_gas = M_gas * np.exp(-dt / t_ge)

    if Mcl<0:
        print('CLUSTER DISSOLVED')
        state['Mcl'] = Mcl
        return False

    # half-mass radius update:
    rh = rh + drh

    if not rh > 0:
        # (possible only if the time step is dt_min and the escape rate is extreme)
        print('CLUSTER DISSOLVED')
        state['Mcl'] = Mcl
        return False

    # galactocentric radius update:
    R_gal = R_gal + dR_gal

    if R_gal<0:
        print('CLUSTER REACHED GALAXY CENTER')
        if not config.get('continue_at_GC', 0):
            state['R_gal'] = R_gal
            return False
        # -GC 1: the cluster stays at the center, without galactic tides from now on
        print('CONTINUING AT THE GALACTIC CENTER WITHOUT GALACTIC TIDES (t = {:.6g} Myr)'.format(t + dt))
        R_gal = 0.0

    # simulation time update:
    t+=dt

    # redshift update:
    try:
        z = redshift_interp(lookback_interp(zCl_form) - t)
    except Exception:
        print('REDSHIFT 0 REACHED')
        state['t'] = t
        return False

    state['N_iter'] = state['N_iter'] + 1
    state['t'] = t; state['z'] = z; state['dt'] = dt
    state['Mcl'] = Mcl; state['rh'] = rh; state['R_gal'] = R_gal
    state['M_gas'] = M_gas
    state['m_avg'] = m_avg

    return True


def print_status(state, config, local_time_initial, simulation_time_initial):
    """Print a summary of the current iteration to stdout.

    Displays cluster properties (t, dt, z, Mcl, rh, R_gal), BH subsystem
    counts (N_BH, N_BBH, N_Triples, N_me, N_tdeBHWD, N_tdeBHstar), and
    timing information (step duration, cumulative runtime). Only prints if
    config['print_info'] == 1.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.
        local_time_initial (float): Wall-clock time at the start of this iteration.
        simulation_time_initial (float): Wall-clock time at the start of the simulation.
    """
    
    # skip printing if disabled:
    if config['print_info']!=1:
        return

    t = state['t']; dt = state['dt']; z = state['z']
    Mcl = state['Mcl']; rh = state['rh']; R_gal = state['R_gal']
    N_BH = state['N_BH']; N_BBH = state['N_BBH']; N_Triples = state['N_Triples']
    N_me = state['N_me']; N_tdeBHWD = state['N_tdeBHWD']; N_tdeBHstar = state['N_tdeBHstar']; N_tdeBBHstar = state['N_tdeBBHstar']
    N_iter = state['N_iter']

    local_time_final = time.time()

    print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    print("Iteration #", N_iter)
    frmt_1 = '%.2f'
    frmt_2 = '%.0f'
    frmt_3 = '%.3f'
    frmt_4 = "%.1f"
    data_1 = {"t[Myr]": [frmt_1%t], "dt[Myr]": [frmt_1%dt], "z": [frmt_1%z], "Mcl[MMsun]": [frmt_1%(Mcl/1e6)], "rh[pc]": [frmt_1%rh], "R_gal[kpc]": [frmt_1%(R_gal/1e3)]}
    data_2 = {"N_BH": [frmt_2%N_BH], "N_BBH": [frmt_2%N_BBH], "N_Triples": [frmt_2%N_Triples], "N_me": [frmt_2%N_me], "N_tdeBHWD": [frmt_2%N_tdeBHWD], "N_tdeBHstar": [frmt_2%N_tdeBHstar], "N_tdeBBHstar": [frmt_2%N_tdeBBHstar]}
    data_3 = {"steptime[ms]": [frmt_4%(np.abs(local_time_final - local_time_initial)*1e3)], "runtime[s]": [frmt_3%np.abs(time.time() - simulation_time_initial)]}
    print(_one_row_table(data_1))
    print(_one_row_table(data_2))
    print(_one_row_table(data_3))
    print("\n")


def _one_row_table(data):
    """Text of print(pd.DataFrame(data, [" "])) for a one-row table of strings.

    Building three DataFrames per iteration cost ~3 ms, comparable to a whole
    simulation step in small clusters. The layout below is pandas' (right-aligned
    columns, a leading space in front of every value, a 1-character index); wider
    tables, which pandas would truncate, are still handed to pandas.
    """
    keys = list(data)
    vals = [data[k][0] for k in keys]
    widths = [max(len(k), len(v) + 1) for k, v in zip(keys, vals)]
    head = ' ' + ' '.join(k.rjust(w) for k, w in zip(keys, widths))
    if len(head) >= 80:
        return str(pd.DataFrame(data, [" "]))
    row = ' ' + ' '.join((' ' + v).rjust(w) for v, w in zip(vals, widths))
    return head + '\n' + row


# (_write_rows: arraycache.write_rows, which the background writer of StreamedRowLog uses too)
_write_rows = write_rows


def _finish_record(record, indicator, results_dir, name, keys, n_rows, prepare=None):
    """Write the output file of a record (evolution, hardening or TDEs) with its first
    n_rows data rows (all if None): a StreamedRowLog writes its remaining rows and
    renames its file; a list of rows or a RowLog (one placeholder row first) is
    written here."""
    if isinstance(record, StreamedRowLog):
        record.finish(n_rows)
        return
    if indicator != 1:
        return
    rows = np.array(record, dtype=float)[1:]   # drop placeholder row
    if n_rows is not None:
        rows = rows[:n_rows]
    if prepare is not None:
        rows = prepare(rows)
    with open(os.path.join(results_dir, name + '.txt'), 'w') as f:
        f.write('# ' + ' '.join(keys) + '\n')
        _write_rows(f, rows)


def write_output(state, config):
    """Export simulation results to the output directory.

    Writes the following files depending on config indicators:
      - outputBHs.pkl: BH masses, spins, and generations at every timestep.
      - tdes.txt: Tidal disruption event parameters.
      - mergers.txt: BBH merger source parameters (plus initial/final cluster state).
            - evolution.txt: Time-dependent cluster and BH subsystem quantities (71 columns).
      - hardening.txt: BBH hardening track details (12 columns).

    Args:
        state (dict): Simulation state at the end of the run.
        config (dict): Configuration dictionary.
    """
    CURRENT_WORKING_DIR = os.getcwd()
    RESULTS_DIR = os.path.join(CURRENT_WORKING_DIR, config['results_folder_name'])
    os.makedirs(RESULTS_DIR, exist_ok=True)

    if config['BOi']==1:
        data_to_save = {
            "t": state['simulation_times'],
            "mBH": state['black_hole_masses'],
            "sBH": state['black_hole_spins'],
            "gBH": state['black_hole_generations'],
            "hBH": state['black_hole_tdes']
        }
        BOF_path = os.path.join(RESULTS_DIR, config['BOF'] + '.pkl')
        with open(BOF_path, "wb") as f:
            dump_streaming(data_to_save, f)
        for key in ('black_hole_masses', 'black_hole_spins', 'black_hole_generations', 'black_hole_tdes'):
            if isinstance(state[key], SnapshotStream):
                state[key].close()

    N_tdeBHWD = state['N_tdeBHWD']
    N_tdeBHstar = state['N_tdeBHstar']
    N_tdeBBHstar = state['N_tdeBBHstar']
    N_me = state['N_me']
    N_iter = state['N_iter']
    N_hardening = state['N_hardening']
    Mcl0 = state['Mcl0']; rh0 = state['rh0']; Z = state['Z']
    zCl_form = state['zCl_form']; R_gal0 = state['R_gal0']
    Mcl = state['Mcl']; rh = state['rh']; R_gal = state['R_gal']

    # every stored TDE row; a whole-binary BBH-star disruption (type 22) adds one row
    # per BH, so there can be more rows than N_tdeBHWD+N_tdeBHstar+N_tdeBBHstar:
    _finish_record(state['tdes'], config['Ti'], RESULTS_DIR, config.get('tdes_file', 'tdes'),
                   tdes_keys, None, prepare=fill_R_mb)

    # initial and final cluster properties, the same on every merger row:
    cluster_columns = ' '+str(Mcl0)+' '+str(rh0)+' '+str(Z)+' '+str(zCl_form)+' '+str(R_gal0)+' '+str(Mcl)+' '+str(rh)+' '+str(R_gal)
    state['merger_suffix'] = cluster_columns     # (restart.json of a restarted run)
    if isinstance(state['mergers'], StreamedRowLog):
        state['mergers'].finish(max(N_me, 0), suffix=cluster_columns)
    elif config['Mi']==1:
        mergers = np.array(state['mergers'], dtype=float)[1:]   # drop placeholder row
        with open(os.path.join(RESULTS_DIR, config['mergers_file'] + '.txt'), 'w') as f_mergers:
            f_mergers.write('# ' + ' '.join(merger_keys) + '\n')
            _write_rows(f_mergers, mergers[:max(N_me, 0)], cluster_columns)

    _finish_record(state['evolution'], config['Ei'], RESULTS_DIR, config.get('evolution_file', 'evolution'),
                   evolution_keys, max(N_iter-1, 0))

    _finish_record(state['hardening'], config['Hi'], RESULTS_DIR, config.get('hardening_file', 'hardening'),
                   hardening_keys, max(N_hardening, 0))

    # populations.txt: the rows of the time steps kept in evolution.txt (its first
    # N_iter - 1 rows):
    if state.get('populations') is not None:
        n_evo = max(N_iter - 1, 0)
        n_pop = 0
        for n_e, n_p in state['population_marks']:
            if n_e - 1 == n_evo:            # (n_e counts the placeholder row)
                n_pop = n_p - 1
        if n_evo > 0 and n_pop == 0 and state['population_marks']:
            raise RuntimeError('populations.txt: rows of the last time steps not found')
        state['populations'].finish(n_pop)

# End of file.
