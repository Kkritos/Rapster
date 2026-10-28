"""A ZLK merger of a hierarchical triple that would happen after redshift 0 has not
happened by the end of the simulation: the triple is kept, nothing is counted."""
import numpy as np

from rapster.triples import evolve_triples
from rapster.functions import lookback_interp, set_escape_factor
from rapster.bhpool import BHPool

ZCL = 3.0
set_escape_factor(2.0)


def _setup():
    BH = BHPool(np.array([8.0, 12.0, 20.0]), np.zeros(3), np.ones(3), np.zeros(3))
    triples = np.zeros((2, 24))
    #               a_in   a_out  e_in e_out m0    m1    m2    s0 s1 s2 g0 g1 g2 incl1 incl2 t  z  ind ch tform zform h0 h1 h2
    triples[1] = [5e-7, 2e-4, 0.1, 0.3, 30.0, 25.0, 15.0, 0, 0, 0, 1, 1, 1, 0.8, 0.75, 0, 0, 7, 3, 0, ZCL, 0, 0, 0]
    # (t_ZLK ~ 1.3 Gyr: shorter than the age of the universe; the next interaction is far away)
    binaries = np.zeros((1, 15))
    return BH, triples, binaries


def _run(t, nc_BH):
    BH, triples, binaries = _setup()
    out = evolve_triples(1, t, 0.0, ZCL, triples, binaries, BH, 10.0, 1, 0, 6, 0, 0, 0, 0,
                         10.0, 10.0, nc_BH, [])
    (_, _, _, _, triples, binaries, BH, _, N_Triples, N_BBH, N_BH, N_me, N_meRe, N_meEj, N_ZLK,
     _, _, _, mergers) = out
    return triples, binaries, BH, N_Triples, N_BBH, N_BH, N_me, N_meRe, N_meEj, N_ZLK, mergers


def test_merger_after_today_keeps_the_triple():
    np.random.seed(3)
    t_today = lookback_interp(ZCL)
    triples, binaries, BH, N_Triples, N_BBH, N_BH, N_me, N_meRe, N_meEj, N_ZLK, mergers = \
        _run(t_today - 1e-9, nc_BH=1e-30)        # the next interaction is far away
    assert N_Triples == 1 and triples.shape[0] == 2           # triple kept
    assert N_me == N_ZLK == N_meRe == N_meEj == 0 and mergers == []
    assert N_BBH == 0 and binaries.shape[0] == 1 and BH.size == 3 and N_BH == 6
    assert BH.size + 2 * N_BBH + 3 * N_Triples == N_BH       # every BH still held


def test_merger_before_today_is_resolved_consistently():
    np.random.seed(3)
    triples, binaries, BH, N_Triples, N_BBH, N_BH, N_me, N_meRe, N_meEj, N_ZLK, mergers = \
        _run(0.0, nc_BH=1e-30)
    assert N_Triples == 0 and triples.shape[0] == 1             # merged
    assert N_me == N_ZLK == 1 == len(mergers) and N_meRe + N_meEj == 1
    assert mergers[0][16] < lookback_interp(ZCL)                # merger time before today
    # N_BH = 6 - 1 (merger) - 1 if the remnant was ejected, and equals the BHs held
    # (singles + 2 per binary; no triple left):
    assert N_BH == 5 - N_meEj and BH.size + 2 * N_BBH == N_BH
