'''
Cosmological coupling of the black holes: M ~ a^k (k = -3 ... 3), applied once per
time step after the BH formation time, to the single BHs (BHPool), the BBHs and the
BH-star pairs.
'''

import numpy as np


def apply_cosmological_coupling(t, tBH_form, z, z_i, z_i_0, k_cosmo, binaries, BH, pairs, m_avg, triples=None):
    """
    Scale masses, spins and semi-major axes by factor = ((1+z_i)/(1+z))**k_cosmo
    (z_i is the redshift of the previous application). Modifies BH, binaries, pairs
    and triples in place. Returns the redshift z_i to use at the next step.
    """
    if t > tBH_form and z <= z_i_0:
        factor = ((1 + z_i) / (1 + z))**k_cosmo

        if binaries.shape[0] > 1:
            binaries[:, 4] *= factor
            binaries[:, 5] *= factor
            # exact BBH semi-major axis scaling (both components grow equally):
            binaries[:, 2] *= factor**(-3)
            binaries[:, 6] = np.minimum(1.0, binaries[:, 6] * (factor**(-2)))  # spin: χ=cJ/(GM²), J conserved
            binaries[:, 7] = np.minimum(1.0, binaries[:, 7] * (factor**(-2)))  # spin: χ=cJ/(GM²), J conserved

        if BH.size > 0:
            BH.set_all_masses(BH.m * factor)
            BH.set_all_spins(np.minimum(1.0, BH.s * (factor**(-2)))  # spin: χ=cJ/(GM²), J conserved)

        if pairs.shape[0] > 1:
            # BH mass fraction of the pair, using CURRENT masses
            # (must be computed before pairs[:, 1] is updated below):
            m_BH_current = pairs[:, 1]
            x = m_BH_current / (m_BH_current + m_avg)

            # exact BH-star semi-major axis scaling (only the BH mass grows):
            pairs[:, 0] *= ((1 - x) + x * factor) / factor**2
            pairs[:, 1] *= factor
            pairs[:, 2] = np.minimum(1.0, pairs[:, 2] * (factor**(-2)))  # spin: χ=cJ/(GM²), J conserved

        if triples is not None and triples.shape[0] > 1:
            # all three BHs grow by f: the inner and the outer orbit both conserve their
            # angular momentum, a -> a f^-3 (mass ratios, hence the hierarchy and
            # stability of the triple, and the eccentricities are unchanged):
            triples[:, 0] *= factor**(-3)
            triples[:, 1] *= factor**(-3)
            triples[:, 4:7] *= factor
            triples[:, 7:10] = np.minimum(1.0, triples[:, 7:10] * (factor**(-2)))  # spin: χ=cJ/(GM²), J conserved

        return z

    return z_i
