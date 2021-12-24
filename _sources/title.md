# Displacement method

:::{note}
I tried to turn a published paper "Efficacious symmetry-adapted atomic displacement method for lattice
dynamical studies", Computer Physics Communications 259 (2021) 107635 into a JupyterBook.
:::

## Abstract

Small displacement methods have been successfully used to calculate the
lattice dynamical properties of crystals.  It involves displacing atoms
by a small amount in order to calculate the induced forces on all atoms
in a supercell for the computation of force constants.  Even though
these methods are widely in use, to our knowledge, there
is no systematic discussion of optimal displacement directions from
the crystal's symmetry point of view nor a rigorous error analysis of
such methods.  Based on the group theory and point group symmetry of a
crystal, we propose displacement directions, with an equivalent concept
of the group of $k$, deduced directly in the
Cartesian coordinates rather than the usual fractional coordinates,
that maintain the theoretical maximum for the triple product $V$
spanned by the three displacements to avoid possible severe roundoff
errors. The proposed
displacement directions are generated from a minimal set
of irreducible atomic displacements
that keep the required independent force calculations to a minimum.
We find the error in the calculated force constants explicitly depends on
the inverse of $V$ and inaccuracy of the forces. Test systems such as Si,
graphene, and orthorhombic Sb2S3 are used to illustrate the
method.  Our symmetry-adapted atomic displacement method is shown to be very robust in treating
low-symmetry cells with a large `aspect ratio' due to huge differences
in lattice parameters, use of a large vacuum height,
or a very oblique unit cell due to unconventional choice of primitive lattice vectors.
It is expected that our atomic displacement
strategy can be used to address higher-order interatomic interactions
to achieve good accuracy and efficiency.
