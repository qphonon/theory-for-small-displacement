# Introduction


Lattice dynamical studies\cite{Born56-book} are important for an 
understanding of the 
phase stability,\cite{VanDeWalle02v74,Gan10v49} ferroelectric transition,\cite{Zhong94v72} Raman and 
infra-red spectroscopies,\cite{Zhao13v13,Chong14v90,Giovanni18v5} 
phonon-mediated superconductivity,\cite{Giustino17v89}
ferroelastic transition,\cite{Togo08v78} 
and
thermodynamics of materials.\cite{Grimvall1999-book,Mujica03v75}
Even though the concept of phonon emerges as a result of a harmonic approximation, its simple extensions
via the quasiharmonic approximation (QHA)\cite{Mounet05v71}
and Gr\"uneisen formalism\cite{Allen20v34,Malica20v32} 
allow thermal properties due to anharmonic effects such as 
lattice thermal conductivities\cite{Toher14v90}, a key quantity that determine the figure of merits for thermoelectrics\cite{Snyder08v7,Madsen16v213,Gorai17v2}, as well as 
thermal expansion  coefficients\cite{Gan15v92,Arnaud16v93,Gan16v94,Romao17v96,Gan18v151,Gan19v31} to be evaluated.  An attempt has been made to extract  the third order interatomic force constants from
the standard phonon calculations through the evaluation of Gr\"uneisen parameters\cite{Lee17v96}.
Phonon databases are being built\cite{Petretto18v5,TogoPhononDB2020-link} with the
help of high-throughput frameworks\cite{Curtarolo13v12,Ong13v68,Pizzi16v111} for data mining and machine learning.

The methods to calculate phonon frequencies
and eigenvectors naturally fall into two distinct approaches. 
One approach is based on the displacement methods\cite{Parlinski97v78,Wang14v185,Wang16v2,Ackland97v9,Togo15v108,Togo2020-github,Alfe09v180,Kresse95v32,Gan10v49} that gain popularity due to their
simplicity in implementations. In these methods,
 the force constants can be deduced from the induced forces via
the Hellman-Feynman theorem\cite{Feynman39v56} when 
a small displacement on an atom in an otherwise perfect supercell is made.
The possibility to calculate exact 
phonon frequencies at commensurate $\VEC{q}$ vectors
using smaller non-diagonal unit cells\cite{Lloyd-Williams15v92} paves the way for
very practical applications of these methods.
The second 
approach is
based on the density-functional perturbation theory.\cite{Baroni01v73,Gonze97v55a}
In this approach, the methods are very versatile because
small unit cells rather than huge supercells
are required. Very accurate energy derivatives could be analytically calculated
within a computer code\cite{Giannozzi09v21,Gonze09v180}.
Both approaches could be used together to 
complement each other\cite{Fu19v100} for the extraction of higher-order 
interatomic interactions.\cite{Parlinski18v98}

A somewhat simple and effective implementation of
a small displacement method is to use
a few pre-selected directions in the fractional coordinates that are to be acted upon by
appropriate space group operations to deduce the displacement
directions as long as the volume $V$
spanned by the actual displacement directions is nonzero. 
This strategy may have been inspired by the fact that
space group operations\cite{ITtable06-book} usually act on the 
positions of atoms in fractional coordinates within a unit cell.
However, as we shall show later, such
implementation may result in a nonoptimal choice of the displacements of atoms that 
affects the accuracy of the force constants and eventually the lattice dynamical
properties. We use a generic
 orthorhombic system of space group $Pnma$ with lattice parameters 
$a$, $b$, and $c$
to show that the volume $V$ may deviate from its ideal value of 1 and 
scale unfavorably
as $2(a/b)^{-1} \rightarrow 0 $ as  $a \gg b$. 
In the case of simulating a graphene sheet using a supercell method, it is 
shown that a large vacuum thickness could reduce the ideal $V=1$ to 
a value that scales as $(c/a)^{-2}$ where $c$ is the vacuum height and $a$ the 
hexagonal lattice parameter for graphene. In this paper
we propose a displacement method that can be applied to any crystal, encompassing the 
entire 32 crystallographic point groups and 230 space groups.
The displacement directions are deduced directly in the Cartesian
coordinates rather than fractional coordinates
that are designed to maintain (i) the theoretical maximum for the 
triple product $V$ spanned by the three displacements
to avoid possible severe roundoff errors, and (ii) a minimal set of irreducible atomic displacements for independent force calculations.
To achieve these aims we rely on the concepts of the star of $k$ and the group of $k$\cite{Dresselhaus2008-book}, defined originally in the 
reciprocal space but extended to real space in this paper.
Various test systems such as Si, graphene, and orthorhombic
\sbs{} are used to illustrate the method.
This paper is organized as follows. In Section~\ref{sec:method}
we provide the full details 
of our displacement method to make judicious atomic displacements for all crystal symmetries.   Section~\ref{sec:FC_error} 
presents an error analysis for the force constants. Results are shown in Section~\ref{sec:results}. Finally we conclude in Section~\ref{sec:conclusion}. Appendix
A illustrates more clearly how $V$ may scale poorly with the lattice parameters or the choice of
oblique unit cell for a few selected cases.

# Methodology

First we define a matrix

$$
A=[ \vec{a}_1|\vec{a}_2 | \vec{a}_3]
$$ (eq:A)
where
the $i$th column of $A$ is taken from
the basic lattice translation vector $\vec{a}_i$ of a crystal.
For simplicity
we use $A$ to describe a primitive cell but it could be easily
extended to deal with a conventional unit cell or even a nonconventional
unit cell.
A space group operator $\{ R | \vec{t} \}$
(in the Seitz notation) corresponds
to a rotation matrix

$$
R_c = A^{-1} R A
$$ (eq:Rc)
in the Cartesian coordinates, which is restricted to
either $p$ or ${\overline p}$ for $p= 1$, $2$, $3$, $4$, and $6$.
$p$ means a $p$-fold rotation in the international
notation, while $ {\overline p}$ means an
improper rotation (i.e., a $p$-fold rotation followed by an inversion).
The 230 space groups are built on top of
the 32 crystallographic
point groups (see e.g., Refs.~[{cite}`Dresselhaus2008-book`.)



Example

$$
\begin{align}
a_{11}& =b_{11}&
  a_{12}& =b_{12}\\
a_{21}& =b_{21}&
  a_{22}& =b_{22}+c_{22}
\end{align}
$$


...

```{bibliography}
```
