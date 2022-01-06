<!---
Example of a block comment.
\begin{figure}
\centering\includegraphics[width=9.2cm,clip]{fig/2020-04-22-mpl-qphonon-det-gives-error-fit.pdf}
\end{figure}
--->

# Introduction

Lattice dynamical studies{cite}`Born1956-book` are important for an 
understanding of the 
phase stability,{cite}`VanDeWalle02v74,Gan10v49` ferroelectric transition,{cite}`Zhong94v72` Raman and 
infra-red spectroscopies,{cite}`Zhao13v13,Chong14v90,Giovanni18v5` 
phonon-mediated superconductivity,{cite}`Giustino17v89`
ferroelastic transition,{cite}`Togo08v78` 
and
thermodynamics of materials.{cite}`Grimvall1999-book,Mujica03v75`
Even though the concept of phonon emerges as a result of a harmonic approximation, its simple extensions
via the quasiharmonic approximation (QHA){cite}`Mounet05v71`
and Gr\"uneisen formalism{cite}`Allen20v34,Malica20v32` 
allow thermal properties due to anharmonic effects such as 
lattice thermal conductivities{cite}`Toher14v90`, a key quantity that determine the figure of merits 
for thermoelectrics{cite}`Snyder08v7,Madsen16v213,Gorai17v2`, as well as 
thermal expansion  coefficients{cite}`Gan15v92,Arnaud16v93,Gan16v94,Romao17v96,Gan18v151,Gan19v31` to be evaluated.  An attempt has been made to extract  the third order interatomic force constants from
the standard phonon calculations through the evaluation of Gr\"uneisen parameters{cite}`Lee17v96`.
Phonon databases are being built{cite}`Petretto18v5,TogoPhononDB2020-link` with the
help of high-throughput frameworks{cite}`Curtarolo13v12,Ong13v68,Pizzi16v111` for data mining and machine learning.

The methods to calculate phonon frequencies
and eigenvectors naturally fall into two distinct approaches. 
One approach is based on the displacement methods{cite}`Parlinski97v78,Wang14v185,Wang16v2,Ackland97v9,Togo15v108,Togo2020-github,Alfe09v180,Kresse95v32,Gan10v49` that gain popularity due to their
simplicity in implementations. In these methods,
 the force constants can be deduced from the induced forces via
the Hellman-Feynman theorem{cite}`Feynman39v56` when 
a small displacement on an atom in an otherwise perfect supercell is made.
The possibility to calculate exact 
phonon frequencies at commensurate $\vec{q}$ vectors
using smaller non-diagonal unit cells{cite}`Lloyd-Williams15v92` paves the way for
very practical applications of these methods.
The second 
approach is
based on the density-functional perturbation theory.{cite}`Baroni01v73,Gonze97v55a`
In this approach, the methods are very versatile because
small unit cells rather than huge supercells
are required. Very accurate energy derivatives could be analytically calculated
within a computer code{cite}`Giannozzi09v21,Gonze09v180`.
Both approaches could be used together to 
complement each other{cite}`Fu19v100` for the extraction of higher-order 
interatomic interactions.{cite}`Parlinski18v98`

A somewhat simple and effective implementation of
a small displacement method is to use
a few pre-selected directions in the fractional coordinates that are to be acted upon by
appropriate space group operations to deduce the displacement
directions as long as the volume $V$
spanned by the actual displacement directions is nonzero. 
This strategy may have been inspired by the fact that
space group operations{cite}`ITtable2006-book` usually act on the 
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
To achieve these aims we rely on the concepts of the star of $k$ and the group of $k${cite}`Dresselhaus2008-book`, defined originally in the 
reciprocal space but extended to real space in this paper.
Various test systems such as Si, graphene, and orthorhombic
\sbs{} are used to illustrate the method.
This paper is organized as follows. In Section 2
we provide the full details 
of our displacement method to make judicious atomic displacements for all crystal symmetries.   Section 3
presents an error analysis for the force constants. Results are shown in Section 4. Finally we conclude in Section 5. The Appendix
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
point groups (see e.g., Refs.~[{cite}`Dresselhaus2008-book,Burns1985-book`]).

In the small displacement method{cite}`Liu14v16` for
a phonon calculation of a crystal,
if the intrinsic symmetries of a 
space group are not
utilized, we have to sequentially
displace all $N_1$ atoms in a primitive cell
embedded in a large $n_1 \times n_2 \times n_3$ supercell of
$n_1 n_2 n_3 N_1$ atoms
along the $x$, $y$, and $z$ Cartesian axes in 
both positive and negative directions,
resulting in
$6N_1$ different 
supercells, 
each
with one atom displaced slightly compared to the unperturbed supercell (we call this the all-displacement method).
Each of the $6N_1$ supercells will be treated independently
where the induced forces are to be calculated so that
the interatomic force constants can be deduced.
Due to a computational cubic scaling with respect to the number of atoms, most density-functional theory (DFT)
implementations{cite}`Payne92v64,Gan01v134` face a severe practical
issue to handle $6N_1$ supercells to evaluate the induced forces. It is therefore important
to reduce the number of calculations as much as possible.

We now review how the force constant matrix $\Psi_{ij}$
 between
the $i$th atom in the primitive cell and the $j$th atom in the
supercell may be calculated.
This is achieved by sequentially displacing the $i$th
atom in the primitive cell by three displacement vectors $\lambda \vec{d}_k^i$ in the 
Cartesian coordinates, $k = 1, 2, 3$, where
$\vec{d}_k^i$ is a unit vector. 
$\lambda$ is the magnitude of the displacement which is typically
$0.010 \sim 0.015$~\AA.
For each displacement vector $\lambda \vec{d}_k^i$,
we calculate the induced forces
on all atoms in the supercell, in particular the induced force
$\vec{F}_k^j$ on the $j$th atom, $k=1,2,3$.
If we make $\vec{F}_k^j$ to 
constitute the $k$th column of a matrix $F^j$ where

$$
F^j = [\vec{F}_1^j | \vec{F}_2^j | \vec{F}_3^j ]
$$ (eq:Fj)
and similarly we
let $\vec{d}_k^i$ to constitute the $k$th column of a displacement matrix $d^i$ where

$$
d^i = [\vec{d}_1^i | \vec{d}_2^i | \vec{d}_3^i] 
$$ (eq:dmat)
then the required force constant $\Phi_{ij}$ can be calculated from

$$
F^j = \lambda \Phi_{ij} d^i
$$ (eq:Phi)
From Eq.{eq}`eq:Phi`
we see that the displacement vectors $\vec{d}_k^i$ do not
need to point along the conventional Cartesian axes but in
any directions as long as $d^i$ is not singular.

A first level of a possible reduction of the number of calculations can be made
if an inequivalent atom in the primitive cell, 
where its position is commonly known 
as the Wyckoff position,
could be mapped under space group operations
to its equivalent atoms in the primitive cell, resulting in a 
star of $k$, which is understood to be a set of atoms.
We note that the concept of the star of $k$ is originally defined in 
the reciprocal space{cite}`Dresselhaus2008-book`.
When handling a star of $k$, we may arbitrarily choose
any atom in a set to be a representative atom (the so-called
inequivalent atom) for the 
purpose of generating the rest of equivalent atoms in the same set. 
The $3\times 3$ force constant matrix $\Phi_{ij}$
between the $i$th inequivalent atom in the
primitive cell and another $j$th atom in the supercell
can be ``copied'' out{cite}`Liu14v16` to 
another $3\times 3$ force constant matrix $\Phi_{i'j'}$ between 
the $i'$th equivalent atom and the $j'$th atom using

$$
\Phi_{i'j'} =  R_c \Phi_{ij} R_c^T
$$ (eq:Phiij)

$R_c$ is as defined in Eq.{eq}`eq:Rc` and $ \{ R | \vec{t} \}  $
maps the $i$th atom to the $i'$th atom, and
the $j$th atom to the $j'$th atom.  $R_c^T$ is the matrix transpose of $R_c$.
With this strategy alone (we call this the 6-displacement method), 
if there are $N_0$ inequivalent atoms in a primitive cell, 
then only $6 N_0$ calculations are needed. 

However, it may still be possible to reduce the 6 calculations for each
of $N_0$ inequivalent atoms using the 
site symmetry{cite}`Burns1985-book` of an inequivalent atom. The
site symmetry (this has the equivalent concept of the group of $k$, see
for example Ref.~[\onlinecite{Altmann-book}]) of an inequivalent
atom is one of the 32 
crystallographic point groups that
leaves the position of an inequivalent atom invariant in periodic
sense.  Note that two different inequivalent atoms from two different 
stars of $k$ may not have the same site symmetries.

Now we shall discuss how the site symmetry can be used
to reduce the number of displacements for an inequivalent atom.
Suppose a displacement $\vec{d}_1^i$ has been applied to an inequivalent
$i$th atom and the forces on all atoms in the supercell have been found.
If an element $\{R|\vec{t}\}$ in the site symmetry of the $i$th
atom is applied to the supercell with the $i$th atom that has been displaced by
$\vec{d}_1^i$, then the net effect of the operation
is to rotate the original displacement
$\vec{d}_1^i$ to become a displacement $\vec{d}_2^i = R_c \vec{d}_1^i$
on the $i$th atom.  The operation $\{R|\vec{t}\}$ also reshuffles
the positions of all atoms in the supercell, as well as to rotate the
induced forces caused by $\vec{d}_1^i$ on all atoms in the supercell. This
crucial observation implies 
that without doing an independent (probably expensive)
induced force calculation
due to $\vec{d}_2^i$, we are able to just use the information due to
$\vec{d}_1^i$ to give us all force information for
a virtual $\vec{d}_2^i$ displacement.  
If there is yet another operation $\{R|\vec{t}\}$ that could
rotate $\vec{d}_1^i$ to $\vec{d}_3^i$, then again an independent displacement
of $\vec{d}_3^i$ does not need to be carried out.  However, in the case when
there is no extra operation in a site symmetry that can
generate an independent $\vec{d}_2^i$ or $\vec{d}_3^i$,
then we have no choice
but to carry out necessary separate displacements and find the induced
forces to fill up the necessary force field for the $\Phi_{ij}$ calculation.


In an elegant implementation, the displacement in the Cartesian
coordinates $\vec{d}^i_k$ can be generated from the directions defined
in the fractional coordinates $\vec{g}_k^i$ where 
$ \vec{d}^i_k =  \frac{A \vec{g}_k^i}{ |  A \vec{g}_k^i  |   }$.  
$\vec{g}_{k}^i$ may be chosen from a set of $S$ 
that consists of nonzero
vectors of the form $ (e_1, e_2, e_3)^T$, where $e_n = 0, \pm
1$ for $n= 1,2,3$ for simplicity.  By systematically applying the elements
of the site symmetry to vectors in $S$, one may find a minimal set
$S_i$ that contains between one to three vectors. 
$S_i$ will generate three $\vec{g}_k^i$ that form a nonzero determinant for $d^i$.
This implementation has a slight drawback that may be illustrated by
two similar
orthorhombic systems Bi$_2$S$_3${cite}`Zhao11v84` and 
\sbs{cite}`Liu14v16,Chong14v90,Gan15v92`
in the $Pnma$ setting, where $a \sim 11.3$~\AA, $b \sim 3.8 $~\AA, and $c \sim
11.1$~\AA. Here $a$ is about three times larger than $b$.  
There are twenty atoms in the primitive cell, with
five inequivalent atoms on the $4c$ Wyckoff site. 
The site symmetry is a group of mirror reflection (a two-element group).
The reflection operator $\overline{2}$
reflects the system
across the $xz$ plane and maps a direction in 
the fractional coordinate 
$ \vec{g}_1^i = (1,-1,0)^T $ to $ \vec{g}_2^i = (1,1,0)^T$.
However, these two directions in the fractional coordinates correspond to
$(a,-b,0)^T$ and $(a,b,0)^T$ in the Cartesian coordinates which
subtend an angle not equal to an ideal
angle of $90^\circ$ since $a \neq b$.
The angle between the two displacements in
the Cartesian coordinates is actually
given by $2\tan^{-1}{\frac{b}{a}}$ (see Appendix A for a detailed discussion).
Hence if $a$ is much larger than $b$, the two vectors $\vec{d}_1^i$
and $\vec{d}_2^i$ are nearly parallel to each other in the Cartesian coordinates.
However the operation is a physical reflection in the $xz$ plane and
hence it is possible to force the angle between the two displaced vectors
to be exactly $90^\circ$, thereby achieving the largest 
determinant for $d^i$ of $1$ for matrix inversion in Eq.{eq}`eq:Phi`.

In another example, a graphene sheet of
a lattice constant $a$ and a vacuum thickness $c$ may use $\vec{g}_1^i =
(1,0,1)^T$ to generate $\vec{g}_2^i = (0, 1, -1)^T$ and $\vec{g}_3^i = 
(-1,-1,1)^T$ with the point group operations. 
This means only one independent force calculation is to be performed.
However, an analysis shows that $V = \det d^i = \frac{ a^2 c \sqrt3  }{ 2(a^2 + c^2)^{3/2}   }$, which 
goes to $\frac{\sqrt3/2}{(c/a)^2}$
for a large $c/a$ ratio. With{cite}`Gan10v81` $c=20$~\AA{} and $a = 2.471$~\AA{}, 
$V = 0.013$. However, with the method to be developed later,
we find that the angles between any two displacements taken from $\vec{d}_k^i$, $k=1,2,3$
can be made $90^\circ$ thereby making $V$ achieves its largest 
value of 1.

To develop a sense of how the inaccuracy of forces and
 $V = \det d^i$ may affect the
accuracy of force constants, we consider the force constants in the
$xy$ plane for the $4c$ site of the
$Pnma$ space group. Here $\vec{d}_1^i = \frac{1}{\sqrt{a^2 + b^2}}(a,-b)^T$ 
and $\vec{d}_2^i = \frac{1}{\sqrt{a^2+ b^2}}(a,b)^T$.
From Eq.{eq}`eq:Phi`, we have


$$
[\vec{F}_1^j | \vec{F}_2^j ] = \lambda [\Phi_{ij}^1| \Phi_{ij}^2] 
\frac{1}{\sqrt{a^2 + b^2}}
  \begin{pmatrix}
       a &  a \\
       -b &  b \\
  \end{pmatrix}
$$ (eq:Fj)

We then have

$$
\Phi_{ij}^2 = \frac{1}{2\lambda} \left( \frac{a^2}{b^2}  +1 \right)^{\frac{1}{2}}
 (\vec{F}_2^j -  \vec{F}_1^j)
$$ (eq:Phiij2)
If $a \gg b$, then
$\vec{d}_1^i$ and $\vec{d}_2^i$ are almost parallel to each other resulting in
very similar forces $\vec{F}_1^j$ and $\vec{F}_2^j$.
If the force $\vec{F}_1^j$ due to $\vec{d}_1^i$ is not 
determined accurate enough, $\vec{F}_2^j$ will inherit the same inaccuracy since it 
is `copied' from $\vec{F}_1^j$ through a symmetry operation, then
$\Phi_{ij}^2$ will be inaccurate due to a large roundoff error that is amplified by 
a large geometry factor $ \left( \frac{a^2}{b^2}  +1 \right)^{\frac{1}{2}}$.

Now we present a method
that will systematically deduce the displacement directions 
directly in the Cartesian coordinates for a
forward difference scheme with the aim of 
maintaining a largest possible magnitude for the determinant of $d^i$ and
a minimum number of independent displacements.
The displacement method
must also maintain the minimal number of independent displacements when a 
central difference scheme is used for an improved accuracy for $\Phi_{ij}$.
For a central difference scheme, we need to displace the atoms
in $-\vec{d}_k^i$, $k = 1, 2, 3$ where $\vec{d}_k^i$ have been chosen 
for a forward difference scheme. We note that two operations
are able to map $\vec{d}_k^i$ to $-\vec{d}_k^i$: one is the inversion
operator, and the other a 
2-fold rotation.


|Crystal | Symmetry  | $n_{\rm FD}$ | $n_{\rm CD}$ | 
|:-------|------------------|----------------|-----------------|
|Triclinic |    $C_1$  | 3 | 6 
| |    $C_i$ | 3 | 3 
|Monoclinic |    $C_2$ |  2 | 3 
| |    $C_{s}$ | 2  | 4 
| |    $C_{2h}$ | 2 | 2 
|Orthorhombic |    $D_2$ | 1 | 2 | 
| |    $C_{2v}$ | 1 | 2 | 
| |    $D_{2h}$ | 1 | 1 | 
|Tetragonal |    $C_4$ | 1 | 2 | 
| |    $S_4$ | 1 | 2 | 
| |    $C_{4h}$ | 1 | 1 | 
| |    $D_4$ | 1 | 1 | 
| |    $C_{4v}$ | 1 | 2 | 
| |    $D_{2d}$ | 1 | 1 | 
| |    $D_{4h}$ | 1 | 1 | 
|Trigonal |    $C_3$ |1 | 2 
| |    $S_6$ | 1 | 1 | 
| |    $D_3$ | 1 | 1 | 
| |    $C_{3v}$ | 1 | 2 | 
| |    $D_{3d}$ | 1 | 1 | 
|Hexagonal  |    $C_{6}$ | 1 | 2 | 
| |    $C_{3h}$ | 1 | 2 | 
| |    $C_{6h}$ | 1 | 1 | 
| |    $D_{6}$ | 1 | 1 | 
| |    $C_{6v}$ | 1 | 2 | 
| |    $D_{3h}$ | 1 | 1 | 
| |    $D_{6h}$ | 1 | 1 | 
|Cubic |    $T$ | 1 | 1 | 
| |    $T_h$ | 1 | 1 | 
| |    $O$ | 1 | 1 | 
| |    $T_d$ | 1 | 1 | 
| |    $O_h$ | 1 | 1 | 

Table 1: $n_{\rm CD}$ ($n_{\rm FD}$) is the minimal number of displacements per atom for a central (forward) difference scheme.


The task may at first seem arduous for all 32 crystallographic 
point groups (see Table 1). However, there are actually 
four distinct cases to consider. (See 
an implementation of our algorithm in fm-forces.f90
  from https://github.com/qphonon/atomic-displacement)
The first case deals with
the triclinic crystals and covers
point groups 1 to 2 ($C_1$ and $C_i$)
that involve a single $1$ or ${\overline 1}$ operation.  
The second case handles the monoclinic crystals and covers
point groups 3 to 5 ($C_2$, $C_s$, $C_{2h}$) 
that involve a single $2$  or ${\overline 2}$ operation. 
The third case deals with the orthorhombic 
crystals and covers point groups 6 to 8  ($D_2$, $C_{2v}$, $D_{2h}$).
These point groups possess
three operations involving 
either $2 $ or ${\overline 2}$ that are perpendicular
to one another.
Finally the fourth case
deals with the remaining four crystal systems, i.e., 
tetragonal, trigonal, hexagonal, and cubic crystals and covers
point groups 9 to 32.
These point groups have a single $3$, ${\overline 3}$, $4$, 
or ${\overline 4}$ operation.

Usually, the principal axis (the axis with the highest symmetry) 
is not pointing along the $z$ axis but
along a vector $\vec{n}$ characterized by 
the polar angle $\theta$ and azimuthal angle $\phi$.
To ease the discussion of the determination
of $\vec{d}_k^i$, we 
transform all operations $R_c$ to $R_c'$  in the site symmetry according to
$R_c' = Q^{-1} R_c Q$, $Q = R_z(\phi) R_y(\theta)$. Here
$R_y (\theta) $ [$R_z (\phi)$]
is a rotation of angle $\theta$ ($\phi$) around the $y$ ($z$) axis.
Under such transformation, the operation with the 
principal axis $R_{\vec{n}}(\eta)$  will be mapped to
$R_z(\eta)$ where the rotation axis is now along the $z$ axis. 
This can be easily seen since $R_{\vec{n}}(\eta)$ is transformed to
$ R_z(\eta) $ via $R_{z}(\eta) = Q^{-1} R_{\vec{n}}(\eta) Q$ or
$ R_{\vec{n}}(\eta) = Q R_{z}(\eta)  Q^{-1}$.
Once we find $\vec{d'}_k^i$ in the rotated frame,
we obtain $\vec{d}_k^i = Q \vec{d'}^{i}_k$.

## The first case 
First we consider the triclinic cell where there are two point groups $C_1 (1)$ and $C_i ({\overline 1})$.
For $C_1 (1)$ the symmetry is so low that
we simply propose to displace the $i$th atom in $x+$, $y+$, and $z+$, resulting in three independent displacements for 
a forward difference scheme. 
In this case, $V = \det d^i$ attains its largest possible value of unity.
For a central difference scheme, we need to displace the atom 
in $x-$, $y-$, and $z-$, resulting
a total of six displacements.
For $ C_i ({\overline 1})$, an
inversion operator exists that could map $\vec{d'}_k^i$ to $-\vec{d'}_k^i$ 
hence we need to do the same number of displacements (i.e., $3$)
for
both the forward and central difference schemes (see the third
 row of Table 1.

## The second case
For the monoclinic cell there are three point groups to consider.  We first
consider $C_2 (2)$, we propose to use a
first displacement $\vec{d'}^i_1 =
(0,-1,1)^T$. Under a 2-fold rotation around the $z$ axis,
we generate a second displacement
$\vec{d'}_2^i = (0, 1,1)^T$ from $\vec{d'}_1^i$.
We exhaust all operations in this point group,
hence we propose a second independent displacement vector
$\vec{d'}_3^i = (1,0,0)^T$.  
This choice makes $V$ to attain its largest possible value of unity.
For a central difference scheme, we need
to displace the $i$th atom in $-\vec{d'}_1^i$ independently
since there is no operation in the site symmetry 
that is able to map $\vec{d'}_1^i$ to its negative. 
However, the 2-fold rotation
is able to map $\vec{d'}_3^i = (1,0,0)^T$ to its negative $(-1,0,0)^T$. 
Hence the total number of
displacements is three for the point group $C_2(2)$ instead of four for a central difference scheme.  
For the second point group $C_{s} (m)$. There is no 
operation that can map $\vec{d'}_1^i$ and $\vec{d'}_3^i$
to their negatives. Hence the number of independent displacements is
four for a central difference scheme. 
The third point group to consider is $C_{2h}$. Since there is an inversion operator
the number of displacements
for the forward and central difference schemes are the same, which is two.

## The third case
For the next three point groups for the orthorhombic cells $D_2 (222)$,
$C_{2v} (mm2)$, and $D_{2h}(mmm)$, we simply need to displace an atom
in $\vec{d'}_1^i = \frac{1}{\sqrt{3}} (1,1,1)$ direction.  This direction is obtained by
maximizing a test vector $(t,u,v)^T$ and its images $(t,-u,-v)$ (under a
$2$-fold rotation around the $x$-axis) and $(-t,u,-v)$ (under a $2$-fold
rotation around the $y$-axis) thus forming a triple product of $ V = 4tuv$. The magnitude of the displacement vector
is constrained according to
$t^2+ u^2 +v^2 = 1$. A standard Lagrange multiplier method 
gives a possible solution 
$\vec{d'}_1^i = \frac{1}{\sqrt{3}} (1,1,1)^T$
that delivers the largest $V = \frac{4}{\sqrt{27}} = 0.7698$.

For $D_2$ and $C_{2v}$ there is no
operation that could map $\vec{d'}_1^i$ to $-\vec{d'}_1^i$ 
hence we need to do
a total of two
displacements for a central difference scheme.  
However, $D_{2h}$ has an inversion operator
that results in the same number (i.e., one)
of independent displacement
for both the forward and central difference schemes.

## The fourth case
For the next 24 point groups from $C_4 (4)$ to $O_h (m{\overline 3}m)$
covering tetragonal, trigonal, hexagonal, and cubic cells, we note
that we have either $3$ or ${\overline 3}$ or $4$ or ${\overline 4}$ that
is able to generate $\vec{d'}_2^i$ and $\vec{d'}_3^i$ from just one
$\vec{d'}_1^i$.
For these point groups, using the Lagrange multiplier method, we find
that a single displacement is able to give a maximum $V$ of
$\frac{4}{\sqrt{27}}$ and $1$ for $C_4$ and $C_3$, respectively,  with a
test vector of $\vec{d'}_1^i = (t,u,\frac{1}{\sqrt{3}})^T$, $t^2 + u^2
= \frac{2}{3}$.  In order to prepare $ \vec{d'}_1^i $ to be mapped under
an existing $C_2$ operation (which is unfortunately 
missing in eight point groups $C_4$, $S_4$, $C_{4v}$, $C_3$, 
$C_{3v}$, $C_6$, $C_{3h}$, and $C_{6v}$)
 to $ - \vec{d'}_1^i $
with the $C_2$ rotation axis pointing along $(n_x, n_y, n_z)$
direction, we can find the values of $t$ and $u$
by solving the simultaneous equations 
$t^2 + u^2 = \frac{2}{3}$ and $n_x t + n_y u + \frac{n_z}{\sqrt{3}} = 0$.
This is equivalent to solving a quadratic equation $(n_x^2 +
n_y^2) u^2 + \frac{2}{\sqrt{3}} n_y n_z u + \frac{1}{3} \left( n_z^2 -
2n_x^2\right) = 0$. 
To see more clearly the nature of solutions to the quadratic equation,
without loss of generality we assume the axis of rotation for
$C_2$ is in the $yz$ plane making an angle 
$\theta$ with the $z$ axis where $(n_x, n_y, n_z) = 
(0,\sin\theta,\cos\theta)$. This gives
two solutions $(t,u,v) = (\pm \sqrt{1- \frac{1}{3\sin^2\theta}}, - \frac{\cos\theta}{\sqrt{3} \sin \theta}, \sqrt{\frac{1}{3}})$.
 For $\theta = \frac{\pi}{2}$, $(t,u,v) = 
( \pm \sqrt{ \frac{2}{3}  }, 0, \sqrt{\frac{1}{3}})$.  
As $\theta$ is decreased from $\frac{\pi}{2}$, 
two solutions will approach one another along the circumference
of a circle and coincide with one another 
when $\theta $ becomes a critical
angle $\theta_c  = \sin^{-1} \sqrt{\frac{1}{3}} = 35.26^\circ$ 
and the solution is $(t,u,v) = ( 0, -\sqrt{\frac{2}{3}}, \sqrt{\frac{1}{3}}   )$.
When  $\theta < \theta_c$, there is no real solution.
 

# Error analysis for the force constants
Now we perform a detailed 
error analysis for the displacement method. We shall 
focus on $\Phi_{ij}$ which is
the $3\times 3$ force constant matrix block between
the $i$th atom and the $j$th atom
based on Eq.{eq}`eq:Phi`.
Since atom pairs are now understood, we 
suppress the
index $i$ in $d^i$ and index $j$ in $F^j$  for  forces in Eq. {eq}`eq:Phi`.
We let the exact force $\vec{F}_k$
to differ from the approximate force $\vec{f}_k$
by an error term $\lambda \vec{\epsilon}_k$ that
measures the intrinsic inaccuracies where

$$
\vec{f}_k = \vec{F}_k +  \lambda\vec{\epsilon}_k, \ \ k = 1, 2, 3
$$ (eq:f_k)

The calculated force constant block with $d$ 
is given by

$$
\Phi_d &=&
[ \vec{f}_1 | \vec{f}_2 | \vec{f}_3] (\lambda d)^{-1} \nonumber \\ 
&=& F(\lambda d)^{-1} + [\vec{\epsilon}_1 | \vec{\epsilon}_2 | \vec{\epsilon}_3]
d^{-1} \nonumber \\ 
&= & \Phi_e +  \epsilon d^{-1}
$$ (eq:Phid)

where $\Phi_e$ is the exact force constant matrix and
$\epsilon = [\vec{\epsilon}_1 | \vec{\epsilon}_2 | \vec{\epsilon}_3]$.
Eq.{eq}`eq:Phid`
clearly shows that the error can be attributed by the 
intrinsic inaccuracy of the force and 
the inverse of $d$, which could potentially be
large if $d$ is near singular.
In the absence of the exact force constant $\Phi_e$ to calculate errors
we rely on

$$
\phi_d  - \phi_{d_0} = \epsilon ( d^{-1} - d_0^{-1})
$$ (eq:Err)
where $d_0$ is chosen to be with the largest possible $V$.
We will study the inaccuracy of 
the calculated force constant as a function of $V$
in the next section.

# Results
To confirm the methodology outline above, we 
carry out density-functional theory (DFT) calculations within the local density
approximation, with projector augmented-wave (PAW) pseudopotential
as implemented in the Vienna Ab initio Simulation Package (VASP)
{cite}`Kresse96v6` on a generic system of silicon. 
We use a 2-atom primitive rhombohedral cell
and a $3\times 3\times 3$ supercell for the phonon calculations.
The cutoff energy for the plane-wave basis set is
307~eV. We use a $3\times 3\times 3$ 
$k$ mesh for the supercell calculations. 
The magnitude of the displacement $\lambda $ is $0.015$~\AA.
Because of the presence of a $C_3$ rotation with its rotation axis
along the $(1,1,1)^T$ direction in the Cartesian
coordinates and another $C_2$ 
operation in Si just one displacement is
sufficient for a central difference scheme for the phonon calculations. 
However, since we want to study the effect of $V$ on the 
phonon dispersions, we will still utilize $C_3$ in
the $(1,1,1)^T$ direction to rotate $\vec{d}_1 = (x,y,y)^T$
that results in just one displacement
if the forward difference scheme is used. However another
displacement in the negative 
direction is needed for a central difference scheme
since the $C_2$ operation is now no longer able to map 
$\vec{d}_1$ to $-\vec{d}_1$.
Under $C_3$ with its rotation axis in $(1,1,1)^T$, we have
$\vec{d}_2 = (y, x,y)$, and $ \vec{d}_3 = (y,y,x)$. This gives
$V = x^3 - 3xy^2 +2y^3$. With the constraint on the
magnitude of $\vec{d}_1$, which gives $2x^2 + y^2 = 1$,
we solve for
$x$ and $y$ using a Newton-Raphson scheme for a predetermined $V$. We notice that 
$V \rightarrow 0$ as
$x$ and $y$ go to $\frac{1}{\sqrt{3}}$. 

We investigate the error of the force constants
as measured by the Frobenius norm  $||\phi_d - \phi_{d_0}||_F $
as a function $V$ using Eq. {eq}`eq:Err`. We use $\det d_0 = 1$. 
Fig 1 shows the error of the force constants
linearly increases with decreasing $V = \det d$ .

For an estimate of $\epsilon $ in Eq. {eq}`eq:Err` without a full
knowledge of the error of the induced forces, we proceed
by assuming $\epsilon= \epsilon_0 J$, which
is characterized by a single
error parameter $\epsilon_0$ and an appropriate $J$ matrix. Our 
aim is to 
estimate  $\epsilon_0$ that measures the inaccuracy of the force by
fitting the RHS of Eq. {eq}`eq:Err`.
We observe that
in the limit when $V $ is very small, the three displacement vectors
$\vec{d}_1$, $\vec{d}_2$, and $\vec{d}_3$ become closer and closer 
to one another and converge to
$\frac{1}{\sqrt{3}}(1,1,1)^T$. 
To maximize error it is therefore reasonable to assume
$\vec{\epsilon}_1 = \epsilon_0 (1,-1,-1)^T$, 
$\vec{\epsilon}_2 = \epsilon_0 (-1,1,-1)^T$, and $
\vec{\epsilon}_3 = \epsilon_0 (-1,-1,1)^T$. Under such assumption and 
with Eq {eq}`eq:Err` we obtain a rather good fit and deduce $\epsilon_0 = 3.7\times 10^{-4} $~eV/\AA$^2$.
This in turn gives a reasonable estimate of the force inaccuracy of $\lambda \epsilon_0 = 5.6\times 10^{-6}$~eV/\AA.
Fig. 2 shows  phonon dispersions
with different $V$ values. 
Noticeable differences in phonon dispersions are observed when $V $ reaches $10^{-6}$.

```{figure} fig/2020-04-22-mpl-qphonon-det-gives-error-fit.png
The error (filled triangle) characterized
by the Frobenius norm of the force constant block $\phi_d  - \phi_{d_0}$ in Eq. {eq}`eq:Err`
as a function $V = \det d$. The fitting form according to the RHS of Eq. {eq}`eq:Err` gives
$\epsilon_0 = 3.7\times 10^{-4}$~eV/\AA$^2$.
```



```{figure} fig/2020-04-24-Si-Phonon-Det-1d-5-and-Det-1d-6.png
(a) The comparison of phonon dispersions obtained with $V = 10^{-5}$ and $V=1$.  (b) The comparison of phonon dispersions obtained with $V = 10^{-6}$ and $V=1$.
```

We select four representative systems to perform more phonon calculations.
For each system, we first carry out a phonon calculation by using
the symmetry-adapted atomic displacements.
Using identical computational parameters 
(i.e., energy cutoff, $k$ mesh, etc), we carry out a second phonon calculation by displacing
each inequivalent atom in the $x+$, $x-$, $y+$, $y-$, $z+$, and $z-$ along the Cartesian axes.
The results for hexagonal \mos{} [space group (SG) \# 194],  trigonal \bise{} (SG \# 166), 
orthorhombic \sbs{} (SG \# 62), and a 2D graphene sheet (SG \# 191) are shown
in Figs. 3, 4, 5, and 6, 
respectively. 
For all calculations the local density approximation is used, 
and the magnitude of the displacement $\lambda $ is $0.015$~\AA{}.
It is seen that the phonon dispersions 
obtained with both displacement methods are essentially the same that indicate
a correct implementation of the proposed methodology.
The efficacy of the symmetry-adapted atomic displacement method is shown in Table 2.


|        | symmetry-adapted | 6-displacement | All-displacemnt | 
|:-------|------------------|----------------|-----------------|
|MoS2    | 3                | 12             |              36 |
|Bi2Se3  | 5                | 18             |              30 |
|Sb2S3   | 20               | 30             |             120 |
|Graphene| 1                | 6              |              12 | 

Table 2:The numbers of atomic displacements required
for the symmetry-adapted, 6-displacement, and all-displacement methods.
The numbers of displacements for 6-displacement and all-displacement methods are $6 n_{\rm ineq}$ and $6 n_{\rm c}$, respectively.
Here $n_{\rm ineq}$ is the number of inequivalent atoms in the unit cell and $n_{\rm c}$ is the number of atoms in the unit cell.

```{figure} fig/2020-08-27-MoS2-Red-is-fmforces-Black-is-6displacement.png
Phonon dispersions of hexagonal \mos{} with $a=3.123$ and $c = 12.087$~\AA{} obtained with 
(a) symmetry-adapted atomic displacements and (b) six displacements for each of the two
inequivalent atoms.
The effect of longitudinal optical (LO) and transverse optical (TO) splitting{cite}`Liu14v16` has been taken account.
A $3\times 3 \times 2$ supercell is used. The cutoff energy for the plane-wave basis set
is $700$~eV. 
A $k$ mesh of $2 \times 2 \times 1$ is 
used for electronic relaxation. 
The selected $\vec{q}$ points (in $\vec{b}_1$, $\vec{b}_2$, and $\vec{b}_3$) are 
$\Gamma=[0,0,0]$, $M= [0,\frac{1}{2},0]$, $L=[0,\frac{1}{2},\frac{1}{2}]$, and $H=[\frac{1}{3},\frac{1}{3},\frac{1}{2}]$.
```

```{figure} fig/2020-08-28-Bi2Se3-fmforces-vs-6disp-LOTO.png
Phonon dispersions of trigonal \bise{} with $a_r=9.621$~\AA, 
$\alpha_r=24.64^\circ$, which corresponds to a conventional hexagonal cell of
$a= 4.105 $ and $c=27.973$~\AA{}, obtained with (a) symmetry-adapted atomic displacements and (b) six displacements for each of the three
inequivalent atoms.
The effect of longitudinal optical (LO) and transverse optical (TO) splitting{cite}`Liu14v16` has been taken account.
The supercell is $4 \times 4 \times 1$ of the 
conventional hexagonal unit cell. The cutoff energy for the plane-wave basis set is $423.2$~eV. A $k$ mesh of $4 \times 4 \times 2$ is 
used for electronic relaxation.
The selected $\vec{q}$ points (in $\vec{b}_1$, $\vec{b}_2$, and $\vec{b}_3$) are 
$\Gamma=[0,0,0]$, $L= [\frac{1}{2},0,0]$, $B=[\eta,\frac{1}{2},1-\eta]$, and $Z=[\frac{1}{2},\frac{1}{2},\frac{1}{2}]$, where $\eta = (1+4 \cos \alpha_r)/ (2+ 4\cos\alpha_r)$.{cite}`Setyawan10v49`
```

```{figure} fig/2020-08-27-Sb2S3-fmforces-vs-6disp-LOTO.png
Phonon dispersions of orthorhombic \sbs{} with $a=11.011$, $b= 3.812 $, 
and $c=10.794$~\AA{}  obtained with (a) symmetry-adapted atomic displacements and (b) six displacements for each of the five
inequivalent atoms.
The effect of longitudinal optical (LO) and transverse optical (TO) splitting{cite}`Liu14v16` has been taken account.
A $2\times 4 \times 2$ supercell is used. 
The cutoff energy for the plane-wave basis set is $323.3$~eV. A $k$ mesh of $2 \times 3 \times 2$ is 
used for electronic relaxation.
The selected $\vec{q}$ points (in $\vec{b}_1$, $\vec{b}_2$, and $\vec{b}_3$) are 
$\Gamma=[0,0,0]$, $X=[\frac{1}{2},0,0]$, $S=[\frac{1}{2},\frac{1}{2},0]$,
$R=[\frac{1}{2},\frac{1}{2},\frac{1}{2}]$,
$T=[0,\frac{1}{2},\frac{1}{2}]$, and $Z=[0,0,\frac{1}{2}]$.  
```

```{figure} fig/2020-08-27-Graphene-Red-fmforces-vs-Black-6displacements.png
Phonon dispersions of 2D graphene sheet with $a=2.462$~\AA{} and a vacuum height of $12$~\AA{} obtained
with  (a) a symmetry-adapted atomic displacement and (b) six displacements for the only
inequivalent atom.
A $4\times 4 \times 1$ supercell is used. The cutoff energy for the plane-wave basis set
is $500$~eV. A $k$ mesh of $6 \times 6 \times 1$ is
used for electronic relaxation.
The selected $\vec{q}$ points (in $\vec{b}_1$, $\vec{b}_2$, and $\vec{b}_3$) are
$\Gamma=[0,0,0]$, $M= [0,\frac{1}{2},0]$, and $K = [\frac{1}{3},\frac{1}{3},0]$.
```

# Conclusion
In summary, we have proposed 
a systematic displacement method that
guarantees a theoretical optimal volume $V$
spanned by the displacement
vectors to minimize possible severe roundoff errors
for lattice dynamical studies. 
The concepts of the star of $k$ and the group of $k$ that 
are originally defined in the reciprocal space have been
extended to real space 
to facilitate a practical implementation of the method that is
readily applied to all 32 crystallographic point groups 
and all 230 space groups. 
The method generates a
minimal set of irreducible displacement vectors to
keep the number of independent force calculations to a minimum.
The error in the calculated
force constants is shown to explicitly depend on the 
inverse of the volume $V$ spanned by the displacement directions and the 
intrinsic accuracy in the induced forces.
This justified the use of the Cartesian coordinates rather 
than the fractional coordinates to optimize $V$.
Several test systems have been used to illustrate the method.
The method is shown to be very effective in dealing
cells with a large aspect ratio due to a huge difference in lattice parameters, cells with a large 
vacuum separation, and 
cells that are very oblique due to an unconventional choice of a primitive cell.
We expect the strategies employed in this paper could be extended to
deal with higher order interatomic interactions for efficiency and accuracy.

# Acknowledgment

We acknowledge fruitful discussions with J. F. Kong.
We thank the National Supercomputing Center, Singapore (NSCC) and A*STAR Computational Resource Center, 
Singapore (ACRC) for computing resources. This work is
supported by RIE2020 Advanced Manufacturing and Engineering (AME) Programmatic Grant No A1898b0043.

# Appendix: Investigation of the determinant of $d$

We consider a generic example of \sbs{} crystal with
space group of  $Pnma$ (SG \# $62$) which has
an orthorhombic cell.
We shall investigate the determinant of $d$ (see Eq. {eq}`eq:dmat`)
as a function of cell dimensions $a$, $b$, and $c$. 
For a general discussion we
choose a primitive cell with $\vec{a}_1 = a\vec{i} + nb \vec{j}$, $\vec{a}_2 = b\vec{j}$, and
$\vec{a}_3 = c \vec{k}$, where $n$ is an integer. $\vec{i}$, $\vec{j}$, and $\vec{k}$ are
unit vectors in the Cartesian directions. First we 
consider $n=0$ that 
corresponds to the conventional primitive cell.

For \sbs{} crystal, all five inequivalent 
atoms occupy the $4c$ Wyckoff position hence
we may just consider any one of them. 
The site symmetry for the $4c$ position 
consists of just two elements:

\begin{equation}
I = 
  \begin{pmatrix}
       1 &  0 &  0 \\
       0 &  1 &  0 \\
       0 &  0 &  1 \\
  \end{pmatrix}
\end{equation}

and

\begin{equation}
m = 
  \begin{pmatrix}
       1 &  0 &  0 \\
       0 & -1 &  0 \\
       0 &  0 &  1 \\
  \end{pmatrix}
\end{equation}
which is a reflection across the $xz$ plane.

$\vec{g}_1$ may be chosen to be $(1,-1,0)^T$ which will be mapped to
$\vec{g}_2 = (1,1,0)^T$ under $m$. 
A second independent displacement is $\vec{g}_3 =
(0,0,1)^T$.

As discussed in the main text, the actual displacement in the Cartesian coordinates is 
$ \vec{d}_k^i =  \lambda \frac{A\vec{g_k}}{ | A\vec{g_k}   | }$, $k=1,2,3$. 
We find 

\begin{eqnarray}
V &=& \det d \\
&=& \frac{ \det A} { |A\vec{g}_1 | |A\vec{g}_2| |A\vec{g}_3|  } 
\\ &=&  \frac{  2 a b   }{  \sqrt{ a^2 + (n-1)^2 b^2    }   \sqrt{  a^2 + (n+1)^2 b^2   }   }
\end{eqnarray}
while the angle $\theta$ between 
the displacements $\vec{d}^i_1$ and $\vec{d}_2^i$ 
satisfies

\begin{equation}
\cos \theta = \frac{  a^2 + (n^2 -1 ) b^2   }{       \sqrt{ a^2 + (n-1)^2 b^2    }   \sqrt{  a^2 + (n+1)^2 b^2   }     }
\end{equation}

When $n=0$, we find $V = \frac{2ab}{a^2 + b^2}$ (therefore $V \rightarrow 2 (a/b)^{-1} $ for $a \gg b$)
 and $\cos \theta =
\frac{a^2 - b^2}{a^2 + b^2}$. The last result
is equivalent to a simpler expression
of $\tan^{-1} \frac{\theta}{2}  = \frac{b}{a}$, consistent
with the fact that the two displacements are $(a, -b,0)^T$ and $(a,b,0)^T$ in
the Cartesian coordinates.
If $a$ is much larger than $b$, then $V$ decreases as $\frac{2b}{a}$
to zero while $\cos \theta$ approaches $1$ since the two displacements
$\vec{d}_1^i$ and $\vec{d}_2^i$ 
are almost parallel.

Next we discuss the effect of $n \ne 0$. 
The site symmetry operations become{cite}`Nespolo16v72` 

$$
I = 
  \begin{pmatrix}
       1 &  0 &  0 \\
       0 &  1 &  0 \\
       0 &  0 &  1 \\
  \end{pmatrix}
$$
and

$$
m = 
  \begin{pmatrix}
       1 &  0 &  0 \\
       -2n & -1 &  0 \\
       0 &  0 &  1 \\
  \end{pmatrix}
$$

If $\vec{g}_1 = (1,0,0)^T$ is used, the $m$ operation on $\vec{g}_1$ gives
$\vec{g}_2 = (1,-2n,0)^T$. A second independent displacement is
$\vec{g}_3 = (0,0,1)^T$.
For a simple analysis, we set $a=b$.
For large $n$, we find $V$ approaches $-\frac{2}{n} $, while
$\cos\theta $ approaches $-1$, which means the displacement $\vec{d}_1^i$
approaches $-\vec{d}_2^i$.

This shows that using displacement directions in the form of nonzero $(e_1,e_2,e_3)^T$
where $e_i= 0,\pm 1$, $i=1,2,3$ may not be an optimal choice.

```{bibliography}
```
