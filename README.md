GaAs-Band-Structure
===================
In the name of God

This is a class project written for Quantum Transport course;
Which is based on the following book and it's matlab code appendix: 
  "Quantum Transport:Atom to Transistor" by Professor Supriyo Datta.

E(k) is calculated by finding the eigenvalues of the matrix in Eq. (5.3.3) for each value of k
along Gamma–X (that is, from k = 0 to k = 2pi/a * [1 0 0]) and
Gamma–L (that is, from k = 0 to k = pi/a * [1 1 1]) directions.
Some lines are degenerate.

-------------------------------------------------
GaAs Crystallographic System : ZinkBlend (2 FCCs)
---> 2 atoms in the primitive unit cell (1rst Brilluin Zone of the Direct Lattice)

Ga:[Ar]3d10,4s2,4p1 ---> Valance Orbirtals of Ga (approx) : 4s,4p,5s
As:[Ar]3d10,4s2,4p3 ---> Valance Orbirtals of As (approx) : 4s,4p,5s

---> Basis Functions:
|Sa>, |Xa>, |S*a>; |Sc>, |Xc>, |S*c>
      |Ya>               |Yc>
      |Za>               |Zc>

----------------
Orthogonality:
<_a|_a>=0
<_c|_c>=0

----------------
Diagonal Elements:
Esa=  <Sa|Sa>
Epa=  <Xa|Xa>
Es*a= <S*a|S*a>
Esc=  <Sc|Sc>
Epc=  <Xc|Xc>
Es*c= <S*c|S*c>

----------------
Overlap Integrals:
Ess=    <Sa|Sc>
Esapc=  <Sa|Xc>
0=      <Sa|S*c>
Epasc=  <Xa|Sc>
Exx=    <Xa|Xc>
Exy=    <Xa|Yc>
Epas*c= <Xa|S*c>
Es*apc= <S*a|Xc>

----------------
The 4 nearest cation atoms (Ga) place vectors: d1,d2,d3,d4.
---> g0,g1,g2,g3

The anion atom (As) placed inside the very tetrahedron.
---> Overlap Integrals (Ess, ... )

<n,Sa|m,Sc>=exp(ik.(dm-dn)) * <0,Sa|0,Sc> =<Sa|Sc>

On Page 118
