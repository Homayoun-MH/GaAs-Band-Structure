%In the name of God, The Compassionate, The Merciful


% This is a class project written for Quantum Transport course;
% Which is based on the following book and it's matlab code appendix: 
  % "Quantum Transport:Atom to Transistor" by Professor Supriyo Datta.

% E(k) is calculated by finding the eigenvalues of the matrix in Eq. (5.3.3) for each value of k
% along Gamma–X (that is, from k = 0 to k = 2pi/a * [1 0 0]) and
% Gamma–L (that is, from k = 0 to k = pi/a * [1 1 1]) directions.
% Some lines on the plot are degenerate.

% -------------------------------------------------
% GaAs Crystallographic System : ZinkBlend (2 FCCs)
% ---> 2 atoms in the primitive unit cell (1rst Brilluin Zone of the Direct Lattice)

% Ga:[Ar]3d10,4s2,4p1 ---> Valance Orbirtals of Ga (approx) : 4s,4p,5s
% As:[Ar]3d10,4s2,4p3 ---> Valance Orbirtals of As (approx) : 4s,4p,5s

% ---> Basis Functions:
% |Sa>, |Xa>, |S*a>; |Sc>, |Xc>, |S*c>
      % |Ya>               |Yc>
      % |Za>               |Zc>

% ----------------
% Orthogonality:
% <_a|_a>=0
% <_c|_c>=0

% ----------------
% Diagonal Elements:
% Esa=  <Sa|Sa>
% Epa=  <Xa|Xa>
% Es*a= <S*a|S*a>
% Esc=  <Sc|Sc>
% Epc=  <Xc|Xc>
% Es*c= <S*c|S*c>

% ----------------
% Overlap Integrals:
% Ess=<Sa|Sc>
% Esapc=<Sa|Xc>
% 0=<Sa|S*c>
% Epasc=<Xa|Sc>
% Exx=<Xa|Xc>
% Exy=<Xa|Yc>
% Epas*c=<Xa|S*c>
% Es*apc=<S*a|Xc>

% ----------------
% The 4 nearest cation atoms (Ga) place vectors: d1,d2,d3,d4.
% ---> g0,g1,g2,g3

% The anion atom (As) placed inside the very tetrahedron.
% ---> Overlap Integrals (Ess, ... )

% <n,Sa|m,Sc>=exp(ik.(dm-dn)) * <0,Sa|0,Sc> =<Sa|Sc>

% On Page 118



clear all;
close all;
clc;

Esa=-8.3431; Epa=1.0414; Esc=-2.6569; Epc=3.6686; Esea=8.5914; Esec=6.7386;
Vss=-6.4513; Vxx=1.9546; Vxy=5.0779; Vsapc=4.4800; Vpasc=5.7839; Vseapc=4.8422; Vpasec=4.8077;
soa=.3787/3; soc=.0129/3;

d1=[0 0 0]/2;
d2=[0 -1 -1]/2;
d3=[-1 0 -1]/2;
d4=[-1 -1 0]/2;



Nt=100;
k_max=zeros(1,Nt+1);



% L:      pi*[1 1 1]
% Gamma:     [0 0 0]
% X:     2pi*[0 1 0]

ks = pi*[1 1 1;
         0 0 0;
         0 2 0];
 


for i=1:2
    
    k0 = ks(i,:);
    dk = (ks(i+1,:)-ks(i,:))/Nt;
    
    E = zeros(Nt+1,20); % Dimension of Hamiltonian is 20*20 ---> 20 eigen values.
    
    for Nk = 0:Nt
        
        k=k0+dk*Nk;
        
        p1=exp(1i*sum(k.*d1));
        p2=exp(1i*sum(k.*d2));
        p3=exp(1i*sum(k.*d3));
        p4=exp(1i*sum(k.*d4));
        
        g0=(p1+p2+p3+p4)/4;
        g1=(p1+p2-p3-p4)/4;
        g2=(p1-p2+p3-p4)/4;
        g3=(p1-p2-p3+p4)/4;
        
        h=[Esa/2 Vss*g0 0 0 0 Vsapc*g1 Vsapc*g2 Vsapc*g3 0 0;
            0 Esc/2 -Vpasc*conj(g1) -Vpasc*conj(g2) -Vpasc*conj(g3) 0 0 0 0 0;
            0 0 Epa/2 0 0 Vxx*g0 Vxy*g3 Vxy*g2 0 -Vpasec*g1;
            0 0 0 Epa/2 0 Vxy*g3 Vxx*g0 Vxy*g1 0 -Vpasec*g2;
            0 0 0 0 Epa/2 Vxy*g2 Vxy*g1 Vxx*g0 0 -Vpasec*g3;
            0 0 0 0 0 Epc/2 0 0 Vseapc*(g1) 0;
            0 0 0 0 0 0 Epc/2 0 Vseapc*(g2) 0;
            0 0 0 0 0 0 0 Epc/2 Vseapc*(g3) 0;
            0 0 0 0 0 0 0 0 Esea/2 0;
            0 0 0 0 0 0 0 0 0 Esec/2];
        H=[h+h' zeros(10);
            zeros(10) h+h'];
        
        hso=zeros(20);
        hso(3,4)=-1i*soa;hso(3,15)=soa;
        hso(4,15)=-1i*soa;
        hso(5,13)=-soa;hso(5,14)=1i*soa;
        hso(6,7)=-1i*soc;hso(6,18)=soc;
        hso(7,18)=-1i*soc;
        hso(8,16)=-soc;hso(8,17)=1i*soc;
        hso(13,14)=1i*soa;
        hso(16,17)=1i*soc;
        Hso=hso+hso';
        
        [V,D]=eig(H+Hso);
        eigst = diag(D);
        E(Nk+1,:) = sort(real(eigst));
        
    end
    
    k_abs=(0:Nt).*sqrt(dk*dk');
    plot(k_max + k_abs./pi,E,'b');
    hold on
    k_max =k_max + ones(1,Nt+1).*Nt*sqrt(dk*dk')/pi;
    
    axis( [0  k_max(Nt)  min(min(E))-2  max(max(E))+2] )
    xlabel('|k|/pi')
    ylabel('Energy (eV)')
    grid on
end
