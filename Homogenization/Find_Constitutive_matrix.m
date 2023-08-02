clear all;clc
syms thick E nu real
TPMS = 'Primitive';
% Compute constitutive matrix
dim = 2;
Coeffs = load(append(TPMS,'_Coeffs_Homogenization.csv'));
E_nu = (E/((1 + nu)*(1 - 2*nu)));
%-----------------------------------
rho = 1.875*thick;
%-----------------------------------
a = E_nu*(1-nu);
b = E_nu*nu;
c = E_nu*(1-2*nu)/2;

a = a*(Coeffs(1,1)*exp(Coeffs(1,2)*rho)-Coeffs(1,1));
b = b*(Coeffs(2,1)*exp(Coeffs(2,2)*rho)-Coeffs(2,1));
c = c*(Coeffs(3,1)*exp(Coeffs(3,2)*rho)-Coeffs(3,1));
%------------------------------------
C2 = [
    a, b, 0;
    b, a, 0;
    0, 0, c];
%------------------------------------
C3 = [
    a,  b,  b,  0,  0,  0;
    b,  a,  b,  0,  0,  0;
    b,  b,  a,  0,  0,  0;
    0,  0,  0,  c,  0,  0;
    0,  0,  0,  0,  c,  0;
    0,  0,  0,  0,  0,  c];
Ce_2D = 0.5*(C2+C2');
Ce_3D = 0.5*(C3+C3');
dCe_2D = diff(Ce_2D, thick);
dCe_3D = diff(Ce_3D, thick);
%--------------
matlabFunction(Ce_2D,"File", append('Ce_2D_',TPMS));
matlabFunction(Ce_3D,"File",append('Ce_3D_',TPMS));
%--------------
matlabFunction(dCe_2D,"File", append('dCe_2D_',TPMS));
matlabFunction(dCe_2D,"File",append('dCe_3D_',TPMS));
