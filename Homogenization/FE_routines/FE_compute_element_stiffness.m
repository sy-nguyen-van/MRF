function FE_compute_element_stiffness()
%
% This function computes the element stiffness matrix fo all elements given
% an elasticity matrix.
% It computes the 'fully-solid' (i.e., unpenalized) matrix.

global FE TPMS

%% Inline Functions
% Jacobian matrix
Jacobian = @(xi,eta,elem) 0.25*[eta-1  xi-1
    1-eta -xi-1
    1+eta  1+xi
    -eta-1 1-xi]'*...
    FE.coords(:,FE.elem_node(:,elem))';
Jacobian8 = @(xi,eta,zeta, elem) 0.125*[-(1-zeta)*(1-eta)  -(1-zeta)*(1-xi)  -(1-eta)*(1-xi)
    (1-zeta)*(1-eta)  -(1-zeta)*(1+xi)  -(1-eta)*(1+xi)
    (1-zeta)*(1+eta)   (1-zeta)*(1+xi)  -(1+eta)*(1+xi)
    -(1-zeta)*(1+eta)   (1-zeta)*(1-xi)  -(1+eta)*(1-xi)
    -(1+zeta)*(1-eta)  -(1+zeta)*(1-xi)   (1-eta)*(1-xi)
    (1+zeta)*(1-eta)  -(1+zeta)*(1+xi)   (1-eta)*(1+xi)
    (1+zeta)*(1+eta)   (1+zeta)*(1+xi)   (1+eta)*(1+xi)
    -(1+zeta)*(1+eta)   (1+zeta)*(1-xi)   (1+eta)*(1-xi)]'* ...
    FE.coords(:,FE.elem_node(:,elem))';

% Gradient of shape function matrix in parent coordinates
G0_N = @(xi,eta,elem) 0.25*[eta-1  xi-1
    1-eta -xi-1
    1+eta  1+xi
    -eta-1 1-xi]';
G0_N8 = @(xi,eta,zeta,elem) 0.125*[-(1-zeta)*(1-eta)  -(1-zeta)*(1-xi)  -(1-eta)*(1-xi)
    (1-zeta)*(1-eta)  -(1-zeta)*(1+xi)  -(1-eta)*(1+xi)
    (1-zeta)*(1+eta)   (1-zeta)*(1+xi)  -(1+eta)*(1+xi)
    -(1-zeta)*(1+eta)   (1-zeta)*(1-xi)  -(1+eta)*(1-xi)
    -(1+zeta)*(1-eta)  -(1+zeta)*(1-xi)   (1-eta)*(1-xi)
    (1+zeta)*(1-eta)  -(1+zeta)*(1+xi)   (1-eta)*(1+xi)
    (1+zeta)*(1+eta)   (1+zeta)*(1+xi)   (1+eta)*(1+xi)
    -(1+zeta)*(1+eta)   (1+zeta)*(1-xi)   (1+eta)*(1-xi)]';

%% Solid Stiffness Matrix Computation
% We will use a 2 point quadrature to integrate the stiffness matrix:
gauss_pt = [-1 1]/sqrt(3);
W = [1 1 1];
num_gauss_pt = length(gauss_pt);

% Material elasticity tensor
% C = FE.material.C;
syms E nu thick real

rho = 1.875*thick;
% inititalize element stiffness matrices
FE.Ke = sym(zeros(FE.n_edof,FE.n_edof,FE.n_elem));
% Store strain-displacement matrix at centroid, which is used to compute
% stresses
%
FE.B0e = sym(zeros(3*(FE.dim-1), FE.n_edof, FE.n_elem));
% loop over elements
% bad_elem = false(FE.n_elem,1); % any element with a negative det_J will be flagged
e=1;
if FE.dim == 2
    if TPMS == 'Primitive'
        C = Ce_2D_Primitive(E,nu,thick);
    end
    % Voigt matrix for von Mises stress computation
    FE.V = [1 -1/2 0; -1/2 1 0; 0 0 3];
    % loop over Gauss Points
    for i = 1:num_gauss_pt
        xi = gauss_pt(i);
        for j = 1:num_gauss_pt
            eta = gauss_pt(j);
            % Compute Jacobian
            J = Jacobian(xi,eta,e);
            det_J = det(J);
            inv_J = J\eye(size(J));
            % Compute shape function derivatives (strain displacement matrix)
            GN = inv_J * G0_N(xi,eta,e);
            B = [GN(1,1) 0 GN(1,2) 0 GN(1,3) 0 GN(1,4) 0
                0 GN(2,1) 0 GN(2,2) 0 GN(2,3) 0 GN(2,4)
                GN(2,1) GN(1,1) GN(2,2) GN(1,2) GN(2,3) GN(1,3) GN(2,4) GN(1,4)];
            FE.Ke(:,:,e) = FE.Ke(:,:,e) + W(i)*W(j)*det_J * B'*C*B;
        end
    end
    %=====================================
    Ke_2D = FE.Ke(:,:,e);
    dKe_2D = diff(Ke_2D,thick);
    matlabFunction(Ke_2D,"File",append('Ke_2D_',TPMS));
    matlabFunction(dKe_2D,"File", append('dKe_2D_',TPMS));
elseif FE.dim == 3
    if TPMS == 'Primitive';
        C = Ce_3D_Primitive(E,nu,thick);
    end
    % Voigt matrix for von Mises stress computation
    FE.V = [1 -1/2 -1/2 0 0 0; -1/2 1 -1/2 0 0 0; -1/2 -1/2 1 0 0 0; ...
        0 0 0 3 0 0; 0 0 0 0 3 0; 0 0 0 0 0 3];
    % loop over Gauss Points
    for i = 1:num_gauss_pt
        xi = gauss_pt(i);
        for j = 1:num_gauss_pt
            eta = gauss_pt(j);
            for k = 1:num_gauss_pt
                zeta = gauss_pt(k);
                % Compute Jacobian
                J = Jacobian8(xi,eta,zeta,e);
                det_J = det(J);
                inv_J = J\eye(size(J));
                % Compute shape function derivatives (strain displacement matrix)
                GN = inv_J * G0_N8(xi,eta,zeta,e);
                B = [GN(1,1) 0 0 GN(1,2) 0 0 GN(1,3) 0 0 GN(1,4) 0 0 GN(1,5) 0 0 GN(1,6) 0 0 GN(1,7) 0 0 GN(1,8) 0 0
                    0 GN(2,1) 0 0 GN(2,2) 0 0 GN(2,3) 0 0 GN(2,4) 0 0 GN(2,5) 0 0 GN(2,6) 0 0 GN(2,7) 0 0 GN(2,8) 0
                    0 0 GN(3,1) 0 0 GN(3,2) 0 0 GN(3,3) 0 0 GN(3,4) 0 0 GN(3,5) 0 0 GN(3,6) 0 0 GN(3,7) 0 0 GN(3,8)
                    GN(2,1) GN(1,1) 0 GN(2,2) GN(1,2) 0 GN(2,3) GN(1,3) 0 GN(2,4) GN(1,4) 0 GN(2,5) GN(1,5) 0 GN(2,6) GN(1,6) 0 GN(2,7) GN(1,7) 0 GN(2,8) GN(1,8) 0
                    0 GN(3,1) GN(2,1) 0 GN(3,2) GN(2,2) 0 GN(3,3) GN(2,3) 0 GN(3,4) GN(2,4) 0 GN(3,5) GN(2,5) 0 GN(3,6) GN(2,6) 0 GN(3,7) GN(2,7) 0 GN(3,8) GN(2,8)
                    GN(3,1) 0 GN(1,1) GN(3,2) 0 GN(1,2) GN(3,3) 0 GN(1,3) GN(3,4) 0 GN(1,4) GN(3,5) 0 GN(1,5) GN(3,6) 0 GN(1,6) GN(3,7) 0 GN(1,7) GN(3,8) 0 GN(1,8)];
                FE.Ke(:,:,e) = FE.Ke(:,:,e) + W(i)*W(j)*W(k)*det_J * B'*C*B;
            end
        end
    end
    %=====================================
    Ke_3D = FE.Ke(:,:,e);
    dKe_3D = diff(Ke_3D,thick);
    matlabFunction(Ke_3D,"File",append('Ke_3D_',TPMS));
    matlabFunction(dKe_3D,"File", append('dKe_3D_',TPMS));
end



