function FE_compute_element_stiffness()
% This function computes the element stiffness matrix fo all elements given
% an elasticity matrix.
% It computes the 'fully-solid' (i.e., unpenalized) matrix.
global FE OPT TPMS SGs
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
if FE.dim == 2
FE.Ce = zeros(3,3,FE.n_elem);
FE.dCe = zeros(3,3,FE.n_elem);

elseif FE.dim == 3
FE.Ce = zeros(6,6,FE.n_elem);   
FE.dCe = zeros(6,6,FE.n_elem);
end
% inititalize element stiffness matrices
FE.Ke = zeros(FE.n_edof,FE.n_edof,FE.n_elem);
% inititalize element stiffness matrices
FE.dKe = zeros(FE.n_edof,FE.n_edof,FE.n_elem);

% Material elasticity tensor
FE.SG_sig_x = zeros(SGs.n_ele_micro,3,FE.n_elem);
FE.SG_sig_y = zeros(SGs.n_ele_micro,3,FE.n_elem);
FE.SG_sig_xy = zeros(SGs.n_ele_micro,3,FE.n_elem);

FE.dSG_sig_x = zeros(SGs.n_ele_micro,3,FE.n_elem);
FE.dSG_sig_y = zeros(SGs.n_ele_micro,3,FE.n_elem);
FE.dSG_sig_xy = zeros(SGs.n_ele_micro,3,FE.n_elem);

%
FE.B0e = zeros(3*(FE.dim-1), FE.n_edof, FE.n_elem);
% loop over elements
bad_elem = false(FE.n_elem,1); % any element with a negative det_J will be flagged
%!!!!!!!!!!!!!!!!!!!!!!!!!
FE_compute_element_stiffness_update();


for e = 1:FE.n_elem
    if FE.dim == 2
        %!!!!!!!!!!!!!!!!!!!!!!!!!
        Ce = FE.Ce(:,:,e);
		
        %--------------------------------------------------------
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
                FE.Ke(:,:,e) = FE.Ke(:,:,e) + W(i)*W(j)*det_J * B'*Ce*B;
            end
        end
      J0 = Jacobian(0,0,e);
      inv_J0 = J0\eye(size(J0));
      GN  = inv_J0 * G0_N(0,0,e); % Shape function gradient at element centroid
      B0 = [GN(1,1) 0 GN(1,2) 0 GN(1,3) 0 GN(1,4) 0
                     0 GN(2,1) 0 GN(2,2) 0 GN(2,3) 0 GN(2,4)
                     GN(2,1) GN(1,1) GN(2,2) GN(1,2) GN(2,3) GN(1,3) GN(2,4) GN(1,4)];
      FE.B0e(:,:,e) = B0;
    elseif FE.dim == 3
        %!!!!!!!!!!!!!!!!!!!!!!!!!
        %!!!!!!!!!!!!!!!!!!!!!!!!!
        Ce = FE.Ce(:,:,e);
        %!!!!!!!!!!!!!!!!!!!!!!!!!
      % Voigt matrix for von Mises stress computation
      FE.V = [1 -1/2 -1/2 0 0 0; -1/2 1 -1/2 0 0 0; -1/2 -1/2 1 0 0 0; ...
              0 0 0 3 0 0; 0 0 0 0 3 0; 0 0 0 0 0 3];        
      loop over Gauss Points
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
                    FE.Ke(:,:,e) = FE.Ke(:,:,e) + W(i)*W(j)*W(k)*det_J * B'*Ce*B;
                end
            end
        end
        J0 = Jacobian8(0,0,0,e);
        inv_J0 = J0\eye(size(J0));            
        GN = inv_J0 * G0_N8(0,0,0,e);
        B0 = [GN(1,1) 0 0 GN(1,2) 0 0 GN(1,3) 0 0 GN(1,4) 0 0 GN(1,5) 0 0 GN(1,6) 0 0 GN(1,7) 0 0 GN(1,8) 0 0
                         0 GN(2,1) 0 0 GN(2,2) 0 0 GN(2,3) 0 0 GN(2,4) 0 0 GN(2,5) 0 0 GN(2,6) 0 0 GN(2,7) 0 0 GN(2,8) 0 
                         0 0 GN(3,1) 0 0 GN(3,2) 0 0 GN(3,3) 0 0 GN(3,4) 0 0 GN(3,5) 0 0 GN(3,6) 0 0 GN(3,7) 0 0 GN(3,8)
                         GN(2,1) GN(1,1) 0 GN(2,2) GN(1,2) 0 GN(2,3) GN(1,3) 0 GN(2,4) GN(1,4) 0 GN(2,5) GN(1,5) 0 GN(2,6) GN(1,6) 0 GN(2,7) GN(1,7) 0 GN(2,8) GN(1,8) 0
                         0 GN(3,1) GN(2,1) 0 GN(3,2) GN(2,2) 0 GN(3,3) GN(2,3) 0 GN(3,4) GN(2,4) 0 GN(3,5) GN(2,5) 0 GN(3,6) GN(2,6) 0 GN(3,7) GN(2,7) 0 GN(3,8) GN(2,8)
                         GN(3,1) 0 GN(1,1) GN(3,2) 0 GN(1,2) GN(3,3) 0 GN(1,3) GN(3,4) 0 GN(1,4) GN(3,5) 0 GN(1,5) GN(3,6) 0 GN(1,6) GN(3,7) 0 GN(1,7) GN(3,8) 0 GN(1,8)];
        FE.B0e(:,:,e) = B0;
    end
    if det_J < 0;bad_elem(e) = true; end
end
if sum(bad_elem) > 0; error('The following elements have nodes in the wrong order:\n%s',sprintf('%i\n',find(bad_elem))); end

