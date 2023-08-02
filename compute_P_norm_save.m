function [g,grad_g] = compute_P_norm()
global FE OPT SGs
sigma_vm = zeros(FE.n_elem, FE.nloads);
p = OPT.parameters.aggregation_parameter;
row_index = [1:size(FE.edofMat(1,:),2)];
% Find A for adjoint method:
% K.Lamda = -A
    for iload = 1:FE.nloads
        Ul = FE.U(:,iload); % Global displacement vector
        F_adjoint= zeros(FE.n_global_dof,1); % Struct for adjoint forces
        for jj=1:FE.n_elem
            C0e = FE.material.C(:,:,jj); % Constitute matrix
            B0e = FE.B0e(:,:,jj); % Gradient of shape function
            Ue = Ul(FE.edofMat(jj,:)); % Element displacement
            %-----------------------
            FE.L_ga= zeros(8,FE.n_global_dof);
            col_index = FE.edofMat(jj,:) ;
            ind = sub2ind(size(FE.L_ga),row_index,col_index);
            FE.L_ga(ind)=1;
            %--------------------
            Traction_Vec = C0e*B0e*Ue; % Element stress
            [sigma_x_true, sigma_y_true, sigma_xy_true, sigma_x_Cases, sigma_y_Cases, sigma_xy_Cases, svm_SGs] = SGs_Stress(OPT.dv(jj), Traction_Vec');
            sigma_vm(jj,iload) = (sum(svm_SGs.^p))^(1/p); %%% Constraint!!!!!!!!!!!!!!!!!!!!!!!!
            adj_e = 0;
            for nn = 1:SGs.n_ele_micro
                sigma_micro = [ sigma_x_Cases(nn,:); sigma_y_Cases(nn,:); sigma_xy_Cases(nn,:)];
                adj_1 = (svm_SGs(nn,1)^(p-2)/2)*(sigma_micro*C0e*B0e*FE.L_ga)'*[2*sigma_x_true(nn,1)-sigma_y_true(nn,1);
                    2*sigma_y_true(nn,1)-sigma_x_true(nn,1);
                    6*sigma_xy_true(nn,1)];
                adj_e = adj_e + adj_1;
            end
            F_adjoint(:,iload)= F_adjoint(:,iload) + adj_e;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    % =================== Compute sigma_vm sensitivities ===========================
    FE.svm = sigma_vm;
    g = (sum(sigma_vm.^p))^(1/p); %%% Constraint!!!!!!!!!!!!!!!!!!!!!!!!
    adj_2 = (sum(sigma_vm.^p))^(1/p-1);
    % Solve pseudo analysis to compute adjoint solution
    FE.dJdu = -adj_2*F_adjoint;
    % Solve Lambda
    FE_solve('adjoint');
    %-----------------------------------------------------------------------
    grad_g = zeros(FE.n_elem,1);
    for iload=1:FE.nloads
        LAMDA = FE.lambda(:,iload);
        for e=1:FE.n_elem
            dK0e = FE.dKe(:,:,e);
            dC0e = FE.dC(:,:,e);
            C0e = FE.material.C(:,:,e); % Constitute matrix
            B0e = FE.B0e(:,:,e);
            Ue = Ul(FE.edofMat(e,:));
            %--------------------
            Traction_Vec = C0e*B0e*Ue; % Element stress
            [sigma_x_true, sigma_y_true, sigma_xy_true, sigma_x_Cases, sigma_y_Cases, sigma_xy_Cases, svm_SGs] = SGs_Stress(OPT.dv(e), Traction_Vec');
            [d_sigma_x_Cases, d_sigma_y_Cases, d_sigma_xy_Cases] = SGs_Stress_Gradient(OPT.dv(e));
            adjoint = 0;
            for nn = 1:SGs.n_ele_micro
                sigma_micro = [ sigma_x_Cases(nn,:); sigma_y_Cases(nn,:);sigma_xy_Cases(nn,:)];
                d_sigma_micro = [ d_sigma_x_Cases(nn,:);  d_sigma_y_Cases(nn,:); d_sigma_xy_Cases(nn,:)];
                adj_3 = (svm_SGs(nn,1)^(p-2)/2)*[2*sigma_x_true(nn,1)-sigma_y_true(nn,1);  2*sigma_y_true(nn,1)-sigma_x_true(nn,1); 6*sigma_xy_true(nn,1)]';
                adj_4 = (sigma_micro*dC0e*B0e*Ue + d_sigma_micro*Traction_Vec);
                adjoint = adjoint + adj_3*adj_4;
            end
            grad_g(e,iload)=  LAMDA(FE.edofMat(e,:))'*dK0e*Ue   + adj_2*adjoint;
        end
    end
    grad_g = sum(grad_g,2);
    
    OPT.approx_h_max = g;
    
    OPT.true_stress_max = max(FE.svm); % Vector with as many components as load casigma_vms
    OPT.grad_stress = grad_g;
    
    
    
    %=================================
    
    
    