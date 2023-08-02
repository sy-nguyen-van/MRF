function [g,grad_g] = compute_P_norm()
global FE OPT SGs P_norm
sig_vm = zeros(FE.n_elem, FE.nloads);
p = OPT.parameters.aggregation_parameter;
row_index = [1:size(FE.edofMat(1,:),2)];
%==============2-P-NORM=========================
switch P_norm
    case 1
        F_adjoint= zeros(FE.n_global_dof,FE.nloads); % Struct for adjoint forces
        for iload = 1:FE.nloads
            U_global = FE.U(:,iload); % Global displacement vector
            for e=1:FE.n_elem
                Ce = FE.Ce(:,:,e); % Constitute matrix
                Be = FE.B0e(:,:,e); % Gradient of shape function
                Ue = U_global(FE.edofMat(e,:)); % Element displacement
                %----------------------------------------------------------
                FE.L_ga= zeros(8,FE.n_global_dof);
                col_index = FE.edofMat(e,:) ;
                ind = sub2ind(size(FE.L_ga),row_index,col_index);
                FE.L_ga(ind)=1;
                %----------------------------------------------------------
                sig = Ce*Be*Ue; % Element stress
                sig_vm(e,iload) = sqrt(sig(1)^2+sig(2)^2+3*sig(3)^2-sig(1)*sig(2)); %Von misig_vmss
                adj_2 = 0.5*(sig_vm(e,iload))^(p-2);
                adj_3 = (Ce*Be*FE.L_ga)';
                %----------------------------------------------------------
                adj_4 =[2*sig(1)-sig(2);2*sig(2)-sig(1);6*sig(3)];
                F_adjoint(:,iload)= F_adjoint(:,iload) + adj_2*adj_3*adj_4;
                %----------------------------------------------------------
            end
        end
        % =================== Compute sig_vm sensitivities ===========================
        FE.svm = sig_vm;
        g = sum(sig_vm.^p).^(1/p); %%% Constraint!!!!!!!!!!!!!!!!!!!!!!!!
        adj_1 = sum(sig_vm.^p).^(1/p-1);
        % Solve pseudo analysis to compute adjoint solution
        FE.dJdu = -adj_1.*F_adjoint;
        % Solve Lambda
        FE_solve('adjoint');
        %-----------------------------------------------------------------------
        grad_g = zeros(FE.n_elem,1);
        for iload=1:FE.nloads
            LAMDA = FE.lambda(:,iload);
            for e=1:FE.n_elem
                dKe = FE.dKe(:,:,e);
                dCe = FE.dCe(:,:,e);
                Ce = FE.Ce(:,:,e); % Constitute matrix
                Be = FE.B0e(:,:,e);
                
                Ue = U_global(FE.edofMat(e,:));
                sig = Ce*Be*Ue; % Element stress
                sig_vm(e,iload) = sqrt(sig(1)^2+sig(2)^2+3*sig(3)^2-sig(1)*sig(2));
                adj_5 = dCe*Be*Ue;
                adj_2 = 0.5*(sig_vm(e,iload))^(p-2);
                adj_4 =[2*sig(1)-sig(2);2*sig(2)-sig(1);6*sig(3)]';
                grad_g(e,1) = grad_g(e,1)+ LAMDA(FE.edofMat(e,:))'*dKe*Ue+  + adj_1(1,iload)*adj_2*adj_4*adj_5;
            end
        end
    case 2
        F_adjoint= zeros(FE.n_global_dof,FE.nloads); % Struct for adjoint forces
        for iload = 1:FE.nloads
            U_global = FE.U(:,iload); % Global displacement vector
            for e=1:FE.n_elem
                Ce = FE.Ce(:,:,e); % Constitute matrix
                Be = FE.B0e(:,:,e); % Gradient of shape function
                Ue = U_global(FE.edofMat(e,:)); % Element displacement
                %-----------------------
                Le= zeros(8,FE.n_global_dof);
                col_index = FE.edofMat(e,:) ;
                ind = sub2ind(size(Le),row_index,col_index);
                Le(ind)=1;
                %----------------------------------------------------------
                CBL = repmat((Ce*Be*Le)',[1,1,SGs.n_ele_micro]);
                %----------------------------------------------------------
                Traction = Ce*Be*Ue; % Element stress
                %----------------------------------------------------------
                %                 sig_x_Cases = FE.SG_sig_x(:,:,e);
                %                 sig_y_Cases = FE.SG_sig_y(:,:,e);
                %                 sig_xy_Cases = FE.SG_sig_xy(:,:,e);
                sig_x_true = FE.SG_sig_x(:,:,e)*Traction;
                sig_y_true = FE.SG_sig_y(:,:,e)*Traction;
                sig_xy_true = FE.SG_sig_xy(:,:,e)*Traction;
                %-------------Von-mises-Stress-----------------
                svm_SGs = sqrt(sig_x_true.^2 + sig_y_true.^2 + 3*sig_xy_true.^2 - sig_x_true.*sig_y_true);
                %----------------------------------------------------------
                sig_vm(e,iload) = (sum(svm_SGs.^p))^(1/p); %%% Constraint!!!!!!!!!!!!!!!!!!!!!!!!
                %----------------------------------------------------------
                %                 d_sig_n = permute(repmat(reshape([2*sig_x_true-sig_y_true, 2*sig_y_true-sig_x_true, 6*sig_xy_true]',[3,1,SGs.n_ele_micro]),[1,3,1]),[2,1,3]);
                %                 %----------------------------------------------------------
                %                 sig_micro = reshape([FE.SG_sig_x(:,:,e)';FE.SG_sig_y(:,:,e)';FE.SG_sig_xy(:,:,e)'],[3,3,SGs.n_ele_micro]);
                %                 %----------------------------------------------------------
                % %                 Result1 = pagemtimes( permute(pagemtimes(sig_micro,CBL),[2,1,3]),d_sig_n);
                %                 svm_n = repmat(reshape((svm_SGs.^(p-2))/2,[1,1,SGs.n_ele_micro]),[3,3,1]);
                %                 aa = sum((svm_n.*sig_micro).*d_sig_n,2);
                % %                 sss=sum(pagemtimes(CBL,aa),3);
                %                 F_adjoint(:,iload)= F_adjoint(:,iload) + sum(pagemtimes(CBL,aa),3);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %----------------------------------------------------------
                d_sig_n = permute(repmat(reshape([2*sig_x_true-sig_y_true, 2*sig_y_true-sig_x_true, 6*sig_xy_true]',[3,1,SGs.n_ele_micro]),[1,3,1]),[2,1,3]);
                %----------------------------------------------------------
                sig_micro = reshape([FE.SG_sig_x(:,:,e)';FE.SG_sig_y(:,:,e)';FE.SG_sig_xy(:,:,e)'],[3,3,SGs.n_ele_micro]);
                %----------------------------------------------------------
                %                 Result1 = pagemtimes( permute(pagemtimes(sig_micro,CBL),[2,1,3]),d_sig_n);
                svm_n = repmat(reshape((svm_SGs.^(p-2))/2,[1,1,SGs.n_ele_micro]),[3,3,1]);
                aa = repmat(permute(sum((svm_n.*sig_micro).*d_sig_n,2),[2,1,3]),[FE.n_global_dof,1,1]);
                %                 sss=sum(pagemtimes(CBL,aa),3);
                F_adjoint(:,iload)= F_adjoint(:,iload) + reshape(sum(sum(CBL.*aa,2),3),[FE.n_global_dof,1]);
                
            end
        end
        % =================== Compute sig_vm sensitivities ===========================
        FE.svm = sig_vm;
        g = sum(sig_vm.^p).^(1/p); %%% Constraint!!!!!!!!!!!!!!!!!!!!!!!!
        adj_2 = (sum(sig_vm.^p)).^(1/p-1);
        % Solve pseudo analysis to compute adjoint solution
        FE.dJdu = -adj_2.*F_adjoint;
        % Solve Lambda
        FE_solve('adjoint');
        %-----------------------------------------------------------------------
        grad_g = zeros(FE.n_elem,1);
        tic
        for iload=1:FE.nloads
            LAMDA = FE.lambda(:,iload);
            for e=1:FE.n_elem
                dKe = FE.dKe(:,:,e);
                dCe = FE.dCe(:,:,e);
                Ce = FE.Ce(:,:,e); % Constitute matrix
                Be = FE.B0e(:,:,e);
                Ue = U_global(FE.edofMat(e,:));
                %------------------------------------------------------------------
                dCeBeUe = repmat((dCe*Be*Ue)',[3,1,SGs.n_ele_micro]);
                %------------------------------------------------------------------
                Traction = Ce*Be*Ue; % Element stress
                %------------------------------------------------------------------
                %                 sig_x_Cases = FE.SG_sig_x(:,:,e);
                %                 sig_y_Cases = FE.SG_sig_y(:,:,e);
                %                 sig_xy_Cases = FE.SG_sig_xy(:,:,e);
                sig_x_true = FE.SG_sig_x(:,:,e)*Traction;
                sig_y_true = FE.SG_sig_y(:,:,e)*Traction;
                sig_xy_true = FE.SG_sig_xy(:,:,e)*Traction;
                Traction_Vec = repmat(Traction',[3,1,SGs.n_ele_micro]);
                %-------------Von-mises-Stress-----------------
                if FE.dim == 2
                    svm_SGs = sqrt(sig_x_true.^2 + sig_y_true.^2 + 3*sig_xy_true.^2 - sig_x_true.*sig_y_true);
                end
                %------------------------------------------------------------------
                %                 d_sig_x_Cases = FE.dSG_sig_x(:,:,e);
                %                 d_sig_y_Cases = FE.dSG_sig_y(:,:,e);
                %                 d_sig_xy_Cases = FE.dSG_sig_xy(:,:,e);
                %------------------------------------------------------------------
                d_sig_n = (repmat(reshape([2*sig_x_true-sig_y_true, 2*sig_y_true-sig_x_true, 6*sig_xy_true]',[3,1,SGs.n_ele_micro]),[1,3,1]));
                %------------------------------------------------------------------
                d_sig_micro = (permute(reshape([FE.dSG_sig_x(:,:,e)';FE.dSG_sig_y(:,:,e)';FE.dSG_sig_xy(:,:,e)'],[3,3,SGs.n_ele_micro]),[2,1,3]));
                %------------------------------------------------------------------
                sig_micro = permute(reshape([FE.SG_sig_x(:,:,e)';FE.SG_sig_y(:,:,e)';FE.SG_sig_xy(:,:,e)'],[3,3,SGs.n_ele_micro]),[2,1,3]);
                %------------------------------------------------------------------
                xx3 = reshape(   sum(sum( d_sig_n.*(d_sig_micro.*Traction_Vec + sig_micro.*dCeBeUe),1),2)  , [SGs.n_ele_micro,1] );
                grad_g(e) = grad_g(e) + LAMDA(FE.edofMat(e,:))'*dKe*Ue  + adj_2(1,iload)*sum(0.5*((svm_SGs.^(p-2))).*xx3,1);
            end
        end
end
% save these values in the OPT structure
slim = OPT.parameters.slimit; % Recall for multiple load cases this is a vector
h = FE.svm./slim;
% Compute aggregate stress function
% switch OPT.parameters.aggregation_type
%     case 'p-norm'
%         h = FE.svm./slim;
% %         dhds = ones(size(h))./slim;
%         phi = h;
% %         dphidh = ones(size(h));
%         P = OPT.parameters.aggregation_parameter;
%         [g, dgdphi] = smooth_max(reshape(phi,[],1), P, 'p-norm');
% %         dgdphi = reshape(dgdphi,size(h));
%         g = g - 1;
% end

% Account for filtering in sensitivities
% grad_g = OPT.H' * grad_g;

switch OPT.parameters.aggregation_type
    case 'p-norm'
        OPT.approx_h_max = g+1;
end
OPT.true_h_max = max(h,[],'all');
OPT.true_stress_max = max(FE.svm); % Vector with as many components as load cases
OPT.grad_stress = grad_g;






