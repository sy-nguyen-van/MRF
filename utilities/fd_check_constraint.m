function fd_check_constraint()
%
% This function performs a finite difference check of the sensitivities of
% the CONSTRAINT function with respect to the bar design variables.
% It is currently setup for one constraint, but it can be easily modified
% for other/more constraints.
global OPT
% ===============================
% FINITE DIFFERENCE SENSITIVITIES
% ===============================
n_dv = OPT.n_dv;
n_con = OPT.functions.n_func-1;
grad_theta_i = zeros(n_dv, n_con);
fd_step = OPT.fd_step_size;
error = zeros(n_dv,n_con);
rel_error = zeros(n_dv,n_con);
dv_0 = OPT.dv;
dv_i = OPT.dv;
[theta_0, ~,grad_theta_0,~] = nonlcon(dv_i);
% Finite differences
disp('Computing finite difference sensitivities...');
% Do this for all design variables or only a few
for i = 1:n_dv
    dv_i(i) = dv_0( i ) + fd_step;
    [theta_i,~,~,~] = nonlcon(dv_i);
    grad_theta_i(i,:) = (theta_i - theta_0)/fd_step;
    error(i,:) = grad_theta_0(i,:) - grad_theta_i(i,:);
    rel_error(i,:) = error(i,:)./abs(theta_0');
    dv_i = dv_0;
end
OPT.dv = dv_0;
[max_error, ind_max_error] = max(error);
[max_rel_error, ind_max_rel_error] = max(rel_error);

disp('Max. ABSOLUTE error is:'); disp(max_error);
disp('It occurs at:');
disp('  variable:'); disp(ind_max_error);

disp('Max. RELATIVE error is:'); disp(max_rel_error);
disp('It occurs at:');
disp('  variable:'); disp(ind_max_rel_error);
OPT.grad_theta_i = grad_theta_i;
OPT.grad_theta_0 = grad_theta_0;
OPT.ratio = (OPT.grad_theta_i./OPT.grad_theta_0);
for j=1:n_con
    figure(2+j);    clf;    hold on
    plot(OPT.grad_theta_i(:,j),'o','LineStyle','-')
    plot(OPT.grad_theta_0(:,j),'.','LineStyle','-')
    legend('fd','analytical')
    title('Constraint function','Interpreter','latex',  'FontSize', 14, 'FontWeight','bold')
    legend('FD','Analytical',  'FontSize', 10, 'FontWeight','bold')
    xlabel('Design variable: v', 'FontSize', 12, 'FontWeight','bold')
    ylabel('dz/dv',  'FontSize', 12, 'FontWeight','bold')
    grid on;
end
