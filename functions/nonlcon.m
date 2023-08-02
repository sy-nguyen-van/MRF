function [g, geq, gradg, gradgeq] = nonlcon(dv)
global  OPT
OPT.dv_old = OPT.dv;
OPT.dv = dv(:);
FE_compute_element_stiffness_update(); % !!! UPDATE stiffness matrix !!!
perform_analysis();
n_con = OPT.functions.n_func-1; % number of constraints
g = zeros(n_con,1);
gradg = zeros(OPT.n_dv,n_con);
for i = 1:n_con    
    g(i) = OPT.functions.f{i+1}.value;
    g = g - OPT.functions.constraint_limit(i);
    gradg(:,i) = OPT.functions.f{i+1}.grad;
end
geq = [];
gradgeq = [];

end