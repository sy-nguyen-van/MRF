function [f, gradf, g, geq, gradg, gradgeq] = obj_nonlcon(dv)
    global  OPT    
    OPT.dv_old = OPT.dv; % save the previous design
    OPT.dv = dv(:); % update the design
    FE_compute_element_stiffness_update(); % !!! UPDATE stiffness matrix !!!
    % Update or perform the analysis only if design has changed
    perform_analysis();
    %-------Evaluating Objective function---------
    f = OPT.functions.f{1}.value;
    gradf = OPT.functions.f{1}.grad;
    gradf=reshape(gradf,[OPT.n_dv],1);
    %-------Evaluating Constraint function---------
    n_con = OPT.functions.n_func-1; % number of constraints
    g = zeros(n_con,1);
    gradg = zeros(OPT.n_dv,n_con);
    for i = 1:n_con        
        g(i,1) = OPT.functions.f{i+1}.value;
        g = g - OPT.functions.constraint_limit(i);
        gradg(:,i) = OPT.functions.f{i+1}.grad;
    end
    geq = [];
    gradgeq = [];  
end
    
    