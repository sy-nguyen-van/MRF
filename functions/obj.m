function [f, gradf] = obj(dv)
    global  OPT    
    OPT.dv_old = OPT.dv; % save the previous design
    OPT.dv = dv(:); % update the design
    
    FE_compute_element_stiffness_update(); % !!! UPDATE stiffness matrix !!!
    
    % Update or perform the analysis only if design has changed
    perform_analysis();
    %-------------------------------
    f = OPT.functions.f{1}.value;
    gradf = OPT.functions.f{1}.grad;
    gradf=reshape(gradf,[OPT.n_dv],1);
 
end
    
    