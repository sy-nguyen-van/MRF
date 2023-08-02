function FE_compute_element_stiffness_update()
global FE OPT TPMS P_norm

if FE.dim == 2
    %!!!!!!!!!!!!!!!!!!!!!!!!!
    switch  TPMS
        case'Primitive'
            Ce_2D_Primitive(OPT.dv(:));
            %------------
            dCe_2D_Primitive(OPT.dv(:));
            %---------
            Ke_2D_Primitive(OPT.dv(:));
            %---------
            dKe_2D_Primitive(OPT.dv(:));
            %-----------------------------
            if P_norm >1
            SGs_Stress(OPT.dv(:));
            end
    end
end


