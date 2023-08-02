function SGs_Stress(dv)
global SGs FE
vec_thick = zeros(SGs.size_coeffs,FE.n_elem);
d_vec_thick = zeros(SGs.size_coeffs,FE.n_elem);
% Calculate the index_interval for all elements in dv at once
index_interval = sum(SGs.Thick_Set_Train <= dv, 2);
% Check if thickness is greater than or equal to max_thick
index_interval(dv >= SGs.max_thick) = index_interval(dv >= SGs.max_thick) - 1;
% Calculate the cubic spline basis element-wise
thick_cubic = dv - [SGs.Thick_Set_Train(index_interval)]';
index_offset = 4*(index_interval' - 1);
col = [1:1:FE.n_elem];
vec_thick(sub2ind(size(vec_thick), index_offset + 1, col)) = thick_cubic.^3; % t^3
vec_thick(sub2ind(size(vec_thick), index_offset + 2, col))  = thick_cubic.^2; % t^2
vec_thick(sub2ind(size(vec_thick), index_offset + 3, col))  = thick_cubic;    % t
vec_thick(sub2ind(size(vec_thick), index_offset + 4, col))  = 1;
% Calculate the gradient element-wise
d_vec_thick(sub2ind(size(vec_thick), index_offset + 1, col)) = 3*thick_cubic.^2; % 3t^2
d_vec_thick(sub2ind(size(vec_thick), index_offset + 2, col)) = 2*thick_cubic;    % 2t
d_vec_thick(sub2ind(size(vec_thick), index_offset + 3, col)) = 1;                % 1
d_vec_thick(sub2ind(size(vec_thick), index_offset + 4, col)) = 0;                % 0
%--------------------------------------------------------
sig_x = SGs.Coeffs_x*vec_thick;
sig_y  = SGs.Coeffs_y*vec_thick;
sig_xy  = SGs.Coeffs_xy*vec_thick;
%-------------------------------------------
sig_x_Cases = reshape(sig_x, [SGs.ele_SGs, 3,  FE.n_elem]);
sig_y_Cases = reshape(sig_y, [SGs.ele_SGs, 3,  FE.n_elem]);
sig_xy_Cases = reshape(sig_xy, [SGs.ele_SGs, 3, FE.n_elem]);
%-------------------------------------------
[FE.SG_sig_x,FE.SG_sig_y,FE.SG_sig_xy] = Stress_Trans(sig_x_Cases,sig_y_Cases,sig_xy_Cases);
%----------Linear combination-Superposition--------------
d_sig_x = SGs.Coeffs_x*d_vec_thick;
d_sig_y  = SGs.Coeffs_y*d_vec_thick;
d_sig_xy  = SGs.Coeffs_xy*d_vec_thick;
%-------------------------------------------
d_sig_x_Cases = reshape(d_sig_x, [SGs.ele_SGs, 3,  FE.n_elem]);
d_sig_y_Cases = reshape(d_sig_y, [SGs.ele_SGs, 3,  FE.n_elem]);
d_sig_xy_Cases = reshape(d_sig_xy, [SGs.ele_SGs, 3,  FE.n_elem]);
%-------------------------------------------
[FE.dSG_sig_x,FE.dSG_sig_y,FE.dSG_sig_xy] = Stress_Trans(d_sig_x_Cases,d_sig_y_Cases,d_sig_xy_Cases);
%------------------------------------------
end

