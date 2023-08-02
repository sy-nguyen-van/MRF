function vec_thick = vec_thick_cubic_spline(output_type, thick, Thick_Set_Train, index_interval, size_coeffs);
vec_thick = zeros(size_coeffs,1);
switch output_type
    case 'Predict'
    vec_thick((index_interval-1)*4+1,1)   = (thick-Thick_Set_Train(index_interval))^3; % t^3
    vec_thick((index_interval-1)*4+2,1)   = (thick-Thick_Set_Train(index_interval))^2; % t^2
    vec_thick((index_interval-1)*4+3,1)   = (thick-Thick_Set_Train(index_interval));   % t
    vec_thick((index_interval-1)*4+4,1)   = 1;
    case 'Gradient'
    vec_thick((index_interval-1)*4+1,1)   = 3*(thick-Thick_Set_Train(index_interval))^2; % 3*t^2
    vec_thick((index_interval-1)*4+2,1)   = 2*(thick-Thick_Set_Train(index_interval)); % 2*t
    vec_thick((index_interval-1)*4+3,1)   = 1; % 1
    vec_thick((index_interval-1)*4+4,1)   = 0; % 0
end
end