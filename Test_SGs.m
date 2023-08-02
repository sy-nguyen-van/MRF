clear all; clc; close all;
%================GET-INPUTS==============================
TPMS = 'Primitive';
Scaler_Name = 'MinMax';
addpath 'Cubic_Spline'
%-------------------
F_X = 10;
F_Y = 1;
F_Z = 10;
F_XY = 1;
F_XZ = 10;
F_YZ =1;
Unit_Size = 1;
opti_methods = 'Cubic_Spline';
thickness = 0.0105;
%-----------Train-Test-No of Thick/Poisson-----------------------
No_Thick_Train = 40;
No_Thick_Test =  15;
No_Thick_Pred = 10;
%----------------Range of Thickness/Poisson-----------------------
min_thick = 0.01; % Minimum thickness-depend on manufacturing's ability
max_thick = 0.05; % Maximum thickness
Thick_Set_Train = linspace(min_thick, max_thick, No_Thick_Train);
Loads = {{'E11'}, {'E22'}, {'E33'}, {'E12'}, {'E13'}, {'E23'}};
name_model = append(TPMS, '_Manual');    % Model TPMS's name
%=========CHECKING VALIDITY OF INPUTS=================
Traction_Vec = [F_X, F_Y, F_Z, F_XY, F_XZ, F_YZ];

% Find the index_interval using the find function
index_interval = find(Thick_Set_Train <= thickness, 1, 'last');
% Check if thickness is greater than or equal to max_thick
if thickness >= max_thick
    index_interval = index_interval - 1;
end
%--------------------------------------------------------
for i_load = 1:6
    Load_Case = string(Loads{i_load});
    File_Name_Coeffs =  append(name_model, '_', Load_Case, '_', opti_methods,'_',Scaler_Name, '_', num2str(No_Thick_Train), '_Thick_Set_');
    %==========================================================================
    Coeffs_x(:,:,i_load) = readmatrix( append(File_Name_Coeffs, 'Coeffs_x', '.csv'));    % Import xlsx file
    Coeffs_y(:,:,i_load) = readmatrix( append(File_Name_Coeffs, 'Coeffs_y', '.csv'));    % Import xlsx file
    Coeffs_xy(:,:,i_load) = readmatrix( append(File_Name_Coeffs, 'Coeffs_xy', '.csv'));   % Import xlsx file
end
%--------------------------------------------------------
for i_load = 1:6
    Load_Case = string(Loads{i_load});
    File_Name_Coeffs =  append(name_model, '_', Load_Case, '_', opti_methods,'_',Scaler_Name, '_', num2str(No_Thick_Train), '_Thick_Set_');
    %==========================================================================
    Scaler_x(:,:,i_load) = readmatrix( append(File_Name_Coeffs, 'Scaler_x', '.csv'));    % Import xlsx file
    Scaler_y(:,:,i_load) = readmatrix( append(File_Name_Coeffs, 'Scaler_y', '.csv'));    % Import xlsx file
    Scaler_xy(:,:,i_load) = readmatrix( append(File_Name_Coeffs, 'Scaler_xy', '.csv'));   % Import xlsx file
end
%--------------------------------------------------------
ele_SGs = size(Coeffs_x(:,:,1),1);
size_coeffs = size(Coeffs_x(:,:,1),2);
%--------------------------------------------------------
output_type = 'Predict';
vec_thick = vec_thick_cubic_spline(output_type, thickness, Thick_Set_Train, index_interval, size_coeffs);
vec_thick = repmat(vec_thick, [ele_SGs,1,6]);
%--------------------------------------------------------
pred_x = sum(Coeffs_x.*vec_thick,2);
pred_y  = sum(Coeffs_y.*vec_thick,2);
pred_xy  = sum(Coeffs_xy.*vec_thick,2);
%---------Scaling----------------------
pred_x =  reshape(pred_x.*(Scaler_x(:,2,:)-Scaler_x(:,1,:)) + Scaler_x(:,1,:), [ ele_SGs,6]);
pred_y =  reshape(pred_y.*(Scaler_y(:,2,:)-Scaler_y(:,1,:)) + Scaler_y(:,1,:), [ ele_SGs,6]);
pred_xy =  reshape(pred_xy.*(Scaler_y(:,2,:)-Scaler_xy(:,1,:)) + Scaler_xy(:,1,:), [ ele_SGs,6]);
%----------Linear combination---------------
pred_x = pred_x.*Traction_Vec;
pred_y = pred_y.*Traction_Vec;
pred_xy = pred_xy.*Traction_Vec;
%----------Superposition---------------
pred_x = sum(pred_x,2);
pred_y = sum(pred_y,2);
pred_xy = sum(pred_xy,2);
pred_svm = sqrt(pred_x.^2 + pred_y.^2 + pred_xy.^2 - pred_x.*pred_y);


