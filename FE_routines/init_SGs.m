function init_SGs()
global SGs TPMS FE
Scaler_Name = 'NoScaler';
%-------------------
Unit_Size = 1;
opti_methods = 'Cubic_Spline';
%-----------Train-Test-No of Thick/Poisson-----------------------
No_Thick_Train = 40;
No_Thick_Test =  15;
No_Thick_Pred = 10;
%----------------Range of Thickness/Poisson-----------------------
SGs.min_thick = 0.01; % Minimum thickness-depend on manufacturing's ability
SGs.max_thick = 0.05; % Maximum thickness
SGs.Thick_Set_Train = [linspace(SGs.min_thick, SGs.max_thick, No_Thick_Train)];
if FE.dim == 2
    Loads = {{'E11'}, {'E22'}, {'E12'}};  % X-Y Normal Traction; XY-Shear Traction
    n_case = length(Loads);
else
    Loads = {{'E11'}, {'E22'}, {'E33'}, {'E12'}, {'E13'}, {'E23'}};
    n_case = length(Loads);
end
SGs.n_case = n_case;
% name_model = append(TPMS, '_Manual');    % Model TPMS's name
% % --------------------------------------------------------
% for i_load = 1:n_case
%     Load_Case = string(Loads{i_load});
%     File_Name_Coeffs =  append(name_model, '_', Load_Case, '_', opti_methods,'_',Scaler_Name, '_', num2str(No_Thick_Train), '_Thick_Set_');
%     %----------------------------------------------------------------------
%     Coeffs_x_icase = load(append(File_Name_Coeffs, 'Coeffs_x', '.csv'));
%     index_iload = (i_load-1)*size(Coeffs_x_icase,1)+1:(i_load-1)*size(Coeffs_x_icase,1)+size(Coeffs_x_icase,1);
%     if i_load == 1
%        Coeffs_x = zeros(size(Coeffs_x_icase,1)*n_case, size(Coeffs_x_icase,2)); 
%        Coeffs_y = zeros(size(Coeffs_x_icase,1)*n_case, size(Coeffs_x_icase,2)); 
%        Coeffs_xy = zeros(size(Coeffs_x_icase,1)*n_case, size(Coeffs_x_icase,2)); 
%     end
%     Coeffs_x(index_iload,:) = Coeffs_x_icase; 
%     %----------------------------------------------------------------------
%     Coeffs_y_icase = readmatrix(append(File_Name_Coeffs, 'Coeffs_y', '.csv'));
%     Coeffs_y(index_iload,:) = Coeffs_y_icase;  
%     %----------------------------------------------------------------------
%     Coeffs_xy_icase = readmatrix(append(File_Name_Coeffs, 'Coeffs_xy', '.csv'));
%     Coeffs_xy(index_iload,:) = Coeffs_xy_icase; 
%     %----------------------------------------------------------------------    
% end
% save('Coeffs_x.mat', 'Coeffs_x') ;
% save('Coeffs_y.mat', 'Coeffs_y') ;
% save('Coeffs_xy.mat', 'Coeffs_xy') ;

%===========================
load Coeffs_x.mat;
load Coeffs_y.mat;
load Coeffs_xy.mat;
%--------------------------------------------------------
SGs.Coeffs_x = Coeffs_x;
SGs.Coeffs_y = Coeffs_y;
SGs.Coeffs_xy = Coeffs_xy;

% %-----------------------
SGs.ele_SGs = size(SGs.Coeffs_x,1)/SGs.n_case;
SGs.size_coeffs = size(SGs.Coeffs_x,2);
SGs.n_ele_micro = SGs.ele_SGs*4;
end

