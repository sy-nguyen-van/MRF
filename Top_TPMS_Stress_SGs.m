clear all; close all; clc;
%% source folders containing scripts not in this folder
addpath 'FE_routines'
addpath 'functions'
addpath 'mesh_utilities'
addpath 'optimization'
addpath 'utilities'
addpath 'plotting'
addpath 'Cubic_Spline'
addpath 'KeCe'
global OPT FE SGs
rho_max = 1;
OPT.LB = 0.001; OPT.UB = rho_max/1.875; OPT.dv0 = 0.015;
global examples TPMS P_norm
P_norm = 2 ;  % 1-No_SGs; 2-SGs_Elementwise
FE.dim = 2;
if FE.dim == 2;
    examples = 'Lbracket2d'; %'Lbracket2d'; 'cracked_plate2d';'doubleLbracket2d';  'Vframe2d';
else
    examples = 'cantilever3d';
end
%% Initialization
TPMS = 'Primitive';
%% Start timer

%% Initialization

init_SGs();

get_inputs();
init_FE();
init_optimization();
perform_analysis();
FD = 1;
if FD == 1
    %% Finite difference check of sensitivities
    run_finite_difference_check();
else
    if OPT.options.save_outputs
        out_folder = OPT.options.outputs_path;
        % Check if folder exists; if not, create:
        if ~exist(out_folder, 'dir')
            mkdir(out_folder);
        end
    end
    %% Optimization
    OPT.history = runmma(OPT.dv,@(x)obj(x),@(x)nonlcon(x));
    %%  Additional postprocessing
    % Compute and report compliance regardless of whether or not it is an
    % optimization function.
    % This allows us to compare the compliance of stress-constrained designs
    [c,~] = compute_compliance();
    fprintf('Compliance of final design = %-12.5e\n', c);
    OPT.c = c;
    
    % Report gray region fraction, which serves as an indication of
    % convergence to 0-1 design.
    fprintf('Gray region fraction of final design = %-12.5e\n', OPT.grf);
    
    fprintf('Largest true max of relaxed stress of final design = %-12.5e\n', ...
        OPT.true_stress_max);
    
    %  Save data structures to .mat file for future recovery
    [folder, baseFileName, ~] = fileparts(OPT.options.mat_output_path);
    mat_filename = fullfile(folder, strcat(baseFileName, '.mat'));
    save(mat_filename, 'OPT');
    
    %% Plot History
    if OPT.options.plot == true
        plot_history(2);
    end
    
    %% Report time
    % toc
    %% Turn off diary logging
    diary off;
    % ================================
    %% Copy outputs to selected folder
    % Update folder name as needed (as well as figure formats)
    if OPT.options.save_outputs
        if OPT.options.plot
            % Density plot
            saveas(1, strcat(out_folder, '/dens.fig'));
            saveas(1, strcat(out_folder, '/dens.pdf'));
            % History plot
            saveas(2, strcat(out_folder, '/hist.fig'));
            saveas(2, strcat(out_folder, '/hist.png'));
            if FE.dim == 2
                offset_fig_n = 2;
                for il=1:FE.nloads
                    saveas(offset_fig_n + il, strcat(out_folder, '/stress', string(il), '.fig'));
                    saveas(offset_fig_n + il, strcat(out_folder, '/stress', string(il), '.pdf'));
                end
            end
        end
%         % OPT data structure
%         save( strcat(out_folder, '/OPT.mat'), 'OPT');
%         % Diary
%         copyfile(diaryname, out_folder);
%         % Input files
%         copyfile(strcat('./', OPT.input_file_name), out_folder);
%         copyfile(FE.mesh_input.bcs_file, out_folder);
    end
    % ================================
end

