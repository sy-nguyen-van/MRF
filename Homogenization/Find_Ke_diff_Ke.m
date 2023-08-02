clear all; close all; clc;
%% source folders containing scripts not in this folder
addpath(genpath('FE_routines'))
addpath(genpath('input_files'))
addpath(genpath('mesh_utilities'))
addpath(genpath('Data'))  
global examples TPMS
examples = 'Lbracket2d'; %'Lbracket2d'; 'cracked_plate2d';'doubleLbracket2d';  'Vframe2d'; 'cantilever3d'
%% Initialization
TPMS = 'Primitive';
get_inputs();
init_FE();
