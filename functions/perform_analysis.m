function [] = perform_analysis()
%
% Filter and penalize densities, solve the finite
% element problem for the displacements and reaction forces, and then
% evaluate the relevant functions.
%
% Perform FE analysis

FE_analysis();
% Evaluate objective and constraint functions

evaluate_relevant_functions()

end

