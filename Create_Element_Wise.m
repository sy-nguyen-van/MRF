clear all;clc; close all
%-------- N x 3 ==> 3 x 1 x N -------
sigma_x_true=[1:1:10]'; sigma_y_true = sigma_x_true; sigma_xy_true = sigma_x_true;
Data = [sigma_x_true, sigma_y_true, sigma_xy_true]';
sigma_true = permute(reshape(Data,[3,1,10]),[2,1,3]);
%-------- N x 3 ==> 3 x 3 x N -------
% sigma_x_Cases = ones(10,3); sigma_y_Cases = 2*ones(10,3); sigma_xy_Cases = 3*ones(10,3);
% Data = [sigma_x_Cases';sigma_y_Cases';sigma_xy_Cases'];
% sigma_micro = permute(reshape(Data,[3,3,10]),[2,1,3]);
% %----------
% Ce = ones(3,3); Be = ones(3,8); Le = ones(8,10); 
% CBL = Ce*Be*Le;
% Result2 = pagemtimes(sigma_micro,CBL);
% Result2 = sum(Result2,'all');





