function  [g,gradg] = compute_P_norm_SGs(dv)
global SGs
p = 8;
[sig, d_sig]= SGs_Stress(dv);
sig_x_Cases = sig(:,1:3);sig_y_Cases = sig(:,4:6);sig_xy_Cases = sig(:,7:9);
d_sig_x_Cases = d_sig(:,1:3);d_sig_y_Cases = d_sig(:,4:6);d_sig_xy_Cases = d_sig(:,7:9);
Traction = [1;1;1];
sig_x = sig_x_Cases*Traction;
sig_y = sig_y_Cases*Traction;
sig_xy = sig_xy_Cases*Traction;
svm = sqrt(sig_x.^2 + sig_y.^2 + 3*sig_xy.^2 - sig_x.*sig_y);
g = (sum(svm.^p))^(1/p);
adj = 0;
for e = 1: SGs.n_ele_micro
    d_sig_e = [d_sig_x_Cases(e,:); d_sig_y_Cases(e,:); d_sig_xy_Cases(e,:)];
    adj = adj + 0.5*(svm(e))^(p-2)*[2*sig_x(e)-sig_y(e);
        2*sig_y(e)-sig_x(e);
        6*sig_xy(e)]'*(d_sig_e*Traction);
end
gradg = (sum(svm.^p))^(1/p-1)*adj;
end