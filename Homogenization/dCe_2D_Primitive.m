function dCe_2D = dCe_2D_Primitive(E,nu,thick)
%DCE_2D_PRIMITIVE
%    DCE_2D = DCE_2D_PRIMITIVE(E,NU,THICK)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    26-Jul-2023 14:44:41

t2 = nu.*2.0;
t3 = nu+1.0;
t4 = nu-1.0;
t8 = thick.*3.491053195219119;
t9 = thick.*3.60540812835199;
t5 = t2-1.0;
t6 = 1.0./t3;
t10 = exp(t8);
t11 = exp(t9);
t7 = 1.0./t5;
t12 = E.*nu.*t6.*t7.*t11.*5.994371621104325e-1;
t14 = E.*t4.*t6.*t7.*t10.*6.377372030858152e-1;
t13 = -t12;
dCe_2D = reshape([t14,t13,0.0,t13,t14,0.0,0.0,0.0,E.*t6.*exp(thick.*2.153371655881386).*4.966991348843022e-1],[3,3]);
