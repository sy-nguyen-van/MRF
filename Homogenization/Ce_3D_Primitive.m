function Ce_3D = Ce_3D_Primitive(E,nu,thick)
%CE_3D_PRIMITIVE
%    CE_3D = CE_3D_PRIMITIVE(E,NU,THICK)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    26-Jul-2023 14:44:41

t2 = nu.*2.0;
t3 = nu+1.0;
t4 = nu-1.0;
t8 = thick.*2.153371655881386;
t10 = thick.*3.491053195219119;
t11 = thick.*3.60540812835199;
t5 = t2-1.0;
t6 = 1.0./t3;
t9 = exp(t8);
t12 = exp(t10);
t13 = exp(t11);
t7 = 1.0./t5;
t14 = t9.*4.613222557543144e-1;
t15 = t13.*1.662605565779405e-1;
t16 = t12.*1.826775953913206e-1;
t17 = t14-4.613222557543144e-1;
t18 = t15-1.662605565779405e-1;
t19 = t16-1.826775953913206e-1;
t20 = (E.*t6.*t17)./2.0;
t21 = E.*nu.*t6.*t7.*t18;
t23 = E.*t4.*t6.*t7.*t19;
t22 = -t21;
Ce_3D = reshape([t23,t22,t22,0.0,0.0,0.0,t22,t23,t22,0.0,0.0,0.0,t22,t22,t23,0.0,0.0,0.0,0.0,0.0,0.0,t20,0.0,0.0,0.0,0.0,0.0,0.0,t20,0.0,0.0,0.0,0.0,0.0,0.0,t20],[6,6]);