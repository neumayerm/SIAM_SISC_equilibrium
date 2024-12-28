% Fast Numerical Techniques for Inverse Problems with Underlying Equilibrium Systems
% 
% Precomputation steps for an exemplary ECT sensor
%
% EMS 2022
% Contact: neumayer@tugraz.at
clear all, close all, clc

addpath .\lib_R

load SENSOR_FEM

profile on
for ii = 1:1
    FEM.x = 1+rand(size(FEM.x));
    [FEM] = func_solve_R(FEM);




    Deltax = rand(size(FEM.x));
    [DY,Dy] = func_Jop_R(FEM,Deltax);


    r = rand(size(FEM.y));
    [Jtr] = func_Jtop_R(FEM,r);

    J_avm = func_getJ_AVM_R(FEM);
    J_green = func_getJ_Green_R(FEM);
    J_Jop = func_getJ_Jop_R(FEM);
    J_Jtop = func_getJ_Jtop_R(FEM);
end
profile viewer