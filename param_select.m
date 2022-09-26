% This file is part of the differentiation toolbox
% developed by Mohammad Rasool Mojallizadeh at INRIA, France
% Email: mohammad-rasool.mojallizadeh@inria.fr
%
% If you use this toolbox please cite the following report:
%
% Mohammad Rasool Mojallizadeh, Bernard Brogliato, Vincent Acary. Discrete-time differentiators:
% design and comparative analysis. 2020. hal-02960923
%
% This toolbox allows studying the behavior of several types of
% differentiators under different realistic conditions.
% 
% Copyright 2020 INRIA.
% 
% This program (the differentiation toolbox) is free software:
% you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%



% Randomized parameter tuning for parameter tuning
    
    % HDSTA Implicit
    mu_i=100*rand;    
    % HDSTA Explicit
    mu_e=100*rand;
    % Quadratic Explicit
    F_e=50*rand;
    alpha_e=50*rand;
    % Quadratic Implicit
    F_i=200*rand;
    alpha_i=200*rand;
    % Linear
    c=100*rand;
    % Toolbox
    r=30*rand;
    % ALIEN
    T_alien=5*rand;
    k_alien=floor(5*rand)+1;
    mu_alien=floor(5*rand)+1;
    %%%%%%%%%%%%%%%%%%%%
    L_E_STD=30*rand;
    L_I_STD=30*rand;
    L_SI_STD=30*rand;
    L_E_AO_STD=30*rand;
    L_I_AO_STD=30*rand;    
    L_SI_AO_STD=30*rand;
    L_E_HD_STD=30*rand;    
    L_I_HD_STD=30*rand;  
    L_EHDD=30*rand;
    L_EGHDD=30*rand;
    L_IHDD=30*rand;
    L_IGHDD=30*rand;    
    %%%%%%%%%%%%%%%%%%
    mu_VGED=20*rand;
    tau_VGED=15*rand;
    wc_VGED=15*rand+1;
    q_VGED=5*rand;
    
    L_SI_HD_STD=30*rand;
    mu_SI_HD_STD=100*rand;
    
    eps_E_STDAC=0.8*rand;
    alpha_E_STDAC=rand;
   
    %%%%%%%%%%%%%%%%%
    omegas_FDFF=100*rand;
    omegaf_FDFF=100*rand;
    ro_FDFF=100*rand;
    gamma_FDFF=100*rand;
    %%%%%%%%%%%%%%%%%
    
    F_AO_FDFF=100*rand;
    eps_AO_FDFF=100*rand;
    ws_AO_FDFF=200*rand;
    wf_AO_FDFF=200*rand;
    a1_AO_FDFF=500*rand;
    ro_AO_FDFF=500*rand;
    
    R_kalman=0.2*rand;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     
%     % HDSTA Implicit
%     mu_i=1000*rand;    
%     % HDSTA Explicit
%     mu_e=12*rand;
%     % Quadratic Explicit
%     F_e=5*rand;
%     alpha_e=50*rand;
%     % Quadratic Implicit
%     F_i=2000*rand;
%     alpha_i=2000*rand;
%     % Linear
%     c=1000000*rand;
%     % Toolbox
%     r=300000*rand;
%     % ALIEN
%     T_alien=5*rand;
%     k_alien=floor(5*rand)+1;
%     mu_alien=floor(5*rand)+1;
%     %%%%%%%%%%%%%%%%%%%%
%     L_E_STD=500*rand;
%     L_I_STD=500*rand;
%     L_SI_STD=500*rand;
%     L_E_AO_STD=500*rand;
%     L_I_AO_STD=500*rand;    
%     L_SI_AO_STD=500*rand;
%     L_E_HD_STD=5*rand;    
%     L_I_HD_STD=500*rand;  
%     L_EHDD=500*rand;
%     L_EGHDD=500*rand;
%     L_IHDD=500*rand;
%     L_IGHDD=500*rand;    
%     %%%%%%%%%%%%%%%%%%
%     mu_VGED=50*rand;
%     tau_VGED=50*rand;
%     wc_VGED=10*rand+1;
%     q_VGED=30*rand;
%     
%     L_SI_HD_STD=500*rand;
%     mu_SI_HD_STD=500*rand;
%     
%     eps_E_STDAC=10*rand;
%     alpha_E_STDAC=10*rand;  
%     R_kalman=10*rand;

