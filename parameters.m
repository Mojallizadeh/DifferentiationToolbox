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



%% Parameters of the differentiators and the simulation

%% Simulation parameters
t_f=10; % Simulation time

%% General parameters
% For sample_length==1 > Only the plot of error will be shown
% For sample_length>1, the cost functions will be shown with 
% respect to "changing_parameter"

sample_length=1; % The cost functions should be calculated for
                   % how many different conditions?
diff_order=1; % Order of differentiation

SNR=30*ones(1,sample_length);
% SNR=linspace(30,120,sample_length);
% Noise type
noise_type="white"; % sin, white, bell
Wn_vector=20*ones(1,sample_length); % Frequency of the sinusoidal noise
% Wn_vector=linspace(0,50,sample_length); % Frequency of the sinusoidal noise
An_vector=0.05*ones(1,sample_length); % Amplitude of the sinusoidal noise
% An_vector=linspace(0,0.5,sample_length); % Amplitude of the sinusoidal noise

h_vector=50e-3*ones(1,sample_length);
% h_vector=linspace(1e-4,100e-3,sample_length);
w_vector=1*ones(1,sample_length); % Frequency of the input signal: can be a vector
% w_vector=linspace(0.1,10,sample_length);

% Quantization?
Quant=0;

% Steady-state performances? Except fot L_infty
ssp=0;

hc=1e-4; % sampling period of the interpolation: hc > min(h)
changing_parameter="SNR"; % h, w, w/h, h/w, param, SNR An Wn
plot_type="N"; % N=normal, Y=semilogy, X=semilogx, XY=loglog

% Initial value of the cost functions
Global=inf*ones(24);
k1=10000;k2=0;k3=0;k4=0;k5=0; % k1=L2bar, k2=L2tilde, k3=Linf, k4=VAR, k5=THD
%  k1=0;k2=0;k3=0;k4=0;k5=1; % k1=L2bar, k2=L2tilde, k3=Linf, k4=VAR, k5=THD

AO_order=3; % Order of the AO differentiator
%% Input signal parameter
amp=1;  % Amplitude of the input signal
phi=0;  % Phase of the input signal
%% Initial condition
x0=[0 1 0 -1 0 1 0 -1]'; % zero init error for sin(t)
% x0=[0 0 0 0 0 0 0 0]';
%%
load parameters_file % Load the optimal parameters

% % %% Custom parameters
%   % HDSTA Implicit
%     mu_i=10;    
%     % HDSTA Explicit
%     mu_e=10;
%     % Quadratic Explicit
%     F_e=100;
%     alpha_e=100;
%     % Quadratic Implicit
%     F_i=6.7806;
%     alpha_i=149.9364;
%     % Linear
%     c=100;
%     % Toolbox
%     r=100;
%     % ALIEN
%     T_alien=100;
%     k_alien=100;
%     mu_alien=100;
%     %%%%%%%%%%%%%%%%%%%%
%     L_E_STD=100;
%     L_I_STD=0.7646;
%     L_SI_STD=100;
%     L_E_AO_STD=100;
%     L_I_AO_STD=100;   
%     L_SI_AO_STD=100;
%     L_E_HD_STD=100;    
%     L_I_HD_STD=100;  
%     L_EHDD=100;
%     L_EGHDD=100;
%     L_IHDD=100;
%     L_IGHDD=100;   
%     %%%%%%%%%%%%%%%%%%
%     mu_VGED=100;
%     tau_VGED=100;
%     wc_VGED=100;
%     q_VGED=100;
%     
%     L_SI_HD_STD=100;
%     mu_SI_HD_STD=100;
%     
%     eps_E_STDAC=100;
%     alpha_E_STDAC=100;

% R_kalman=0.000785460630004231920298085611875649192370474338531494140625; % Quant

% R_kalman=0.0000000000000000000000000000000087388370254397048797664760087879992861665356596832258499419629959627;

% R_kalman=1e-6;
% R_kalman=3.5920e-08;
% R_kalman=3.8107e-09;  %Third-order

% R_kalman=7.8546e-4;

% 
% omegas_FDFF=99.9819;
% omegaf_FDFF=564.3257;
% gamma_FDFF=194.0534;
% ro_FDFF=659.3234;
% 
% %%%%%%%%%%%%%%%
% 
% F_AO_FDFF=198.8190;
% eps_AO_FDFF=887.4548;
% ws_AO_FDFF=5.4169;
% wf_AO_FDFF=1684.9944;
% a1_AO_FDFF=4359.4855;
% ro_AO_FDFF=4451.7299;
% 



% %%%%%%%%%%%%%%%%%  no noise- test h
% % %% Custom parameters
%   % HDSTA Implicit
%     mu_i=202.4230;    
%     % HDSTA Explicit
%     mu_e=410.1379;
%     % Quadratic Explicit
%     F_e=5.1201;
%     alpha_e=3.5689;
%     % Quadratic Implicit
%     F_i=1999.7762;
%     alpha_i=41.8930;
%     % Linear
%     c=999991;
%     % Toolbox
%     r=86.9927;
%     % ALIEN
%     T_alien=0.0341;
%     k_alien=4;
%     mu_alien=4;
%     %%%%%%%%%%%%%%%%%%%%
%     L_E_STD=0.8383;
%     L_I_STD=137.0017;
%     L_SI_STD=0.8338;
%     L_E_AO_STD=1.4859;
%     L_I_AO_STD=0.6714;   
%     L_SI_AO_STD=1.0277;
%     L_E_HD_STD=0.0410;    
%     L_I_HD_STD=181.0637;  
%     L_EHDD=1.3011;
%     L_EGHDD=1.3438;
%     L_IHDD=0.6346;
%     L_IGHDD=0.6468;   
% 
%     
%     mu_VGED=1.4324;
%     tau_VGED=138.3938;
%     wc_VGED=17.2907;
%     q_VGED=2.8214;
%     
%     L_SI_HD_STD=34.6716;
%     mu_SI_HD_STD=0.9485;
%     
%     eps_E_STDAC=0.01;
%     alpha_E_STDAC=0.8270;
% 
%     R_kalman=1e-6;
% 
% 
% omegas_FDFF=999.9819;
% omegaf_FDFF=564.3257;
% gamma_FDFF=194.0534;
% ro_FDFF=659.3234;
% 
% 
% F_AO_FDFF=198.8190;
% eps_AO_FDFF=887.4548;
% ws_AO_FDFF=5.4169;
% wf_AO_FDFF=1684.9944;
% a1_AO_FDFF=4359.4855;
% ro_AO_FDFF=4451.7299;
% 
% %%%%%%%%%%%%%



% %%%%%%%%%%%%%%%%% Third-order
% % %% Custom parameters
%   % HDSTA Implicit
%     mu_i=   81.7318;
%     % HDSTA Explicit
%     mu_e=0.1932;
%     % Quadratic Explicit
%     F_e=3.8498;
%     alpha_e=0.4715;
%     % Quadratic Implicit
%     F_i=15.2169;
%     alpha_i=3.3914;
%     % Linear
%     c=16.9388;
%     % Toolbox
%     r=1.2621;
%     % ALIEN
%     T_alien=0.8497;
%     k_alien=2;
%     mu_alien=1;
%     %%%%%%%%%%%%%%%%%%%%
%     L_E_STD=0.7532;
%     L_I_STD=0.7728;
%     L_SI_STD=0.7281;
%     L_E_AO_STD=2.4167;
%     L_I_AO_STD= 4.0095;
%     L_SI_AO_STD=2.2759;
%     L_E_HD_STD=   1.1782;
%     L_I_HD_STD= 2.4025;
%     L_EHDD=1.9111;
%     L_EGHDD=1.9649;
%     L_IHDD=4.2086;
%     L_IGHDD=   3.3526;
% 
% %     
% %     mu_VGED=
% %     tau_VGED=
% %     wc_VGED=
% %     q_VGED=
% %     
% %     L_SI_HD_STD=
% %     mu_SI_HD_STD=
%     
% %     eps_E_STDAC=
% %     alpha_E_STDAC=
% 
%     R_kalman= 3.8107e-7;
% 
% 
% %     omegas_FDFF=
% %     omegaf_FDFF=
% %     gamma_FDFF=
% %     ro_FDFF=
% 
% 
%     F_AO_FDFF= 9.8148;
%     eps_AO_FDFF= 38.2115;
%     ws_AO_FDFF= 20.6008;
%     wf_AO_FDFF= 159.0598;
%     a1_AO_FDFF=129.0693;
%     ro_AO_FDFF=149.3592;
% 
% %%%%%%%%%%%%%



% 
% %%%%%%%%%%%%%%%% 
% % %% Custom parameters
%   % HDSTA Implicit
%     mu_i=   81.7318;
%     % HDSTA Explicit
%     mu_e=0.1932;
%     % Quadratic Explicit
%     F_e=3.8498;
%     alpha_e=0.4715;
%     % Quadratic Implicit
%     F_i=15.2169;
%     alpha_i=3.3914;
%     % Linear
%     c=16.9388;
%     % Toolbox
%     r=1.2621;
%     % ALIEN
%     T_alien=0.8497;
%     k_alien=2;
%     mu_alien=1;
%     %%%%%%%%%%%%%%%%%%%%
%     L_E_STD=0.7532;
%     L_I_STD=0.7728;
%     L_SI_STD=0.7281;
%     L_E_AO_STD=2.4167;
%     L_I_AO_STD= 4.0095;
%     L_SI_AO_STD=2.2759;
%     L_E_HD_STD=   1.1782;
%     L_I_HD_STD= 2.4025;
%     L_EHDD=1.9111;
%     L_EGHDD=1.9649;
%     L_IHDD=4.2086;
%     L_IGHDD=   3.3526;
% 
% %     
% %     mu_VGED=
% %     tau_VGED=
% %     wc_VGED=
% %     q_VGED=
% %     
% %     L_SI_HD_STD=
% %     mu_SI_HD_STD=
%     
% %     eps_E_STDAC=
% %     alpha_E_STDAC=
% 
%     R_kalman= 1;
% 
% 
%     omegas_FDFF=100
%     omegaf_FDFF=100
%     gamma_FDFF=100
%     ro_FDFF=100
% 
% 
%     F_AO_FDFF= 100;
%     eps_AO_FDFF= 100;
%     ws_AO_FDFF= 100;
%     wf_AO_FDFF= 500;
%     a1_AO_FDFF=500;
%     ro_AO_FDFF=500;

%%%%%%%%%%%%%













