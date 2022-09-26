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




% MATLAB script for comparison between differentiators

%% Initializations
initialization % Clear screen, remove variables and close plots

%% Parameters
parameters % Parameters of the differentiators and simulation

%% Outer-loop for calculating objective functions for each condition

for i=1:length(h_vector) % Different conditions

    h=h_vector(i); % Select a sampling period
    w=w_vector(i); % Select an input frequency
	len=floor(t_f/h); % length of each simulation
%     if h~=0.05
%         T=15*h;
%     end

    %% Calculate time vectors
    t=0:h:t_f-h; % Simulation time vector
    tc=0:hc:max(t); % Time vector for the interpolated vectors

    %% Pre allocation to increase the calculation speed
    u=zeros(24,len);
    e_interpolated=zeros(24,length(tc));
    
    %% Generate input signal and its differentiation
    [f,df,ddf,dddf]=signal_generator(len,w,amp,phi,SNR(i),h,noise_type,An_vector(i),Wn_vector(i),Quant);
    diff=[df;ddf;dddf]; % Vector of differentiation
	reset_flag=1; % To reset the internal variables of the differentiators
    
    %% Parameter selection for parameter tuning only
    if changing_parameter == "param" %&& i>10
    param_select
    end
    
	%% iner-loop for differentiation
	for k=1:floor(len)+1   % k is the time step
        
        [u(1,k),CT(1)] = Euler(h,f(k),reset_flag,diff_order,x0);   
        
        [u(2,k),CT(2)] = LF(f(k),c,h,reset_flag,diff_order,x0);
        
        [u(3,k),CT(3)] = E_STD(L_E_STD,h,f(k),reset_flag,diff_order,x0);
        
        [u(4,k),CT(4)] = I_STD(L_I_STD,h,f(k),reset_flag,diff_order,x0);
        
        [u(5,k),CT(5)] = SI_STD_sch1(L_SI_STD,h,f(k),reset_flag,diff_order,x0);
        
        [u(6,k),CT(6)] = E_HD_STD(L_E_HD_STD,mu_e,h,f(k),reset_flag,diff_order,x0);
        
        [u(7,k),CT(7)] = I_HD_STD(L_I_HD_STD,mu_i,h,f(k),reset_flag,diff_order,x0);
        
        [u(8,k),CT(8)] = E_QD(f(k),F_e,alpha_e,h,reset_flag,diff_order,x0);  

        [u(9,k),CT(9)] = I_QD(f(k),F_i,alpha_i,h,reset_flag,diff_order,x0);     
        
        [u(10,k),CT(10)] = alien(T_alien,k_alien,mu_alien,h,f(k),reset_flag,diff_order,x0);   
 
        [u(11,k),CT(11)] = toolbox_diff(f(k),r,h,diff_order,reset_flag,x0);
         
        [u(12,k),CT(12)] = E_AO_STD(f(k),h,L_E_AO_STD,reset_flag,diff_order,AO_order,x0);    
        
        [u(13,k),CT(13)] = I_AO_STD(f(k),h,L_I_AO_STD,reset_flag,diff_order,AO_order,x0);
        
        [u(14,k),CT(14)] = SI_AO_STD_sch1(f(k),h,L_SI_AO_STD,reset_flag,diff_order,AO_order,x0); 
        
        [u(15,k),CT(15)] = E_HDD(f(k),h,L_EHDD,reset_flag,diff_order,x0);
        
        [u(16,k),CT(16)] = E_GHDD(f(k),h,L_EGHDD,reset_flag,diff_order,x0);
        
        [u(17,k),CT(17)] = I_HDD(f(k),h,L_IHDD,reset_flag,diff_order,x0);
        
        [u(18,k),CT(18)] = I_GHDD(f(k),h,L_IGHDD,reset_flag,diff_order,x0);       
        
        [u(19,k),CT(19)] = VGED(mu_VGED,tau_VGED,wc_VGED,q_VGED,h,f(k),reset_flag,diff_order,x0);
        
        [u(20,k),CT(20)] = SI_HD_STD(L_SI_HD_STD,mu_SI_HD_STD,h,f(k),reset_flag,diff_order,x0);
        
        [u(21,k),CT(21)] = E_STDAC(alpha_E_STDAC,eps_E_STDAC,h,f(k),reset_flag,diff_order,x0);
        
        [u(22,k),CT(22)] = I_FDFF(omegas_FDFF,omegaf_FDFF,gamma_FDFF,ro_FDFF,f(k),h,reset_flag,diff_order,x0);
        
        [u(23,k),CT(23)] = I_AO_FDFF(F_AO_FDFF,eps_AO_FDFF,ws_AO_FDFF,wf_AO_FDFF,a1_AO_FDFF,ro_AO_FDFF,f(k),h,reset_flag,diff_order,x0);
        
        [u(24,k),CT(24)] = kalman_diff(R_kalman,h,f(k),reset_flag,diff_order,x0);

        reset_flag=0; % To avoid resetting the internal variables of the differentiators anymore
        
	end
    [M,N]=size(u);
    
	%% Calculate errors
    
    e=u-repelem(diff(diff_order,:),M,1);
    
    
	%% Interpolating the error
    
    for j=1:M
        e_interpolated(j,:)=interp1q(t',e(j,:)',tc')';
    end
    
	%% cost functions
    % L2 norm of the error
    if ssp==1
        LV=floor(len/2):len+1; % vector for calculating norms in steady-state
        for j=1:M
            L2(j,i)=norm(e(j,LV))/len ;           
        end
    % Interpolated L2 norm of the error    
        for j=1:M
            L2_interpolated(j,i)=norm(e_interpolated(j,LV))/len;
        end    
    
    % Variance of the error    
    
        for j=1:M
            VAR(j,i)=variation(u(j,LV));
        end  
    
    % THD of the error    
    
        for j=1:M
            THD(j,i)=thd_index(u(j,LV),h,t_f-min(LV)*h+h,w);
        end  

    else
        
        for j=1:M
            L2(j,i)=norm(e(j,:))/len;     
        end
        

    % Interpolated L2 norm of the error    
	for j=1:M
            L2_interpolated(j,i)=norm(e_interpolated(j,:))/len;
	end    
    
    % Variance of the error    
        for j=1:M
            VAR(j,i)=variation(u(j,:));
        end     
    
    % THD of the error    
        for j=1:M
            THD(j,i)=thd_index(u(j,:),h,t_f,w);
            if THD(j,i)==inf
                THD(j,i)=0;
            end
        end  
             
    end

    % L_inf norm of the error  
	for j=1:M
            L_inf(j,i)=norm(e(j,:),inf);
	end  
     
    
    %% Parameter tuning
	if changing_parameter == "param"

        % Calculating objective function for current parameters
        for j=1:M

            local(j)=[k1 k2 k3 k4 k5]*[L2(j,i)   L2_interpolated(j,i)   L_inf(j,i)   VAR(j,i)   THD(j,i)]';
        end
    
        
        %Euler
        if local(1)<Global(1)
            Global(1)=local(1);
        end
        
        %LF
        if local(2)<Global(2)
            Global(2)=local(2);            
            Global_c=c;
        end       
        
        %E_STD
        if local(3)<Global(3)
            Global(3)=local(3);  
            L_E_STD_global=L_E_STD;
        end   
        
        %I_STD
        if local(4)<Global(4)
            Global(4)=local(4);    
            L_I_STD_global=L_I_STD;            
        end     
        
        %SI_STD
        if local(5)<Global(5)
            Global(5)=local(5);    
            L_SI_STD_global=L_SI_STD;            
        end     
        
        %E-HD-STD
        if local(6)<Global(6)
            Global(6)=local(6);
            Global_mu_e=mu_e;
            L_E_HD_STD_global=L_E_HD_STD;        
        end
        
        %I-HD-STD        
        if local(7)<Global(7)
            Global(7)=local(7);            
            Global_mu_i=mu_i;
            L_I_HD_STD_global=L_I_HD_STD;            
        end

        %E-QD
        if local(8)<Global(8)
            Global(8)=local(8);            
            Global_F_e=F_e;
            Global_alpha_e=alpha_e;
        end
        
        %I-QD
        if local(9)<Global(9)
            Global(9)=local(9);            
            Global_F_i=F_i;
            Global_alpha_i=alpha_i;
        end
        

        %ALIEN
        if local(10)<Global(10)
            Global(10)=local(10);            
             Global_T_alien=T_alien;   
             Global_k_alien=k_alien;
             Global_mu_alien=mu_alien;
        end
   
        %HD
        if local(11)<Global(11)
            Global(11)=local(11); 
            Global_r=r;
        end
        
        %E-AO_STD
        if local(12)<Global(12)
            Global(12)=local(12);   
            L_E_AO_STD_global=L_E_AO_STD;            
        end
        
        %I-AO_STD
        if local(13)<Global(13)
            Global(13)=local(13);       
            L_I_AO_STD_global=L_I_AO_STD;                 
        end    
        
       %SI-AO_STD
        if local(14)<Global(14)
            Global(14)=local(14);  
            L_SI_AO_STD_global=L_SI_AO_STD;                 
        end   
        
        % EHDD
        if local(15)<Global(15)
            Global(15)=local(15);  
            L_EHDD_global=L_EHDD;                 
        end  
            
        % EGHDD
        if local(16)<Global(16)
            Global(16)=local(16);  
            L_EGHDD_global=L_EGHDD;                 
        end   
        
        % IHDD
        if local(17)<Global(17)
            Global(17)=local(17);  
            L_IHDD_global=L_IHDD;                 
        end  
            
        % IGHDD
        if local(18)<Global(18)
            Global(18)=local(18);  
            L_IGHDD_global=L_IGHDD;                 
        end  
        
        %
        if local(19)<Global(19) 
            Global(19)=local(19);
            mu_VGED_global=mu_VGED;
            tau_VGED_global=tau_VGED;
            wc_VGED_global=wc_VGED;
            q_VGED_global=q_VGED;
        end       
        
        if local(20)<Global(20)
            Global(20)=local(20);  
            L_SI_HD_STD_global=L_SI_HD_STD;
            mu_SI_HD_STD_global=mu_SI_HD_STD;         
        end   
        
        if local(21)<Global(21)
            Global(21)=local(21);  
            alpha_E_STDAC_global=alpha_E_STDAC;
            eps_E_STDAC_global=eps_E_STDAC;         
        end    
        
        if local(22)<Global(22)
            Global(22)=local(22);  
            omegas_FDFF_global=omegas_FDFF;
            omegaf_FDFF_global=omegaf_FDFF;  
            ro_FDFF_global=ro_FDFF;
            gamma_FDFF_global=gamma_FDFF;
        end          
        
        if local(23)<Global(23)
            Global(23)=local(23);  
            F_AO_FDFF_global=F_AO_FDFF;
            eps_AO_FDFF_global=eps_AO_FDFF;
            ws_AO_FDFF_global=ws_AO_FDFF;
            wf_AO_FDFF_global=wf_AO_FDFF;
            a1_AO_FDFF_global=a1_AO_FDFF;
            ro_AO_FDFF_global=ro_AO_FDFF;
        end    
        
        if local(24)<Global(24)
            Global(24)=local(24);  
            R_kalman_global=R_kalman;
        end           
        
	else
        
        for j=1:M
            Global(j,i)=[k1 k2 k3 k4 k5]*[L2(j,i)   L2_interpolated(j,i)   L_inf(j,i)   VAR(j,i)   THD(j,i)]';
        end
        
	end
    
    
    if mod(i,100)==0 && changing_parameter=="param"
        clc
        fprintf('-------  Optimization Progress= %.0f %%  ------- \n', 100*i/sample_length)
    end
    
    if sample_length == 1  

        % Latex friendly:
        fprintf('Euler  & %.4f   & %.4f     & %.4f  & %.4f &   %.4f & %.2f $\\beta$ \\\\ \n', L2(1,:),L2_interpolated(1,:),L_inf(1,:),VAR(1,:),THD(1,:),CT(1)/min(CT))
        fprintf('LF  & %.4f   & %.4f    & %.4f  & %.4f &  %.4f & %.2f  $\\beta$ \\\\ \n', L2(2,:),L2_interpolated(2,:),L_inf(2,:),VAR(2,:),THD(2,:),CT(2)/min(CT))
        fprintf('E-STD  & %.4f   & %.4f     & %.4f & %.4f  &  %.4f & %.2f $\\beta$  \\\\ \n', L2(3,:),L2_interpolated(3,:),L_inf(3,:),VAR(3,:),THD(3,:),CT(3)/min(CT))
        fprintf('I-STD  & %.4f   & %.4f     & %.4f & %.4f  &  %.4f & %.2f $\\beta$  \\\\ \n', L2(4,:),L2_interpolated(4,:),L_inf(4,:),VAR(4,:),THD(4,:),CT(4)/min(CT))
        fprintf('SI-STD  & %.4f   & %.4f     & %.4f & %.4f  &  %.4f & %.2f $\\beta$  \\\\ \n', L2(5,:),L2_interpolated(5,:),L_inf(5,:),VAR(5,:),THD(5,:),CT(5)/min(CT))
        fprintf('E-HD-STD   & %.4f     & %.4f  & %.4f &  %.4f & %.4f & %.2f $\\beta$ \\\\ \n', L2(6,:),L2_interpolated(6,:),L_inf(6,:),VAR(6,:),THD(6,:),CT(6)/min(CT))
        fprintf('I-HD-STD   & %.4f     & %.4f  & %.4f &  %.4f & %.4f & %.2f $\\beta$ \\\\ \n', L2(7,:),L2_interpolated(7,:),L_inf(7,:),VAR(7,:),THD(7,:),CT(7)/min(CT))
        fprintf('E-QD  & %.4f   & %.4f     & %.4f & %.4f  &  %.4f & %.2f $\\beta$  \\\\ \n', L2(8,:),L2_interpolated(8,:),L_inf(8,:),VAR(8,:),THD(8,:),CT(8)/min(CT))
        fprintf('I-QD  & %.4f   & %.4f    & %.4f & %.4f  &  %.4f & %.2f  $\\beta$ \\\\ \n', L2(9,:),L2_interpolated(9,:),L_inf(9,:),VAR(9,:),THD(9,:),CT(9)/min(CT))
        fprintf('ALIEN  & %.4f   & %.4f     & %.4f  &  %.4f & %.4f & %.2f $\\beta$  \\\\ \n', L2(10,:),L2_interpolated(10,:),L_inf(10,:),VAR(10,:),THD(10,:),CT(10)/min(CT))
        fprintf('HD  & %.4f   & %.4f     & %.4f  &  %.4f & %.4f & %.2f  $\\beta$ \\\\ \n', L2(11,:),L2_interpolated(11,:),L_inf(11,:),VAR(11,:),THD(11,:),CT(11)/min(CT))
        fprintf('E-AO-STD  & %.4f     & %.4f   & %.4f & %.4f  &  %.4f & %.2f $\\beta$  \\\\ \n', L2(12,:),L2_interpolated(12,:),L_inf(12,:),VAR(12,:),THD(12,:),CT(12)/min(CT))
        fprintf('I-AO-STD  & %.4f     & %.4f   & %.4f & %.4f  &  %.4f & %.2f $\\beta$  \\\\ \n', L2(13,:),L2_interpolated(13,:),L_inf(13,:),VAR(13,:),THD(13,:),CT(13)/min(CT))
        fprintf('SI-AO-STD  & %.4f     & %.4f   & %.4f & %.4f  &  %.4f & %.2f  $\\beta$ \\\\ \n', L2(14,:),L2_interpolated(14,:),L_inf(14,:),VAR(14,:),THD(14,:),CT(14)/min(CT))
        fprintf('E-HDD  & %.4f     & %.4f   & %.4f & %.4f  &  %.4f & %.2f  $\\beta$ \\\\ \n', L2(15,:),L2_interpolated(15,:),L_inf(15,:),VAR(15,:),THD(15,:),CT(15)/min(CT))
        fprintf('E-GHDD  & %.4f     & %.4f   & %.4f & %.4f  &  %.4f & %.2f  $\\beta$ \\\\ \n', L2(16,:),L2_interpolated(16,:),L_inf(16,:),VAR(16,:),THD(16,:),CT(16)/min(CT))
        fprintf('I-HDD  & %.4f     & %.4f   & %.4f & %.4f  &  %.4f & %.2f  $\\beta$ \\\\ \n', L2(17,:),L2_interpolated(17,:),L_inf(17,:),VAR(17,:),THD(17,:),CT(17)/min(CT))
        fprintf('I-GHDD  & %.4f     & %.4f   & %.4f & %.4f  &  %.4f & %.2f  $\\beta$ \\\\ \n', L2(18,:),L2_interpolated(18,:),L_inf(18,:),VAR(18,:),THD(18,:),CT(18)/min(CT))
        fprintf('VGED  & %.4f     & %.4f   & %.4f & %.4f  &  %.4f & %.2f  $\\beta$ \\\\ \n', L2(19,:),L2_interpolated(19,:),L_inf(19,:),VAR(19,:),THD(19,:),CT(19)/min(CT))
        fprintf('SI-HD-STD  & %.4f     & %.4f   & %.4f & %.4f  &  %.4f & %.2f  $\\beta$ \\\\ \n', L2(20,:),L2_interpolated(20,:),L_inf(20,:),VAR(20,:),THD(20,:),CT(20)/min(CT))
        fprintf('E-STDAC  & %.4f     & %.4f   & %.4f & %.4f  &  %.4f & %.2f  $\\beta$ \\\\ \n', L2(21,:),L2_interpolated(21,:),L_inf(21,:),VAR(21,:),THD(21,:),CT(21)/min(CT))
        fprintf('I-FDFF  & %.4f     & %.4f   & %.4f & %.4f  &  %.4f & %.2f  $\\beta$ \\\\ \n', L2(22,:),L2_interpolated(22,:),L_inf(22,:),VAR(22,:),THD(22,:),CT(22)/min(CT))
        fprintf('I-AO-FDFF  & %.4f     & %.4f   & %.4f & %.4f  &  %.4f & %.2f  $\\beta$ \\\\ \n', L2(23,:),L2_interpolated(23,:),L_inf(23,:),VAR(23,:),THD(23,:),CT(23)/min(CT))        
        fprintf('Kalman  & %.4f     & %.4f   & %.4f & %.4f  &  %.4f & %.2f  $\\beta$ \\\\ \n', L2(24,:),L2_interpolated(24,:),L_inf(24,:),VAR(24,:),THD(24,:),CT(24)/min(CT))        
        betta=1000*min(CT)
    end
    
    
    
       
end


%% Printing the optimal parameters and saving them to a file
if changing_parameter == "param"

    mu_e=Global_mu_e; mu_i=Global_mu_i;

    F_e=Global_F_e; alpha_e=Global_alpha_e;
    
    F_i=Global_F_i; alpha_i=Global_alpha_i;

    c=Global_c;
    
    r=Global_r;

    k_alien=Global_k_alien; mu_alien=Global_mu_alien; T_alien=Global_T_alien;
    
    L_E_STD=L_E_STD_global;
    L_I_STD=L_I_STD_global;
    L_SI_STD=L_SI_STD_global;
    L_E_AO_STD=L_E_AO_STD_global;
    L_I_AO_STD=L_I_AO_STD_global;  
    L_SI_AO_STD=L_SI_AO_STD_global;
    L_E_HD_STD=L_E_HD_STD_global; 
    L_I_HD_STD= L_I_HD_STD_global; 
    
    mu_VGED=mu_VGED_global;
    tau_VGED=tau_VGED_global;
    wc_VGED=wc_VGED_global;
    q_VGED=q_VGED_global;

    L_EHDD=L_EHDD_global;
    L_EGHDD=L_EGHDD_global;
    
    L_IHDD=L_IHDD_global;
    L_IGHDD=L_IGHDD_global; 
    
    L_SI_HD_STD=L_SI_HD_STD_global;
    mu_SI_HD_STD=mu_SI_HD_STD_global;  
    
    alpha_E_STDAC=alpha_E_STDAC_global;
    eps_E_STDAC=eps_E_STDAC_global;      
    
    omegas_FDFF=omegas_FDFF_global;
    omegaf_FDFF=omegaf_FDFF_global;
    ro_FDFF=ro_FDFF_global;
    gamma_FDFF=gamma_FDFF_global;
    
    
    F_AO_FDFF=F_AO_FDFF_global;
    eps_AO_FDFF=eps_AO_FDFF_global;
    ws_AO_FDFF=ws_AO_FDFF_global;
    wf_AO_FDFF=wf_AO_FDFF_global;
    a1_AO_FDFF=a1_AO_FDFF_global;
    ro_AO_FDFF=ro_AO_FDFF_global;
    
    R_kalman=R_kalman_global;

    save parameters_file F_e r alpha_e F_i alpha_i c  T_alien k_alien mu_alien mu_e mu_i L_E_STD L_I_STD  L_SI_STD L_E_AO_STD L_I_AO_STD L_SI_AO_STD L_E_HD_STD...
        L_I_HD_STD L_EHDD L_EGHDD L_IHDD L_IGHDD mu_VGED tau_VGED wc_VGED q_VGED L_SI_HD_STD mu_SI_HD_STD alpha_E_STDAC eps_E_STDAC...
        omegas_FDFF omegaf_FDFF ro_FDFF gamma_FDFF F_AO_FDFF eps_AO_FDFF ws_AO_FDFF wf_AO_FDFF a1_AO_FDFF ro_AO_FDFF R_kalman

    fprintf('--------------------------------------------------- \n')
    fprintf('Euler \\cref{E_Euler} & No parameter & %.4f  \\\\ \n', Global(1))
    fprintf('LF \\cref{E_linear} & $c$=%.4f & %.4f  \\\\ \n', c,Global(2))
    fprintf('E-STD \\cref{E_STA_explicit} & $L$=%.4f & %.4f \\\\ \n', L_E_STD,Global(3))
    fprintf('I-STD (\\cref{F_I_STD_flowchar}) & $L$=%.4f & %.4f  \\\\  \n', L_I_STD,Global(4))
    fprintf('SI-STD \\cref{E_STD_1_semi} & $L$=%.4f & %.4f  \\\\  \n', L_SI_STD,Global(5))
    fprintf('E-HD-STD \\cref{E_HDSTA_explicit} &  $L$=%.4f, $\\mu$=%.4f & %.4f  \\\\  \n', L_E_HD_STD,mu_e,Global(6))
    fprintf('I-HD-STD (\\cref{F_I_HD_STD_flowchar}) &  $L$=%.4f, $\\mu$=%.4f & %.4f  \\\\  \n', L_I_HD_STD,mu_i,Global(7))
    fprintf('E-QD \\cref{E_quadratic_explicit} & $F$=%.4f, $\\alpha$= %.4f  &  %.4f  \\\\   \n', F_e,alpha_e,Global(8))
    fprintf('I-QD \\cref{E_quadratic_implicit_algorithm} & $F$=%.4f, $\\alpha$=%.4f & %.4f  \\\\   \n', F_i,alpha_i,Global(9))
    fprintf('ALIEN \\cref{E_ALIEN} & $T$=%.4f , $\\kappa$=%.0f , $\\mu$=%.0f & %.4f  \\\\   \n', T_alien,k_alien,mu_alien,Global(10))
    fprintf('HD \\cref{E_system_obs_HD_disc}$^*$ & $r$=%.4f & %.4f  \\\\   \n', r,Global(11))
    fprintf('E-AO-STD \\cref{E_STA_AO_explicit}$^{**}$ & $L$=%.4f  & %.4f  \\\\   \n',L_E_AO_STD,Global(12))
    fprintf('I-AO-STD (\\cref{F_I_AO_STD_flowchar}) $^{**}$& $L$=%.4f  & %.4f  \\\\   \n', L_I_AO_STD,Global(13))
    fprintf('SI-AO-STD \\cref{E_AO_STD_1_semi}$^{**}$& $L$=%.4f  & %.4f  \\\\   \n', L_SI_AO_STD,Global(14))
    fprintf('E-HDD \\cref{E_HDD_E}$^{**}$ & $L$=%.4f  & %.4f  \\\\   \n', L_EHDD,Global(15))
    fprintf('E-GHDD \\cref{E_new_obs_non}$^{**}$ & $L$=%.4f  & %.4f  \\\\   \n', L_EGHDD,Global(16))
    fprintf('I-HDD (\\cref{F_flowchart_I_HDD})$^{**}$ & $L$=%.4f  & %.4f  \\\\   \n', L_IHDD,Global(17))
    fprintf('I-GHDD (\\cref{F_flowchart_I_GHDD})$^{**}$ & $L$=%.4f  & %.4f  \\\\   \n', L_IGHDD,Global(18))    
    fprintf('VGED \\cref{E_VGED_explicit} & $\\mu$=%.4f, $\\tau$=%.4f, $\\omega_c$=%.4f, $q$=%.4f  & %.4f  \\\\   \n', mu_VGED,tau_VGED,wc_VGED,q_VGED,Global(19))
    fprintf('SI-HD-STD \\cref{E_semi_HDSTD} &  $L$=%.4f, $\\mu$=%.4f & %.4f  \\\\  \n', L_SI_HD_STD,mu_SI_HD_STD,Global(20))
    fprintf('E-STDAC \\cref{E_adaptive_reich} &  $\\alpha$=%.4f, $\\epsilon$=%.4f & %.4f  \\\\  \n', alpha_E_STDAC,eps_E_STDAC,Global(21))
    fprintf('I-FDFF (\\Cref{F_SHMD_flowchar}) &  $\\omega_s$=%.4f, $\\omega_f$=%.4f,  $\\rho$=%.4f, $\\gamma$=%.4f & %.4f  \\\\  \n', omegas_FDFF,omegaf_FDFF,ro_FDFF,gamma_FDFF,Global(22))
    fprintf('I-AO-FDFF (\\Cref{F_SHMD_flowchar}) &  $F$=%.4f, $\\epsilon$=%.4f  $\\omega_s$=%.4f, $\\omega_f$=%.4f, $a_1$=%.4f, $\\rho$=%.4f & %.4f  \\\\  \n', F_AO_FDFF,eps_AO_FDFF,ws_AO_FDFF,wf_AO_FDFF,a1_AO_FDFF,ro_AO_FDFF,Global(23))
    fprintf('Kalman (\\Cref{S_kalman}) &  $R$=%.4f & %.4f  \\\\  \n', R_kalman,Global(24))
end


%Showing the plots
if changing_parameter~="param"
        plotting(changing_parameter,plot_type,h_vector,w_vector,L2,L2_interpolated,L_inf,VAR,THD,Global,e,u,diff,t,SNR,diff_order,sample_length,Wn_vector,An_vector)
end

% N_max_arbitrary
% N_total_arbitrary

% load data_Levant
% load data_Levant_same
% 
% L_inf_arbitrary_e_3=L_inf(12,:);
% L_inf_arbitrary_i_3=L_inf(13,:);
% 
% % % % 
% % save data_Levant   L_inf_arbitrary_e_2 L_inf_arbitrary_i_2 %   L_inf_arbitrary_e_3 L_inf_arbitrary_i_3   L_inf_arbitrary_e_4 L_inf_arbitrary_i_4    L_inf_arbitrary_e_5 L_inf_arbitrary_i_5    L_inf_arbitrary_e_6 L_inf_arbitrary_i_6    L_inf_arbitrary_e_7 L_inf_arbitrary_i_7 h_vector
% save data_Levant_same   L_inf_arbitrary_e_1 L_inf_arbitrary_i_1    L_inf_arbitrary_e_2 L_inf_arbitrary_i_2   L_inf_arbitrary_e_3 L_inf_arbitrary_i_3 %   L_inf_arbitrary_e_5 L_inf_arbitrary_i_5    L_inf_arbitrary_e_6 L_inf_arbitrary_i_6    L_inf_arbitrary_e_7 L_inf_arbitrary_i_7 h_vector
% % save data_Levant L_inf_STA_e  L_inf_STA_i  L_inf_arbitrary_e_2 L_inf_arbitrary_i_2 L_inf_arbitrary_e_3 L_inf_arbitrary_i_3 L_inf_arbitrary_e_4 L_inf_arbitrary_i_4 L_inf_arbitrary_e_5 L_inf_arbitrary_i_5 L_inf_arbitrary_e_6 L_inf_arbitrary_i_6 L_inf_arbitrary_e_7 L_inf_arbitrary_i_7 h_vector
