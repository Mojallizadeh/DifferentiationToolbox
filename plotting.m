% This file is part of the differentiation toolbox
% developed by Mohammad Rasool Mojallizadeh at INRIA, France
% Email: mohammad-rasool.mojallizadeh@inria.fr
%
% If you use this toolbox please cite the following report:
%
% Mohammad Rasool Mojallizadeh, Bernard Brogliato, Vincent Acary. Discrete-time differentiators:
% design and comparative analysis. 2020. hal-02960923
%
% This toolbox allows studying the behavior of several types of differentiator
% under different realistic conditions.
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

function [] = plotting(changing_parameter,plot_type,h_vector,w_vector,L2,L2_interpolated,L_inf,VAR,THD,Global,e,u,diff,t,SNR,diff_order,sample_length,Wn_vector,An_vector)

%% Plot all the information

    if sample_length > 1

        if changing_parameter=='h'
            samp = h_vector;
            label="Sampling period ($h$)";
        elseif changing_parameter=="w"
            samp = w_vector;
            label="$\omega$ (rad/s)";
        elseif changing_parameter=="w/h"
            samp = (w_vector)./(h_vector);
            label="$\omega/h$ (rad/s)";
        elseif changing_parameter=="h/w"
            samp = (h_vector)./(w_vector);
            label="$h/\omega$ (1/rad)";
        elseif changing_parameter=="SNR"
            samp = SNR;
            label="SNR";
        elseif changing_parameter=="An" 
            samp=An_vector;
            label="Noise amplitude";
        elseif changing_parameter=="Wn" 
            samp=Wn_vector;
            label="Noise frequency (rad/s)";            
        else
            disp('Error')
            return
        end

        sgtitle("$\bar{L}_2$",'Interpreter','latex','FontSize',18) 
        subplot(2,2,1)
        ploter(samp,L2(3,:),plot_type,'k')
        hold on
        ploter(samp,L2(6,:),plot_type,'r')
        hold on   
        ploter(samp,L2(8,:),plot_type,'g')
        hold on        
        ploter(samp,L2(12,:),plot_type,'b')  
        hold on
        ploter(samp,L2(15,:),plot_type,'c')   
        hold on
        ploter(samp,L2(16,:),plot_type,'m')             
        legend("E-STD","E-HD-STD","E-QD","E-AO-STD","E-HDD","E-GHDD",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(A)','FontSize',24)         
        
        subplot(2,2,2)
        ploter(samp,L2(4,:),plot_type,'k')
        hold on
        ploter(samp,L2(7,:),plot_type,'r')
        hold on     
        ploter(samp,L2(9,:),plot_type,'g')
        hold on        
        ploter(samp,L2(13,:),plot_type,'b')  
        hold on     
        ploter(samp,L2(17,:),plot_type,'c')
        hold on        
        ploter(samp,L2(18,:),plot_type,'m')   
        hold on
        ploter(samp,L2(22,:),plot_type,[0.8500, 0.3250, 0.0980])  
        hold on
        ploter(samp,L2(23,:),plot_type,[0.4940, 0.1840, 0.5560])                
        legend("I-STD","I-HD-STD","I-QD","I-AO-STD","I-HDD","I-GHDD","I-FDFF","I-AO-FDFF",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(B)','FontSize',24)         
        
        
        subplot(2,2,3)
        ploter(samp,L2(1,:),plot_type,'r')
        hold on
        ploter(samp,L2(2,:),plot_type,'g')
        hold on
        ploter(samp,L2(11,:),plot_type,'b')  
        hold on
        ploter(samp,L2(19,:),plot_type,'c') 
        hold on
        ploter(samp,L2(21,:),plot_type,'m')          
        legend("Euler","LF","HD","VGED","E-STDAC",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(C)','FontSize',24)         
        
        subplot(2,2,4)
        ploter(samp,L2(10,:),plot_type,'r')  
        hold on
        ploter(samp,L2(5,:),plot_type,'g')
        hold on    
        ploter(samp,L2(14,:),plot_type,'b')  
        hold on    
        ploter(samp,L2(20,:),plot_type,'c')
        hold on    
        ploter(samp,L2(24,:),plot_type,'m')          
        legend("ALIEN","SI-STD","SI-AO-STD","SI-HD-STD","Kalman",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(D)','FontSize',24)         


        figure
        
        sgtitle("$\tilde{L}_2$",'Interpreter','latex','FontSize',18) 
        subplot(2,2,1)
        ploter(samp,L2_interpolated(3,:),plot_type,'k')
        hold on    
        ploter(samp,L2_interpolated(6,:),plot_type,'r')
        hold on   
        ploter(samp,L2_interpolated(8,:),plot_type,'g')
        hold on        
        ploter(samp,L2_interpolated(12,:),plot_type,'b') 
        hold on
        ploter(samp,L2_interpolated(15,:),plot_type,'c')   
        hold on
        ploter(samp,L2_interpolated(16,:),plot_type,'m')          
        legend("E-STD","E-HD-STD","E-QD","E-AO-STD","E-HDD","E-GHDD",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(A)','FontSize',24)         
        
        subplot(2,2,2)
        ploter(samp,L2_interpolated(4,:),plot_type,'k')
        hold on
        ploter(samp,L2_interpolated(7,:),plot_type,'r')
        hold on     
        ploter(samp,L2_interpolated(9,:),plot_type,'g')
        hold on        
        ploter(samp,L2_interpolated(13,:),plot_type,'b') 
        hold on     
        ploter(samp,L2_interpolated(17,:),plot_type,'c')
        hold on        
        ploter(samp,L2_interpolated(18,:),plot_type,'m') 
        hold on
        ploter(samp,L2_interpolated(22,:),plot_type,[0.8500, 0.3250, 0.0980])  
        hold on
        ploter(samp,L2_interpolated(23,:),plot_type,[0.4940, 0.1840, 0.5560])          
        legend("I-STD","I-HD-STD","I-QD","I-AO-STD","I-HDD","I-GHDD","I-FDFF","I-AO-FDFF",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(B)','FontSize',24)         
        
        
        subplot(2,2,3)
        ploter(samp,L2_interpolated(1,:),plot_type,'r')
        hold on
        ploter(samp,L2_interpolated(2,:),plot_type,'g')
        hold on
        ploter(samp,L2_interpolated(11,:),plot_type,'b')
        hold on
        ploter(samp,L2_interpolated(19,:),plot_type,'c') 
        hold on
        ploter(samp,L2_interpolated(21,:),plot_type,'m')           
        legend("Euler","LF","HD","VGED","E-STDAC",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(C)','FontSize',24)         
        
        subplot(2,2,4)
        ploter(samp,L2_interpolated(10,:),plot_type,'r')  
        hold on
        ploter(samp,L2_interpolated(5,:),plot_type,'g')
        hold on          
        ploter(samp,L2_interpolated(14,:),plot_type,'b')
        hold on          
        ploter(samp,L2_interpolated(20,:),plot_type,'c')
        hold on          
        ploter(samp,L2_interpolated(24,:),plot_type,'m')            
        legend("ALIEN","SI-STD","SI-AO-STD","SI-HD-STD","Kalman",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(D)','FontSize',24)         

        figure
        
        sgtitle("$L_\infty$",'Interpreter','latex','FontSize',18) 
        subplot(2,2,1)
        ploter(samp,L_inf(3,:),plot_type,'k')
        hold on    
        ploter(samp,L_inf(6,:),plot_type,'r')
        hold on   
        ploter(samp,L_inf(8,:),plot_type,'g')
        hold on        
        ploter(samp,L_inf(12,:),plot_type,'b')  
        hold on
        ploter(samp,L_inf(15,:),plot_type,'c') 
        hold on
        ploter(samp,L_inf(16,:),plot_type,'m')        
        legend("E-STD","E-HD-STD","E-QD","E-AO-STD","E-HDD","E-GHDD",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(A)','FontSize',24)         
        
        
        subplot(2,2,2)
        ploter(samp,L_inf(4,:),plot_type,'k')
        hold on
        ploter(samp,L_inf(7,:),plot_type,'r')
        hold on     
        ploter(samp,L_inf(9,:),plot_type,'g')
        hold on        
        ploter(samp,L_inf(13,:),plot_type,'b')  
        hold on     
        ploter(samp,L_inf(17,:),plot_type,'k')
        hold on        
        ploter(samp,L_inf(18,:),plot_type,'m')    
        hold on
        ploter(samp,L_inf(22,:),plot_type,[0.8500, 0.3250, 0.0980])  
        hold on
        ploter(samp,L_inf(23,:),plot_type,[0.4940, 0.1840, 0.5560])           
        legend("I-STD","I-HD-STD","I-QD","I-AO-STD","I-HDD","I-GHDD","I-FDFF","I-AO-FDFF",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(B)','FontSize',24)         
        
        
        subplot(2,2,3)
        ploter(samp,L_inf(1,:),plot_type,'r')
        hold on
        ploter(samp,L_inf(2,:),plot_type,'g')
        hold on
        ploter(samp,L_inf(11,:),plot_type,'b') 
        hold on
        ploter(samp,L_inf(19,:),plot_type,'c')
        hold on
        ploter(samp,L_inf(21,:),plot_type,'m')          
        legend("Euler","LF","HD","VGED","E-STDAC",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(C)','FontSize',24)         
        
        subplot(2,2,4)
        ploter(samp,L_inf(10,:),plot_type,'r')  
        hold on
        ploter(samp,L_inf(5,:),plot_type,'g')
        hold on           
        ploter(samp,L_inf(14,:),plot_type,'b') 
        hold on           
        ploter(samp,L_inf(20,:),plot_type,'c') 
        hold on           
        ploter(samp,L_inf(24,:),plot_type,'m')             
        legend("ALIEN","SI-STD","SI-AO-STD","SI-HD-STD","Kalman",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(D)','FontSize',24)         
        
        
       figure


        sgtitle("Var",'Interpreter','latex','FontSize',18) 
        subplot(2,2,1)
        ploter(samp,VAR(3,:),plot_type,'k')
        hold on
        ploter(samp,VAR(6,:),plot_type,'r')
        hold on   
        ploter(samp,VAR(8,:),plot_type,'g')
        hold on        
        ploter(samp,VAR(12,:),plot_type,'b')  
        hold on
        ploter(samp,VAR(15,:),plot_type,'c')    
        hold on
        ploter(samp,VAR(16,:),plot_type,'m')          
        legend("E-STD","E-HD-STD","E-QD","E-AO-STD","E-HDD","E-GHDD",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(A)','FontSize',24)         
        
        subplot(2,2,2)
        ploter(samp,VAR(4,:),plot_type,'k')
        hold on
        ploter(samp,VAR(7,:),plot_type,'r')
        hold on     
        ploter(samp,VAR(9,:),plot_type,'g')
        hold on        
        ploter(samp,VAR(13,:),plot_type,'b') 
        hold on     
        ploter(samp,VAR(17,:),plot_type,'c')
        hold on        
        ploter(samp,VAR(18,:),plot_type,'m')   
        hold on
        ploter(samp,VAR(22,:),plot_type,[0.8500, 0.3250, 0.0980])   
        hold on
        ploter(samp,VAR(23,:),plot_type,[0.4940, 0.1840, 0.5560])             
        legend("I-STD","I-HD-STD","I-QD","I-AO-STD","I-HDD","I-GHDD","I-FDFF","I-AO-FDFF",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(B)','FontSize',24)         
        
        
        subplot(2,2,3)
        ploter(samp,VAR(1,:),plot_type,'r')
        hold on
        ploter(samp,VAR(2,:),plot_type,'g')
        hold on
        ploter(samp,VAR(11,:),plot_type,'b')  
        hold on
        ploter(samp,VAR(19,:),plot_type,'c') 
        hold on
        ploter(samp,VAR(21,:),plot_type,'m')         
        legend("Euler","LF","HD","VGED","E-STDAC",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(C)','FontSize',24)         
        
        subplot(2,2,4)
        ploter(samp,VAR(10,:),plot_type,'r') 
        hold on         
        ploter(samp,VAR(5,:),plot_type,'g')
        hold on        
        ploter(samp,VAR(14,:),plot_type,'b')
        hold on        
        ploter(samp,VAR(20,:),plot_type,'c')
        hold on        
        ploter(samp,VAR(24,:),plot_type,'m')         
        legend("ALIEN","SI-STD","SI-AO-STD","SI-HD-STD","Kalman",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(D)','FontSize',24)         
        
        figure

       
        sgtitle("THD",'Interpreter','latex','FontSize',18) 
        subplot(2,2,1)
        ploter(samp,THD(3,:),plot_type,'k')
        hold on    
        ploter(samp,THD(6,:),plot_type,'r')
        hold on   
        ploter(samp,THD(8,:),plot_type,'g')
        hold on        
        ploter(samp,THD(12,:),plot_type,'b')   
        hold on
        ploter(samp,THD(15,:),plot_type,'c')   
        hold on
        ploter(samp,THD(16,:),plot_type,'m')        
        legend("E-STD","E-HD-STD","E-QD","E-AO-STD","E-HDD","E-GHDD",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)        
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(A)','FontSize',24)         
        
        subplot(2,2,2)
        ploter(samp,THD(4,:),plot_type,'k')
        hold on
        ploter(samp,THD(7,:),plot_type,'r')
        hold on     
        ploter(samp,THD(9,:),plot_type,'g')
        hold on        
        ploter(samp,THD(13,:),plot_type,'b')
        hold on     
        ploter(samp,THD(17,:),plot_type,'c')
        hold on        
        ploter(samp,THD(18,:),plot_type,'m')
        hold on
        ploter(samp,THD(22,:),plot_type,[0.8500, 0.3250, 0.0980])  
        hold on
        ploter(samp,THD(23,:),plot_type,[0.4940, 0.1840, 0.5560])           
        legend("I-STD","I-HD-STD","I-QD","I-AO-STD","I-HDD","I-GHDD","I-FDFF","I-AO-FDFF",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)        
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(B)','FontSize',24)         
        
        
        subplot(2,2,3)
        ploter(samp,THD(1,:),plot_type,'r')
        hold on
        ploter(samp,THD(2,:),plot_type,'g')
        hold on
        ploter(samp,THD(11,:),plot_type,'b')
        hold on
        ploter(samp,THD(19,:),plot_type,'c') 
        hold on
        ploter(samp,THD(21,:),plot_type,'m')            
        legend("Euler","LF","HD","VGED","E-STDAC",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)        
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(C)','FontSize',24)         
        
        subplot(2,2,4)
        ploter(samp,THD(10,:),plot_type,'r')
        hold on
        ploter(samp,THD(5,:),plot_type,'g')
        hold on       
        ploter(samp,THD(14,:),plot_type,'b')  
        hold on       
        ploter(samp,THD(20,:),plot_type,'c') 
        hold on       
        ploter(samp,THD(24,:),plot_type,'m')           
        legend("ALIEN","SI-STD","SI-AO-STD","SI-HD-STD","Kalman",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)        
        set(gca,'FontSize',18)    
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(D)','FontSize',24)         

        %--------------------------------
        

        figure

        sgtitle("$J$",'Interpreter','latex','FontSize',18) 
        subplot(2,2,1)
        ploter(samp,Global(3,:),plot_type,'k')
        hold on  
        ploter(samp,Global(6,:),plot_type,'r')
        hold on   
        ploter(samp,Global(8,:),plot_type,'g')
        hold on        
        ploter(samp,Global(12,:),plot_type,'b')  
        hold on
        ploter(samp,Global(15,:),plot_type,'c')  
        hold on
        ploter(samp,Global(16,:),plot_type,'m')          
        legend("E-STD","E-HD-STD","E-QD","E-AO-STD","E-HDD","E-GHDD",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)        
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(A)','FontSize',24)         
        
        subplot(2,2,2)
        ploter(samp,Global(4,:),plot_type,'k')
        hold on
        ploter(samp,Global(7,:),plot_type,'r')
        hold on     
        ploter(samp,Global(9,:),plot_type,'g')
        hold on        
        ploter(samp,Global(13,:),plot_type,'b') 
        hold on     
        ploter(samp,Global(17,:),plot_type,'c')
        hold on        
        ploter(samp,Global(18,:),plot_type,'m')  
        hold on
        ploter(samp,Global(22,:),plot_type,[0.8500, 0.3250, 0.0980])   
        hold on
        ploter(samp,Global(23,:),plot_type,[0.4940, 0.1840, 0.5560])            
        legend("I-STD","I-HD-STD","I-QD","I-AO-STD","I-HDD","I-GHDD","I-FDFF","I-AO-FDFF",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)        
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(B)','FontSize',24)         
        
        
        subplot(2,2,3)
        ploter(samp,Global(1,:),plot_type,'r')
        hold on
        ploter(samp,Global(2,:),plot_type,'g')
        hold on
        ploter(samp,Global(11,:),plot_type,'b') 
        hold on
        ploter(samp,Global(19,:),plot_type,'c') 
        hold on
        ploter(samp,Global(21,:),plot_type,'m')           
        legend("Euler","LF","HD","VGED","E-STDAC",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)        
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(C)','FontSize',24)         
        
        subplot(2,2,4)
        ploter(samp,Global(10,:),plot_type,'r')  
        hold on
        ploter(samp,Global(5,:),plot_type,'g')
        hold on             
        ploter(samp,Global(14,:),plot_type,'b')
        hold on             
        ploter(samp,Global(20,:),plot_type,'c')         
        hold on             
        ploter(samp,Global(24,:),plot_type,'m')         
        legend("ALIEN","SI-STD","SI-AO-STD","SI-HD-STD","Kalman",'FontSize',18)
        grid on
        xlabel(label,'Interpreter','latex','FontSize',18)        
        set(gca,'FontSize',18)  
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(D)','FontSize',24)         
        
        
        %--------------------------------
        
        
        
    else
        t(end+1)=t(end);
        sgtitle("Differentiation error",'FontSize',18)      
        subplot(2,2,1)
        ploter(t,e(3,:),plot_type,'k')
        hold on
        ploter(t,e(6,:),plot_type,'r')
        hold on
        ploter(t,e(8,:),plot_type,'g')
        hold on
        ploter(t,e(12,:),plot_type,'b')    
        hold on
        ploter(t,e(15,:),plot_type,'c') 
        hold on
        ploter(t,e(16,:),plot_type,'m')        
        legend("E-STD","E-HD-STD","E-QD","E-AO-STD","E-HDD","E-GHDD",'FontSize',18)
        xlabel('Time (s)','Interpreter','latex','FontSize',18)
        grid on
        set(gca,'FontSize',18)   
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(A)','FontSize',24)         
        
        subplot(2,2,2)
        ploter(t,e(4,:),plot_type,'k')
        hold on
        ploter(t,e(7,:),plot_type,'r')
        hold on        
        ploter(t,e(9,:),plot_type,'g')
        hold on
        ploter(t,e(13,:),plot_type,'b') 
        hold on        
        ploter(t,e(17,:),plot_type,'c')
        hold on
        ploter(t,e(18,:),plot_type,'m')    
        hold on     
        ploter(t,e(22,:),plot_type,[0.8500, 0.3250, 0.0980])  
        hold on
        ploter(t,e(23,:),plot_type,[0.4940, 0.1840, 0.5560])               
        legend("I-STD","I-HD-STD","I-QD","I-AO-STD","I-HDD","I-GHDD","I-FDFF","I-AO-FDFF",'FontSize',18)
        xlabel('Time (s)','Interpreter','latex','FontSize',18)
        grid on
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(B)','FontSize',24)         
        
        subplot(2,2,3)
        ploter(t,e(1,:),plot_type,'r')
        hold on
        ploter(t,e(2,:),plot_type,'g')
        hold on
        ploter(t,e(11,:),plot_type,'b')
        hold on
        ploter(t,e(19,:),plot_type,'c') 
        hold on
        ploter(t,e(21,:),plot_type,'m')            
        legend("Euler","LF","HD","VGED","E-STDAC",'FontSize',18)
        xlabel('Time (s)','Interpreter','latex','FontSize',18)
        grid on
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(C)','FontSize',24)         
        
        subplot(2,2,4)
        ploter(t,e(10,:),plot_type,'r')
        hold on
        ploter(t,e(5,:),plot_type,'g')
        hold on
        ploter(t,e(14,:),plot_type,'b')  
        hold on
        ploter(t,e(20,:),plot_type,'c')
        hold on
        ploter(t,e(24,:),plot_type,'m')          
        legend("ALIEN","SI-STD","SI-AO-STD","SI-HD-STD","Kalman",'FontSize',18)
        xlabel('Time (s)','Interpreter','latex','FontSize',18)
        grid on
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(D)','FontSize',24)         
        
        figure

        sgtitle("Output",'FontSize',18)      
        subplot(2,2,1)
        ploter(t,u(3,:),plot_type,'r')
        hold on    
        ploter(t,u(6,:),plot_type,'g')
        hold on
        ploter(t,u(8,:),plot_type,'b')
        hold on
        ploter(t,u(12,:),plot_type,'c')
        hold on             
        ploter(t,u(15,:),plot_type,'m')    
        hold on     
        ploter(t,u(16,:),plot_type,[.5 0 .5])    
        hold on       
        ploter_dash(t,diff(diff_order,:),plot_type,'k')        
        legend("E-STD","E-HD-STD","E-QD","E-AO-STD","E-HDD","E-GHDD","Real",'FontSize',18)
        xlabel('Time (s)','Interpreter','latex','FontSize',18)
        grid on
        set(gca,'FontSize',18)      
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(A)','FontSize',24)
        
        subplot(2,2,2)
        ploter(t,u(4,:),plot_type,'r')
        hold on
        ploter(t,u(7,:),plot_type,'g')
        hold on        
        ploter(t,u(9,:),plot_type,'b')
        hold on     
        ploter(t,u(13,:),plot_type,'c')
        hold on 
        ploter(t,u(17,:),plot_type,'m')
        hold on     
        ploter(t,u(18,:),plot_type,[.5 0 .5])
        hold on     
        ploter(t,u(22,:),plot_type,[0.8500, 0.3250, 0.0980])   
        hold on 
        ploter(t,u(23,:),plot_type,[0.4940, 0.1840, 0.5560])   
        hold on           
        ploter_dash(t,diff(diff_order,:),plot_type,'k')        
        legend("I-STD","I-HD-STD","I-QD","I-AO-STD","I-HDD","I-GHDD","I-FDFF","I-AO-FDFF","Real",'FontSize',18)
        xlabel('Time (s)','Interpreter','latex','FontSize',18)
        grid on
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(B)','FontSize',24)        
        
        subplot(2,2,3)
        ploter(t,u(1,:),plot_type,'r')
        hold on
        ploter(t,u(2,:),plot_type,'g')
        hold on
        ploter(t,u(11,:),plot_type,'b')  
        hold on
        ploter(t,u(19,:),plot_type,'c')          
        hold on 
        ploter(t,u(21,:),plot_type,'m')          
        hold on            
        ploter_dash(t,diff(diff_order,:),plot_type,'k')        
        legend("Euler","LF","HD","VGED","E-STDAC","Real",'FontSize',18)
        xlabel('Time (s)','Interpreter','latex','FontSize',18)
        grid on
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(C)','FontSize',24)        
        
        subplot(2,2,4)
        ploter(t,u(10,:),plot_type,'r')
        hold on     
        ploter(t,u(5,:),plot_type,'g')
        hold on          
        ploter(t,u(14,:),plot_type,'b')
        hold on
        ploter(t,u(20,:),plot_type,'c')
        hold on  
        ploter(t,u(24,:),plot_type,'m')
        hold on           
        ploter_dash(t,diff(diff_order,:),plot_type,'k')        
        legend("ALIEN","SI-STD","SI-AO-STD","SI-HD-STD","Kalman","Real",'FontSize',18)
        xlabel('Time (s)','Interpreter','latex','FontSize',18)
        grid on
        set(gca,'FontSize',18)
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(end),ylim(end),'(D)','FontSize',24)        
    end

end

