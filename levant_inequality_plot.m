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
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%


clc;clear;close all

tf=10;
h_vector=[1e-4 1e-3 1e-2 2e-2 3e-2 4e-2 5e-2 6e-2 7e-2 8e-2 9e-2 1e-1];
x0=[0 1 0 -1 0 1 0 -1];

L=20;
L_E_STD=L;
L_I_STD=L;
lambdan=1.1;
diff_order=1;

AO_order=3;

for j=1:length(h_vector)
    h=h_vector(j);
    len=floor(tf/h);
    
    reset_flag=1;
    for k=1:len
        
        f(k)=sin(h*k);
        df(k)=cos(h*k);
        
        u_E(k) = E_STD(L_E_STD,h,f(k),reset_flag,diff_order,x0); 
        
        u_I(k) = I_STD(L_I_STD,h,f(k),reset_flag,diff_order,x0);
        reset_flag=0;
        
    end
    
    Linf_E(j)=norm(u_E-df,inf);
    Linf_I(j)=norm(u_I-df,inf);
    clear u_E u_I df
    b(j)=L*h*lambdan;    
end



subplot(2,2,1)
loglog(h_vector,Linf_E,'color','r','LineWidth',2)
hold on
loglog(h_vector,Linf_I,'color','b','LineWidth',2)
hold on
loglog(h_vector,b,'color','k','LineWidth',2,'LineStyle','--')
grid on
set(gca,'FontSize',24)
xlabel('$h$','Interpreter','latex','FontSize',28)
ylabel('$L_\infty(z_{1,k}-f_k^{(1)})$','Interpreter','latex','FontSize',28)

legend("E-STD","I-STD","error band",'FontSize',18)

ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'$f(t)=$sin($t$)','FontSize',24,'Interpreter','latex') 
text(xlim(end),ylim(end),'(A)','Interpreter','latex','FontSize',24)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0=[3 0 0 0];

for j=1:length(h_vector)
    h=h_vector(j);
    len=floor(tf/h);
    
    reset_flag=1;
    for k=1:len
        
        f(k)=5*(h*k)^2+3;
        df(k)=10*(h*k);
        
        u_E(k) = E_STD(L_E_STD,h,f(k),reset_flag,diff_order,x0); 
        
        u_I(k) = I_STD(L_I_STD,h,f(k),reset_flag,diff_order,x0);
        reset_flag=0;
        
    end
    
    Linf_E(j)=norm(u_E-df,inf);
    Linf_I(j)=norm(u_I-df,inf);
    clear u_E u_I df
    b(j)=L*h*lambdan;    
end



subplot(2,2,2)
loglog(h_vector,Linf_E,'color','r','LineWidth',2)
hold on
loglog(h_vector,Linf_I,'color','b','LineWidth',2)
hold on
loglog(h_vector,b,'color','k','LineWidth',2,'LineStyle','--')
grid on
set(gca,'FontSize',24)
xlabel('$h$','Interpreter','latex','FontSize',28)
ylabel('$L_\infty(z_{1,k}-f_k^{(1)})$','Interpreter','latex','FontSize',28)

legend("E-STD","I-STD","error band",'FontSize',18)

ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'$f(t)=5t^2+3$','FontSize',24,'Interpreter','latex')
text(xlim(end),ylim(end),'(B)','Interpreter','latex','FontSize',24)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x0=[3 6 0 0];

for j=1:length(h_vector)
    h=h_vector(j);
    len=floor(tf/h);
    
    reset_flag=1;
    for k=1:len
        
        f(k)=2*sin(k*h)+3*cos(2*k*h)+4*(k*h);
        df(k)=2*cos(k*h)-6*sin(2*k*h)+4;
        
        u_E(k) = E_STD(L_E_STD,h,f(k),reset_flag,diff_order,x0); 
        
        u_I(k) = I_STD(L_I_STD,h,f(k),reset_flag,diff_order,x0);
        reset_flag=0;
        
    end
    
    Linf_E(j)=norm(u_E-df,inf);
    Linf_I(j)=norm(u_I-df,inf);
    clear u_E u_I df
    b(j)=L*h*lambdan;    
end


subplot(2,2,3)
loglog(h_vector,Linf_E,'color','r','LineWidth',2)
hold on
loglog(h_vector,Linf_I,'color','b','LineWidth',2)
hold on
loglog(h_vector,b,'color','k','LineWidth',2,'LineStyle','--')
grid on
set(gca,'FontSize',24)
xlabel('$h$','Interpreter','latex','FontSize',28)
ylabel('$L_\infty(z_{1,k}-f_k^{(1)})$','Interpreter','latex','FontSize',28)

legend("E-STD","I-STD","error band",'FontSize',18)

ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'$f(t)=$2sin($t$)+3cos($2t$)+$4t$','FontSize',24,'Interpreter','latex') 

text(xlim(end),ylim(end),'(C)','Interpreter','latex','FontSize',24)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x0=[0 0 0 0];

for j=1:length(h_vector)
    h=h_vector(j);
    len=floor(tf/h);
    
    reset_flag=1;
    for k=1:len
        
        f(k)=(k*h)^3;
        df(k)=3*(k*h)^2;
        
        u_E(k) = E_STD(L_E_STD,h,f(k),reset_flag,diff_order,x0); 
        
        u_I(k) = I_STD(L_I_STD,h,f(k),reset_flag,diff_order,x0);
        reset_flag=0;
        
    end
    
    Linf_E(j)=norm(u_E-df,inf);
    Linf_I(j)=norm(u_I-df,inf);
    clear u_E u_I df
    b(j)=L*h*lambdan;    
end


subplot(2,2,4)
loglog(h_vector,Linf_E,'color','r','LineWidth',2)
hold on
loglog(h_vector,Linf_I,'color','b','LineWidth',2)
hold on
loglog(h_vector,b,'color','k','LineWidth',2,'LineStyle','--')
grid on
set(gca,'FontSize',24)
xlabel('$h$','Interpreter','latex','FontSize',28)
ylabel('$L_\infty(z_{1,k}-f_k^{(1)})$','Interpreter','latex','FontSize',28)

legend("E-STD","I-STD","error band",'FontSize',18)

ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'$f(t)=t^3$','FontSize',24,'Interpreter','latex') 

text(xlim(end),ylim(end),'(D)','Interpreter','latex','FontSize',24)