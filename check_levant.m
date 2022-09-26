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
L_I_STD=20000000;
hd=[1e-5 1e-4 1e-3 1e-2];
output_order=1;
differentiator_order=6;
for i=1:length(hd)
        reset_flag=1;
        h=hd(i);
t=0:h:tf;
f=ones(1,length(t));
% f=10*t.^2+5.*t+3;
% f=sin(t);
fr=zeros(1,length(t));
% fr=cos(t);
% fr=20.*t+5*ones(1,length(t));




    for j=1:length(f)

    [y(j),x(j)] = I_AO_STD_experimental(f(j),h,L_I_STD,reset_flag,output_order,differentiator_order);
    reset_flag=0;

    end
    



e0=abs(x-f);
e0=max(e0(floor(end/2):end))
e0_log(i)=log10(e0);
e1=abs(y-fr);
e1=max(e1(floor(end/2):end))
e1_log(i)=log10(e1);

clear eo e1 f fr x y
end


subplot(3,2,1)
plot(log10(hd),e0_log,'LineWidth',2)
set(gca,'FontSize',18)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'(A)','FontSize',18) 
text(xlim(end),ylim(end),'$f(t)=1$','Interpreter','latex','FontSize',18)
grid on
% xlabel('log($h$)','Interpreter','latex','FontSize',18)
% ylabel('log($L_\infty(|z_0-f_0|)$)','Interpreter','latex','FontSize',18)
text(xlim(end),ylim(end),'log($L_\infty(|z_0-f_0|)$)','Interpreter','latex','FontSize',24) 
subplot(3,2,2)
plot(log10(hd),e1_log,'LineWidth',2)
set(gca,'FontSize',18)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'(B)','FontSize',18)
text(xlim(end),ylim(end),'log($L_\infty(|z_1-f_1|)$)','Interpreter','latex','FontSize',24) 
% text('log($L_\infty(|z_0-f_0|)$)','Interpreter','latex','FontSize',18)
grid on
% xlabel('log($h$)','Interpreter','latex','FontSize',18)
% ylabel('log($L_\infty(|z_1-f^{(1)}_0|)$)','Interpreter','latex','FontSize',18)


slop1=(e1_log(end)-e1_log(1))/(hd(end)-hd(1))
slop0=(e0_log(end)-e0_log(1))/(hd(end)-hd(1))


for i=1:length(hd)
        reset_flag=1;
        h=hd(i);
t=0:h:tf;
% % f=ones(1,length(t));
f=10*t.^2+5.*t+3;
% f=sin(t);
% fr=zeros(1,length(t));
% fr=cos(t);
fr=20.*t+5*ones(1,length(t));




    for j=1:length(f)

    [y(j),x(j)] = E_STD_experimental(L_I_STD,h,f(j),reset_flag);
    reset_flag=0;

    end
    



e0=abs(x-f);
e0=max(e0(floor(end/2):end))
e0_log(i)=log10(e0);
e1=abs(y-fr);
e1=max(e1(floor(end/2):end))
e1_log(i)=log10(e1);

clear eo e1 f fr x y
end


subplot(3,2,3)
plot(log10(hd),e0_log,'LineWidth',2)
set(gca,'FontSize',18)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'(C)','FontSize',18) 
text(xlim(end),ylim(end),'$f(t)=10*t.^2+5.*t+3$','Interpreter','latex','FontSize',24) 
grid on
% xlabel('log($h$)','Interpreter','latex','FontSize',18)
% ylabel('log($L_\infty(|z_0-f_0|)$)','Interpreter','latex','FontSize',18)
subplot(3,2,4)
plot(log10(hd),e1_log,'LineWidth',2)
set(gca,'FontSize',18)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'(D)','FontSize',18)
grid on
% xlabel('log($h$)','Interpreter','latex','FontSize',18)
% ylabel('log($L_\infty(|z_1-f^{(1)}_0|)$)','Interpreter','latex','FontSize',18)


slop1=(e1_log(end)-e1_log(1))/(hd(end)-hd(1))
slop0=(e0_log(end)-e0_log(1))/(hd(end)-hd(1))


for i=1:length(hd)
        reset_flag=1;
        h=hd(i);
t=0:h:tf;
% % f=ones(1,length(t));
% f=10*t.^2+5.*t+3;
f=sin(t);
% fr=zeros(1,length(t));
fr=cos(t);
% fr=20.*t+5*ones(1,length(t));




    for j=1:length(f)

    [y(j),x(j)] = E_STD_experimental(L_I_STD,h,f(j),reset_flag);
    reset_flag=0;

    end
    



e0=abs(x-f);
e0=max(e0(floor(end/2):end))
e0_log(i)=log10(e0);
e1=abs(y-fr);
e1=max(e1(floor(end/2):end))
e1_log(i)=log10(e1);

clear eo e1 f fr x y
end

% hd=log10(hd);

subplot(3,2,5)
plot(log10(hd),e0_log,'LineWidth',2)
set(gca,'FontSize',18)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'(E)','FontSize',18) 
text(xlim(end),ylim(end),'$f(t)$=sin($t$)','Interpreter','latex','FontSize',24) 
grid on
xlabel('log($h$)','Interpreter','latex','FontSize',18)
% % ylabel('log($L_\infty(|z_0-f_0|)$)','Interpreter','latex','FontSize',18)
subplot(3,2,6)
plot(log10(hd),e1_log,'LineWidth',2)
set(gca,'FontSize',18)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'(F)','FontSize',18) 
% text(xlim(end),ylim(end),'(A)','FontSize',24) 
grid on
xlabel('log($h$)','Interpreter','latex','FontSize',18)
% ylabel('log($L_\infty(|z_1-f^{(1)}_0|)$)','Interpreter','latex','FontSize',18)

slop1=(e1_log(end)-e1_log(1))/(hd(end)-hd(1))
slop0=(e0_log(end)-e0_log(1))/(hd(end)-hd(1))

% subplot(2,1,1)
% plot(x)
% hold on
% plot(f)
% subplot(2,1,2)
% plot(y)
% hold on
% plot(fr)



