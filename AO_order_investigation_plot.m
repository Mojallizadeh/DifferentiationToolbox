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


clear;close all;clc
% SNR=30;
h=[1 5 10 20 30 40 50 60 70 80];
%%%%%%%%%%%%%%%%%%%%%%  1st
L2b_1=[];
L2t_1=[];
L_inf_1=[];
var_1=[];
thd_1=[];
j_1=[];
param_L_1=[];
time_1=[];
%------------------------- 2nd
L2b_2=[0.0004 0.0012 0.0021 0.0029 0.0043 0.0057 0.0068 0.0068 0.0074 0.0112];
L2t_2=[12.7054 16.6977 21.0628 20.5783 24.2004 28.3282 29.6993 27.1834 27.0031 38.8531];
L_inf_2=[0.0963 0.1344 0.1648 0.1839 0.1934 0.2125 0.2767 0.1928 0.2167 0.3435];
var_2=[14.5381 12.0183 11.2501 9.7918 9.8873 9.0275 9.0045 8.0849 8.1296 8.6361];
thd_2=[0.1441 0.6175 1.1126 2.0103 5.7753 3.6432 4.3687 5.0020 6.1492 12.4149];
j_2=[27.5242 29.5863 33.8015 32.8579 40.4841 41.7856 44.0245 41.1447 42.2400 61.3722];
param_L_2=[4.7469 3.5960 3.2988 2.7140 2.8757 2.4423 2.6199 2.2906 2.2906 2.4653];
time_2=[145.6470 30.1650 15.9660 9.7370 6.6420 5.1560 3.3270 3.2440 2.9690 3.3560];
%------------------------- 3rd
L2b_3=[0.0004 0.0010 0.0017 0.0028 0.0045 0.0062 0.0079 0.0075 0.0092 0.0101];
L2t_3=[11.3210 13.5549 16.4762 19.4245 25.44103 32 34.7026 29.8810 33.8137 34.7423]; %--------------
L_inf_3=[0.0785 0.1249 0.1560 0.1671 0.2088 0.2470 0.3387 0.2528 0.2772 0.2661];
var_3=[12.2185 12.0921 11.0959 10.5177 9.2140 8.2807 8.8245 8.3368 8.1188 7.7077];
thd_3=[0.1457 0.6120 1.1136 2.0934 5.6198 3.7739 4.5780 5.2075 5.9702 12.0586];
j_3=[23.7995 26.4799 29.0073 32.4803 40.9292 43.8193 49.2291 44.4283 49.0952 55.7858];
param_L_3=[3.9213 5.8064 4.9887 4.7430 3.8781 2.5504 3.8921 3.6514 2.9112 2.1139];
time_3=[177.0440 35.9640 17.6860 9.0690 6.3230 4.9220 4.1750 3.9510 2.8630 0.6640];
%------------------------- 4th
L2b_4=[0.0004 0.0010 0.0019 0.0032 0.0047 0.0053 0.0069 0.0085 0.0075 0.0103];      
L2t_4=[12.5415 14.4935 18.7757 21.9235 26.9743 25.4737 29.8450 33.7360 27.0916 34.9496];
L_inf_4=[0.0966 0.1358 0.1909 0.2111 0.2415 0.2372 0.3046 0.2638 0.2828 0.2910];
var_4=[16.7846 16.5693 14.0845 13.5114 10.6126 11.1516 10.6947 8.8704 8.5149 9.5643]; 
thd_4=[0.1459 0.6207 1.1239 2.0851 5.5855 3.8330 4.8532 5.1531 6.0554 13.4193];
j_4=[29.6083 31.9222 34.3640 38.0465 43.8877 41.2218 46.3866 48.8721 42.6914 59.2541];
param_L_4=[6.0451 6.4572 4.7559 4.6037 3.0266 4.6028 5.0567 2.2616 3.2091 5.0240];
time_4=[220.1070 42.3380 20.6010 10.0730 9.6300 5.1190 4.7030 3.5860 3.5280 3.2540] ;
%------------------------- 5th
L2b_5=[0.0005 0.0013 0.0025 0.0033 0.0050 0.0066 0.0073 0.0097 0.0111 0.0120]; 
L2t_5=[15.8936 18.2833 24.7670  22.8021 28.0114 31.8638 31.5172 38.1450 39.7424 38.6349];
L_inf_5=[0.1115 0.1621 0.2119  0.2386 0.2782 0.2789 0.2922 0.3066 0.3252 0.3676];  
var_5=[18.2313 18.6170 14.9753  14.2347 14.3154 13.2046 11.2847 10.7768 10.9103 11.4695] ;    
thd_5=[0.1436 0.6078 1.1398  2.0999 6.0028 3.9258 4.5792 5.3617 6.7363 12.5530];               
j_5=[34.4302 37.7999 41.3426  39.7035 49.1057 49.9287 48.3995 55.5561 58.8290 64.2232];
param_L_5=[1.8839 2.0304 1.2418  1.5103 1.4747 2.0320 1.3574 1.2010 1.9552 2.8289];    
time_5=[272.3870 57.4860 26.3100 13.5680 8.6410 6.7670  5.3040 5.0880 4.5910 3.0390];
%------------------------- 6th
L2b_6=[0.0004 0.0012 0.0020 0.0036 0.0045 0.0062 0.0064 0.0089 0.0111 0.0113]; 
L2t_6=[12.2283 17.1676 19.8814 25.2579 25.3032 29.6888 27.6356 34.6659 40.1068 38.4951];
L_inf_6=[0.1017 0.1754 0.1844 0.2407 0.2303 0.3161 0.2424 0.2920 0.3158 0.3046]; 
var_6=[17.2913 15.6457 14.7694 13.1930 10.8777 12.9430 10.3336 9.8555 11.0150 9.1765] ;
thd_6=[0.1426 0.6046 1.1076 2.0374 5.3484 3.6891 4.4965 5.2217 6.1238 12.0134];
j_6=[29.8026 33.7150 36.1429 41.0907 42.2054 47.2527  43.3461 50.9226 58.6726 61.1160];
param_L_6=[0.2702 0.1948 0.1948 0.1449 0.1012 0.2029 0.1580 0.1292 0.3425 0.1163]; 
time_6=[300.6560 64.4500 30.9700 13.8360 10.0410 7.3940 5.9030 5.8470 4.5860 4.1800];
%------------------------- 7th
L2b_7=[0.0005 0.0020 0.0038 0.0058 0.0076 0.0094 0.0136 0.0130 0.0176 0.0235];
L2t_7=[14.8194 28.3582 37.5522 40.5068 42.3667 45.2692 58.8494 49.9539 63.9521 81.6813];
L_inf_7=[0.1374 0.2674 0.3484 0.4300 0.3748 0.4107 0.5681 0.4697 0.5302 0.6396];  
var_7=[38.5226 30.5271 24.9442 25.4504 20.7315 18.9824 19.4108 15.9218 14.2684 13.4303];
thd_7=[0.1473 0.6335 1.1347 2.1404 4.8849 3.9334 5.2177 5.6675 6.6496 12.1265];
j_7=[53.6736 59.9875 64.3573 69.1122 69.1181 69.5315 85.4031 73.3135 87.1575 110.2311];
param_L_7=[1.1064 0.5600 0.2764 0.4473 0.2529 0.2529 0.4620 0.3163 0.2491 0.1728];
time_7=[347.5870 66.0220 33.7330 16.7750 10.5600 10.4840 6.8150 6.0840 6.1340 5.3320];











subplot(4,2,1)

% loglog(h,L2b_1,'LineWidth',2)
hold on
loglog(h,L2b_2,'LineWidth',2)
hold on
loglog(h,L2b_3,'LineWidth',2)
hold on
loglog(h,L2b_4,'LineWidth',2)
hold on
loglog(h,L2b_5,'LineWidth',2)
hold on
loglog(h,L2b_6,'LineWidth',2)
hold on
loglog(h,L2b_7,'LineWidth',2)

legend("2nd order","3rd order","4th order","5th order","6th order","7th order",'FontSize',18)
ylabel('$\bar{L}_2$','Interpreter','latex','FontSize',18)
grid on
set(gca,'FontSize',18)

subplot(4,2,2)

% loglog(h,L2t_1,'LineWidth',2)
hold on
loglog(h,L2t_2,'LineWidth',2)
hold on
loglog(h,L2t_3,'LineWidth',2)
hold on
loglog(h,L2t_4,'LineWidth',2)
hold on
loglog(h,L2t_5,'LineWidth',2)
hold on
loglog(h,L2t_6,'LineWidth',2)
hold on
loglog(h,L2t_7,'LineWidth',2)

legend("2nd order","3rd order","4th order","5th order","6th order","7th order",'FontSize',18)
ylabel('$\tilde{L}_2$','Interpreter','latex','FontSize',18)
grid on
set(gca,'FontSize',18)

subplot(4,2,3)

% loglog(h,L_inf_1,'LineWidth',2)
hold on
loglog(h,L_inf_2,'LineWidth',2)
hold on
loglog(h,L_inf_3,'LineWidth',2)
hold on
loglog(h,L_inf_4,'LineWidth',2)
hold on
loglog(h,L_inf_5,'LineWidth',2)
hold on
loglog(h,L_inf_6,'LineWidth',2)
hold on
loglog(h,L_inf_7,'LineWidth',2)

legend("2nd order","3rd order","4th order","5th order","6th order","7th order",'FontSize',18)
ylabel('$L_\infty$','Interpreter','latex','FontSize',18)
grid on
set(gca,'FontSize',18)

subplot(4,2,4)

% loglog(h,var_1,'LineWidth',2)
hold on
loglog(h,var_2,'LineWidth',2)
hold on
loglog(h,var_3,'LineWidth',2)
hold on
loglog(h,var_4,'LineWidth',2)
hold on
loglog(h,var_5,'LineWidth',2)
hold on
loglog(h,var_6,'LineWidth',2)
hold on
loglog(h,var_7,'LineWidth',2)

legend("2nd order","3rd order","4th order","5th order","6th order","7th order",'FontSize',18)
ylabel('VAR','Interpreter','latex','FontSize',18)
grid on
set(gca,'FontSize',18)

subplot(4,2,5)

% loglog(h,thd_1,'LineWidth',2)
hold on
loglog(h,thd_2,'LineWidth',2)
hold on
loglog(h,thd_3,'LineWidth',2)
hold on
loglog(h,thd_4,'LineWidth',2)
hold on
loglog(h,thd_5,'LineWidth',2)
hold on
loglog(h,thd_6,'LineWidth',2)
hold on
loglog(h,thd_7,'LineWidth',2)

legend("2nd order","3rd order","4th order","5th order","6th order","7th order",'FontSize',18)
ylabel('THD','Interpreter','latex','FontSize',18)
grid on
set(gca,'FontSize',18)


subplot(4,2,6)

% loglog(h,j_1,'LineWidth',2)
hold on
loglog(h,j_2,'LineWidth',2)
hold on
loglog(h,j_3,'LineWidth',2)
hold on
loglog(h,j_4,'LineWidth',2)
hold on
loglog(h,j_5,'LineWidth',2)
hold on
loglog(h,j_6,'LineWidth',2)
hold on
loglog(h,j_7,'LineWidth',2)

legend("2nd order","3rd order","4th order","5th order","6th order","7th order",'FontSize',18)
ylabel('$J$','Interpreter','latex','FontSize',18)
grid on
set(gca,'FontSize',18)


subplot(4,2,7)

% loglog(h,time_1,'LineWidth',2)
hold on
loglog(h,time_2,'LineWidth',2)
hold on
loglog(h,time_3,'LineWidth',2)
hold on
loglog(h,time_4,'LineWidth',2)
hold on
loglog(h,time_5,'LineWidth',2)
hold on
loglog(h,time_6,'LineWidth',2)
hold on
loglog(h,time_7,'LineWidth',2)

legend("2nd order","3rd order","4th order","5th order","6th order","7th order",'FontSize',18)
ylabel('Calculation time (s)','Interpreter','latex','FontSize',18)
xlabel('$h$ (ms)','Interpreter','latex','FontSize',18)
grid on
set(gca,'FontSize',18)

subplot(4,2,8)

% loglog(h,param_L_1,'LineWidth',2)
hold on
loglog(h,param_L_2,'LineWidth',2)
hold on
loglog(h,param_L_3,'LineWidth',2)
hold on
loglog(h,param_L_4,'LineWidth',2)
hold on
loglog(h,param_L_5,'LineWidth',2)
hold on
loglog(h,param_L_6,'LineWidth',2)
hold on
loglog(h,param_L_7,'LineWidth',2)

legend("2nd order","3rd order","4th order","5th order","6th order","7th order",'FontSize',18)
ylabel('$L$','Interpreter','latex','FontSize',18)
xlabel('$h$ (ms)','Interpreter','latex','FontSize',18)
grid on
set(gca,'FontSize',18)









