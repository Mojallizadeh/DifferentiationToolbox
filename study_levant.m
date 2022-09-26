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

clear;close all;clc

t_f=50;
L=24;
sample_length=10;
diff_order=1;
h_vector=linspace(1e-6,1e-3,sample_length);
% h_vector=1e-5;

for i=1:length(h_vector)
    
    h=h_vector(i);
% h=1e-3;
    reset_flag=1;
    
    f=zeros(1,floor(t_f/h)+1);
    df=zeros(1,floor(t_f/h)+1);
    u=zeros(1,floor(t_f/h)+1);
    
    for k=1:floor(t_f/h+1)

        f(k)=(h*k)^4-5*(h*k)^2+2*(h*k); % Signal
        df(k)=4*(h*k)^3-10*(h*k)+2; % Real 1st-order differentiation
% f(k)=1;
% df(k)=0;
%         ddf(k)=12*(h*k)^2-10; % Real 2nd-order differentiation
%         dddf(k)=24*(h*k); % Real 3rd-order differentiation

        u(k)=SI_AO_STD_3_optimized_levant(f(k),h,L,reset_flag,diff_order);
        reset_flag=0;

    end
    
    Linf(i)=min(  abs(  u(end)-df(end)  )  );
    
%     Linf(i)=max(abs(df(floor(8*end/10):end)-u(floor(8*end/10):end)));

end

loglog(h_vector,Linf)
grid on

figure 
plot(u-df)