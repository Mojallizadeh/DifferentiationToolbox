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

function [f,df,ddf,dddf] = signal_generator(len,w,amp,phi,SNR,h,type,An,Wn,Quant)
%% Generate the input signal and its first three differentiation
% A sinusoidal signal is considered to investigate the THD effectively
    f_t=zeros(1,len+1);
    f=zeros(1,len+1);
    df=zeros(1,len+1);
    ddf=zeros(1,len+1);
    dddf=zeros(1,len+1);
    
    for k=1:floor(len+1)
        
        f_t(k)=amp*sin(w*h*k+phi); % Signal
        df(k)=amp*w*cos(w*h*k+phi); % Real 1st-order differentiation
        ddf(k)=-amp*w^2*sin(w*h*k+phi); % Real 2nd-order differentiation
        dddf(k)=-amp*w^3*cos(w*h*k+phi); % Real 3rd-order differentiation
    end
    
    %% Add noise
    if type=="white"
        f=awgn(f_t,SNR,'measured',67);
    elseif type=="sin"
        for k=1:floor(len+1)
            f(k)=f_t(k)+An*sin(Wn*h*k);
        end
        
    elseif type=="bell"
       n=awgn(f_t,SNR,'measured',67)-f_t;
        Z1=1/5*bandpass(n,[4 10],1/h);
        Z2=1/5*bandpass(n,[5 9],1/h);
        Z3=1/5*bandpass(n,[6 8],1/h);
        Z4=1/5*bandpass(n,[6.5 7.5],1/h);
        Z5=1/5*bandpass(n,[6.75 7.25],1/h);
        f=Z1+Z2+Z3+Z4+Z5+f_t;

    else
        disp('Undefined noise')
    end
    
    if Quant==1
        f=quant(f,0.1); % For quantization
%         f=quant_custom(f); 
    end
    
end

