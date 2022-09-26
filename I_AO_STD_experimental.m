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


function [u,z00] = I_AO_STD_experimental(f_k,h,L,reset_flag,output_order,differentiator_order)
% Arbitrary differentiator order  
    persistent z0 z1 z2 z3 z4 z5 z6 z7

    lambda0=12;lambda1=84.14;lambda2=281.37;lambda3=455.40;lambda4=295.74;lambda5=88.78;lambda6=14.13;lambda7=1.1;

    if reset_flag==1
        z0=0;z1=1;z2=0;z3=-1;z4=0;z5=-1;z6=0;z7=1;
    end
    tic
    
    if differentiator_order==7
        lambda0=12;lambda1=84.14;lambda2=281.37;lambda3=455.40;lambda4=295.74;lambda5=88.78;lambda6=14.13;lambda7=1.1;
        a1=h*lambda0*L^(1/8);
        a2=h^2*lambda1*L^(2/8);
        a3=h^3*lambda2*L^(3/8);
        a4=h^4*lambda3*L^(4/8);
        a5=h^5*lambda4*L^(5/8);
        a6=h^6*lambda5*L^(6/8);
        a7=h^7*lambda6*L^(7/8);
        a8=h^8*lambda7*L;
    
        bk=-z0-h*z1-h^2*z2-h^3*z3-h^4*z4-h^5*z5-h^6*z6-h^7*z7+f_k;
    
        if bk<-a8
            sol=roots([1  a1  a2   a3  a4 a5 a6  a7    a8+bk]);
            for i=1:length(sol)
                if isreal(sol(i))
                    XX=sol(i);
                end
            end
            X=abs(XX);
        
        elseif bk>a8
            sol=roots([-1  -a1  -a2   -a3  -a4 -a5 -a6 -a7  bk-a8]);
            for i=1:length(sol)
                if isreal(sol(i))
                    XX=sol(i);
                end
            end
            X=-abs(XX);
        
        else
            X=0;
        end
    
        s=X;
    
        if bk<-a8
        
            z7=h*(-lambda7*L)+z7;
            z6=h*(-lambda6*L^(7/8)*abs(s)^1 +z7  )+z6;
            z5=h*(-lambda5*L^(6/8)*abs(s)^2 +z6  )+z5;
            z4=h*(-lambda4*L^(5/8)*abs(s)^3 +z5  )+z4;
            z3=h*(-lambda3*L^(4/8)*abs(s)^4 +z4  )+z3;
            z2=h*(-lambda2*L^(3/8)*abs(s)^5 +z3  )+z2;
            z1=h*(-lambda1*L^(2/8)*abs(s)^6 +z2  )+z1;
            z0=h*(-lambda0*L^(1/8)*abs(s)^7 +z1   )+z0;
        
        elseif bk>a8
           
            z7=h*(lambda7*L)+z7;
            z6=h*(lambda6*L^(7/8)*abs(s)^1 +z7  )+z6;  
            z5=h*(lambda5*L^(6/8)*abs(s)^2 +z6  )+z5;   
            z4=h*(lambda4*L^(5/8)*abs(s)^3 +z5  )+z4;
            z3=h*(lambda3*L^(4/8)*abs(s)^4 +z4  )+z3;
            z2=h*(lambda2*L^(3/8)*abs(s)^5 +z3  )+z2;
            z1=h*(lambda1*L^(2/8)*abs(s)^6 +z2  )+z1;
            z0=h*(lambda0*L^(1/8)*abs(s)^7 +z1  )+z0;
        
        else
          
            z7=z7+h*(bk);
            z6=z6+h*z7;  
            z5=z5+h*z6; 
            z4=z4+h*z5;
            z3=z3+h*z4;
            z2=z2+h*z3;
            z1=z1+h*z2;
            z0=z0+h*z1;
        
        end
    
    elseif differentiator_order==6
        lambda0=10;lambda1=47.69;lambda2=110.08;lambda3=101.96;lambda4=43.65;lambda5=9.91;lambda6=1.1;
        a1=h*lambda0*L^(1/7);
        a2=h^2*lambda1*L^(2/7);
        a3=h^3*lambda2*L^(3/7);
        a4=h^4*lambda3*L^(4/7);
        a5=h^5*lambda4*L^(5/7);
        a6=h^6*lambda5*L^(6/7);
        a7=h^7*lambda6*L;
    
        bk=-z0-h*z1-h^2*z2-h^3*z3-h^4*z4-h^5*z5-h^6*z6+f_k;
    
        if bk<-a7
            sol=roots([1  a1  a2   a3  a4 a5 a6   a7+bk]);
            for i=1:length(sol)
                if isreal(sol(i))
                    XX=sol(i);
                end
            end
            X=abs(XX);
        
        elseif bk>a7
        sol=roots([-1  -a1  -a2   -a3  -a4 -a5 -a6  bk-a7]);
        for i=1:length(sol)
            if isreal(sol(i))
                XX=sol(i);
            end
        end
        X=-abs(XX);
        
        else
        X=0;
        end
    
        s=X;
    
        if bk<-a7
            z6=h*(-lambda6*L)+z6;
            z5=h*(-lambda5*L^(6/7)*abs(s)^1 +z6  )+z5;
            z4=h*(-lambda4*L^(5/7)*abs(s)^2 +z5  )+z4;
            z3=h*(-lambda3*L^(4/7)*abs(s)^3 +z4  )+z3;
            z2=h*(-lambda2*L^(3/7)*abs(s)^4 +z3  )+z2;
            z1=h*(-lambda1*L^(2/7)*abs(s)^5 +z2  )+z1;
            z0=h*(-lambda0*L^(1/7)*abs(s)^6 +z1   )+z0;    
        elseif bk>a7
            z6=h*(lambda6*L)+z6;
            z5=h*(lambda5*L^(6/7)*abs(s)^1 +z6  )+z5;  
            z4=h*(lambda4*L^(5/7)*abs(s)^2 +z5  )+z4;
            z3=h*(lambda3*L^(4/7)*abs(s)^3 +z4  )+z3;
            z2=h*(lambda2*L^(3/7)*abs(s)^4 +z3  )+z2;
            z1=h*(lambda1*L^(2/7)*abs(s)^5 +z2  )+z1;
            z0=h*(lambda0*L^(1/7)*abs(s)^6 +z1  )+z0;
        else
            z6=z6+h*(bk);
            z5=z5+h*z6;    
            z4=z4+h*z5;
            z3=z3+h*z4;
            z2=z2+h*z3;
            z1=z1+h*z2;
            z0=z0+h*z1;
        end
        
    elseif differentiator_order==5
        lambda0=7;lambda1=23.72;lambda2=32.24;lambda3=20.26;lambda4=6.75;lambda5=1.1;        
        a1=h*lambda0*L^(1/6);
        a2=h^2*lambda1*L^(2/6);
        a3=h^3*lambda2*L^(3/6);
        a4=h^4*lambda3*L^(4/6);
        a5=h^5*lambda4*L^(5/6);
        a6=h^6*lambda5*L;

        bk=-z0-h*z1-h^2*z2-h^3*z3-h^4*z4-h^5*z5+f_k;
    
        if bk<-a6
            sol=roots([1  a1  a2   a3  a4 a5   a6+bk]);
            for i=1:length(sol)
                if isreal(sol(i))
                    XX=sol(i);
                end
            end
            X=abs(XX);
        
        elseif bk>a6
            sol=roots([-1  -a1  -a2   -a3  -a4 -a5  bk-a6]);
            for i=1:length(sol)
                if isreal(sol(i))
                    XX=sol(i);
                end
            end
            X=-abs(XX);
        else
            X=0;
        end
    
        s=X;
    
        if bk<-a6
            z5=h*(-lambda5*L)+z5;
            z4=h*(-lambda4*L^(5/6)*abs(s)^1 +z5  )+z4;
            z3=h*(-lambda3*L^(4/6)*abs(s)^2 +z4  )+z3;
            z2=h*(-lambda2*L^(3/6)*abs(s)^3 +z3  )+z2;
            z1=h*(-lambda1*L^(2/6)*abs(s)^4 +z2  )+z1;
            z0=h*(-lambda0*L^(1/6)*abs(s)^5 +z1   )+z0;
        elseif bk>a6
            z5=h*(lambda5*L)+z5;
            z4=h*(lambda4*L^(5/6)*abs(s)^1 +z5  )+z4;  
            z3=h*(lambda3*L^(4/6)*abs(s)^2 +z4  )+z3;
            z2=h*(lambda2*L^(3/6)*abs(s)^3 +z3  )+z2;
            z1=h*(lambda1*L^(2/6)*abs(s)^4 +z2  )+z1;
            z0=h*(lambda0*L^(1/6)*abs(s)^5 +z1  )+z0;
        else
            z5=z5+h*(bk);
            z4=z4+h*z5;    
            z3=z3+h*z4;
            z2=z2+h*z3;
            z1=z1+h*z2;
            z0=z0+h*z1;
        end
    
	elseif differentiator_order==4
        lambda0=5;lambda1=10.03;lambda2=9.3;lambda3=4.57;lambda4=1.1;
        a1=h*lambda0*L^(1/5);
        a2=h^2*lambda1*L^(2/5);
        a3=h^3*lambda2*L^(3/5);
        a4=h^4*lambda3*L^(4/5);
        a5=h^5*lambda4*L;
    
        bk=-z0-h*z1-h^2*z2-h^3*z3-h^4*z4+f_k;
    
        if bk<-a5
            sol=roots([1  a1  a2   a3  a4   a5+bk]);
            for i=1:length(sol)
                if isreal(sol(i)) && sol(i)>0
                    XX=sol(i);
                end
            end
                X=XX;

        elseif bk>a5
        sol=roots([-1  -a1  -a2   -a3  -a4    bk-a5]);
            for i=1:length(sol)
                if isreal(sol(i)) && sol(i)>0
                    XX=sol(i);
                end
            end
            X=-XX;
        
        else
            X=0;
        end
    
        s=X;
    
        if bk<-a5

            z4=h*(-lambda4*L)+z4;
            z3=h*(-lambda3*L^(4/5)*abs(s) +z4  )+z3;        
            z2=h*(-lambda2*L^(3/5)*abs(s)^2 +z3  )+z2;      
            z1=h*(-lambda1*L^(2/5)*abs(s)^3 +z2  )+z1;         
            z0=h*(-lambda0*L^(1/5)*abs(s)^4 +z1   )+z0;         
         
        elseif bk>a5
        
            z4=h*(lambda4*L)+z4;
            z3=h*(lambda3*L^(4/5)*abs(s) +z4  )+z3;        
            z2=h*(lambda2*L^(3/5)*abs(s)^2 +z3  )+z2;        
            z1=h*(lambda1*L^(2/5)*abs(s)^3 +z2  )+z1;       
            z0=h*(lambda0*L^(1/5)*abs(s)^4 +z1  )+z0;        
        
        else
        
            z4=z4+h*(bk);
            z3=z3+h*z4;        
            z2=z2+h*z3;        
            z1=z1+h*z2;      
            z0=z0+h*z1;
        end
    
	elseif differentiator_order==3
        lambda0=3;lambda1=4.16;lambda2=3.06;lambda3=1.1;        
        a1=h*lambda0*L^(1/4);
        a2=h^2*lambda1*L^(2/4);
        a3=h^3*lambda2*L^(3/4);
        a4=h^4*lambda3*L;
    
        bk=-z0-h*z1-h^2*z2-h^3*z3+f_k;
    
        if bk<-a4
           sol=roots([1  a1  a2   a3    a4+bk]);
            for i=1:length(sol)
                if isreal(sol(i))
                    XX=sol(i);
                end
            end
            X=abs(XX);

        elseif bk>a4
           sol=roots([-1  -a1  -a2   -a3    bk-a4]);
            for i=1:length(sol)
                if isreal(sol(i))
                    XX=sol(i);
                end
            end
            X=-abs(XX);

        else
            X=0;
        end

        s=X;

        if bk<-a4

            z3=h*(-lambda3*L)+z3;
            z2=h*(-lambda2*L^(3/4)*abs(s)^1 +z3  )+z2;        
            z1=h*(-lambda1*L^(2/4)*abs(s)^2 +z2  )+z1;       
            z0=h*(-lambda0*L^(1/4)*abs(s)^3 +z1   )+z0;
        elseif bk>a4
            z3=h*(lambda3*L)+z3;
            z2=h*(lambda2*L^(3/4)*abs(s)^1 +z3  )+z2;        
            z1=h*(lambda1*L^(2/4)*abs(s)^2 +z2  )+z1;        
            z0=h*(lambda0*L^(1/4)*abs(s)^3 +z1  )+z0;        
        else
            z3=z3+h*(bk);
            z2=z2+h*z3;
            z1=z1+h*z2;
            z0=z0+h*z1;
        end
        
    elseif differentiator_order==2
        lambda0=2;lambda1=2.12;lambda2=1.1;        
        a1=h*lambda0*L^(1/3);
        a2=h^2*lambda1*L^(2/3);
        a3=h^3*lambda2*L;

        bk=-z0-h*z1-h^2*z2+f_k;

        if bk<-a3
           sol=roots([1  a1  a2    a3+bk]);
            for i=1:length(sol)
                if isreal(sol(i))
                    XX=sol(i);
                end
            end
            X=abs(XX);

        elseif bk>a3
           sol=roots([-1  -a1  -a2     bk-a3]);
            for i=1:length(sol)
                if isreal(sol(i))
                    XX=sol(i);
                end
            end
            X=-abs(XX);

        else
            X=0;
        end

        s=X;

        if bk<-a3
            z2=h*(-lambda2*L)+z2;
            z1=h*(-lambda1*L^(2/3)*abs(s)^1 +z2  )+z1;        
            z0=h*(-lambda0*L^(1/3)*abs(s)^2 +z1   )+z0;       
        elseif bk>a3
            z2=h*(lambda2*L)+z2;
            z1=h*(lambda1*L^(2/3)*abs(s)^1 +z2  )+z1;        
            z0=h*(lambda0*L^(1/3)*abs(s)^2 +z1  )+z0;        

        else
            z2=z2+h*(bk);
            z1=z1+h*z2;
            z0=z0+h*z1;
	end
    
    else
        disp("Error in E-AO-STD")
    end

    if output_order==0
        u=z0;
    elseif output_order==1
        u=z1;
    elseif output_order==2
        u=z2;
    elseif output_order==3
        u=z3;
    elseif output_order==4
        u=z4;     
    elseif output_order==5
        u=z5;     
    elseif output_order==6
        u=z6;     
    elseif output_order==7
        u=z7;          
    else
        disp("Error in E-AO-STD")
    end

    z00=z0; 
end

