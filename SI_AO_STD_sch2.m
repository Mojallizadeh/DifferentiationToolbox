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

function [u,calculation_time] = SI_AO_STD_sch2(f_k,h,L,reset_flag,output_order,differentiator_order)
% Arbitrary differentiator order 3
%   Version 29-11-2019  Ok

persistent z0 z1 z2 z3 z4 z5 z6 z7 CT

    if reset_flag==1
        z0=0;z1=1;z2=0;z3=-1;z4=0;z5=1;z6=0;z7=-1;CT=0;
    end
    tic
    
	if differentiator_order==7
        lambda0=12;lambda1=84.14;lambda2=281.37;lambda3=455.40;lambda4=295.74;lambda5=88.78;lambda6=14.13;lambda7=1.1;
        bk=f_k  -z0-h^1*z1-h^2*z2-h^3*z3-h^4*z4-h^5*z5-h^6*z6-h^7*z7  ...
            +h*lambda0*L^((1)/(8))*sign(z0-f_k)*abs(z0-f_k)^((7)/(8))...
            +h^2*lambda1*L^((2)/(8))*sign(z0-f_k)*abs(z0-f_k)^((6)/(8))...
            +h^3*lambda2*L^((3)/(8))*sign(z0-f_k)*abs(z0-f_k)^((5)/(8))...
            +h^4*lambda3*L^((4)/(8))*sign(z0-f_k)*abs(z0-f_k)^((4)/(8))...
            +h^5*lambda4*L^((5)/(8))*sign(z0-f_k)*abs(z0-f_k)^((3)/(8))...
            +h^6*lambda5*L^((6)/(8))*sign(z0-f_k)*abs(z0-f_k)^((2)/(8))...
            +h^7*lambda6*L^((7)/(8))*sign(z0-f_k)*abs(z0-f_k)^((1)/(8));        
        if bk<-h^(8)*lambda7*L
            
            z7=h*(-lambda7*L)+z7;
            z6=h*(-lambda6*L^(7/8)*abs(z0-f_k)^(1/8)*sign(z0-f_k)+z7)+z6;
            z5=h*(-lambda5*L^(6/8)*abs(z0-f_k)^(2/8)*sign(z0-f_k)+z6)+z5;
            z4=h*(-lambda4*L^(5/8)*abs(z0-f_k)^(3/8)*sign(z0-f_k)+z5)+z4;
            z3=h*(-lambda3*L^(4/8)*abs(z0-f_k)^(4/8)*sign(z0-f_k)+z4)+z3;
            z2=h*(-lambda2*L^(3/8)*abs(z0-f_k)^(5/8)*sign(z0-f_k)+z3)+z2;
            z1=h*(-lambda1*L^(2/8)*abs(z0-f_k)^(6/8)*sign(z0-f_k)+z2)+z1;
            z0=h*(-lambda0*L^(1/8)*abs(z0-f_k)^(7/8)*sign(z0-f_k)+z1)+z0;            
            
        elseif bk>h^(8)*lambda7*L
            z7=h*(lambda7*L)+z7;
            z6=h*(-lambda6*L^(7/8)*abs(z0-f_k)^(1/8)*sign(z0-f_k)+z7)+z6;
            z5=h*(-lambda5*L^(6/8)*abs(z0-f_k)^(2/8)*sign(z0-f_k)+z6)+z5;
            z4=h*(-lambda4*L^(5/8)*abs(z0-f_k)^(3/8)*sign(z0-f_k)+z5)+z4;
            z3=h*(-lambda3*L^(4/8)*abs(z0-f_k)^(4/8)*sign(z0-f_k)+z4)+z3;
            z2=h*(-lambda2*L^(3/8)*abs(z0-f_k)^(5/8)*sign(z0-f_k)+z3)+z2;
            z1=h*(-lambda1*L^(2/8)*abs(z0-f_k)^(6/8)*sign(z0-f_k)+z2)+z1;
            z0=h*(-lambda0*L^(1/8)*abs(z0-f_k)^(7/8)*sign(z0-f_k)+z1)+z0;                
            
        else
            
            z7=bk/h+z7;
            z6=h*(-lambda6*L^(7/8)*abs(z0-f_k)^(1/8)*sign(z0-f_k)+z7)+z6;
            z5=h*(-lambda5*L^(6/8)*abs(z0-f_k)^(2/8)*sign(z0-f_k)+z6)+z5;
            z4=h*(-lambda4*L^(5/8)*abs(z0-f_k)^(3/8)*sign(z0-f_k)+z5)+z4;
            z3=h*(-lambda3*L^(4/8)*abs(z0-f_k)^(4/8)*sign(z0-f_k)+z4)+z3;
            z2=h*(-lambda2*L^(3/8)*abs(z0-f_k)^(5/8)*sign(z0-f_k)+z3)+z2;
            z1=h*(-lambda1*L^(2/8)*abs(z0-f_k)^(6/8)*sign(z0-f_k)+z2)+z1;
            z0=h*(-lambda0*L^(1/8)*abs(z0-f_k)^(7/8)*sign(z0-f_k)+z1)+z0;  
            
        end

	elseif differentiator_order==6
        lambda0=10;lambda1=47.69;lambda2=110.08;lambda3=101.96;lambda4=43.65;lambda5=9.91;lambda6=1.1;
        bk=f_k  -z0-h^1*z1-h^2*z2-h^3*z3-h^4*z4-h^5*z5-h^6*z6 ...
            +h*lambda0*L^((1)/(7))*sign(z0-f_k)*abs(z0-f_k)^((6)/(7))...
            +h^2*lambda1*L^((2)/(7))*sign(z0-f_k)*abs(z0-f_k)^((5)/(7))...
            +h^3*lambda2*L^((3)/(7))*sign(z0-f_k)*abs(z0-f_k)^((4)/(7))...
            +h^4*lambda3*L^((4)/(7))*sign(z0-f_k)*abs(z0-f_k)^((3)/(7))...
            +h^5*lambda4*L^((5)/(7))*sign(z0-f_k)*abs(z0-f_k)^((2)/(7))...
            +h^6*lambda5*L^((6)/(7))*sign(z0-f_k)*abs(z0-f_k)^((1)/(7));  
        
         if bk<-h^(7)*lambda6*L
        
            z6=h*(-lambda6*L)+z6;
            z5=h*(-lambda5*L^(6/7)*abs(z0-f_k)^(1/7)*sign(z0-f_k)+z6)+z5;
            z4=h*(-lambda4*L^(5/7)*abs(z0-f_k)^(2/7)*sign(z0-f_k)+z5)+z4; 
            z3=h*(-lambda3*L^(4/7)*abs(z0-f_k)^(3/7)*sign(z0-f_k)+z4)+z3; 
            z2=h*(-lambda2*L^(3/7)*abs(z0-f_k)^(4/7)*sign(z0-f_k)+z3)+z2; 
            z1=h*(-lambda1*L^(2/7)*abs(z0-f_k)^(5/7)*sign(z0-f_k)+z2)+z1; 
            z0=h*(-lambda0*L^(1/7)*abs(z0-f_k)^(6/7)*sign(z0-f_k)+z1)+z0;  
        
         elseif bk>h^(7)*lambda6*L
            z6=h*(lambda6*L)+z6;
            z5=h*(-lambda5*L^(6/7)*abs(z0-f_k)^(1/7)*sign(z0-f_k)+z6)+z5;
            z4=h*(-lambda4*L^(5/7)*abs(z0-f_k)^(2/7)*sign(z0-f_k)+z5)+z4; 
            z3=h*(-lambda3*L^(4/7)*abs(z0-f_k)^(3/7)*sign(z0-f_k)+z4)+z3; 
            z2=h*(-lambda2*L^(3/7)*abs(z0-f_k)^(4/7)*sign(z0-f_k)+z3)+z2; 
            z1=h*(-lambda1*L^(2/7)*abs(z0-f_k)^(5/7)*sign(z0-f_k)+z2)+z1; 
            z0=h*(-lambda0*L^(1/7)*abs(z0-f_k)^(6/7)*sign(z0-f_k)+z1)+z0;               
             
         else
             
            z6=bk/h+z6;
            z5=h*(-lambda5*L^(6/7)*abs(z0-f_k)^(1/7)*sign(z0-f_k)+z6)+z5;
            z4=h*(-lambda4*L^(5/7)*abs(z0-f_k)^(2/7)*sign(z0-f_k)+z5)+z4; 
            z3=h*(-lambda3*L^(4/7)*abs(z0-f_k)^(3/7)*sign(z0-f_k)+z4)+z3; 
            z2=h*(-lambda2*L^(3/7)*abs(z0-f_k)^(4/7)*sign(z0-f_k)+z3)+z2; 
            z1=h*(-lambda1*L^(2/7)*abs(z0-f_k)^(5/7)*sign(z0-f_k)+z2)+z1; 
            z0=h*(-lambda0*L^(1/7)*abs(z0-f_k)^(6/7)*sign(z0-f_k)+z1)+z0;              
         end
         
	elseif differentiator_order==5
        lambda0=7;lambda1=23.72;lambda2=32.24;lambda3=20.26;lambda4=6.75;lambda5=1.1;   
         bk=f_k  -z0-h^1*z1-h^2*z2-h^3*z3-h^4*z4-h^5*z5...
            +h*lambda0*L^((1)/(6))*sign(z0-f_k)*abs(z0-f_k)^((5)/(6))...
            +h^2*lambda1*L^((2)/(6))*sign(z0-f_k)*abs(z0-f_k)^((4)/(6))...
            +h^3*lambda2*L^((3)/(6))*sign(z0-f_k)*abs(z0-f_k)^((3)/(6))...
            +h^4*lambda3*L^((4)/(6))*sign(z0-f_k)*abs(z0-f_k)^((2)/(6))...
            +h^5*lambda4*L^((5)/(6))*sign(z0-f_k)*abs(z0-f_k)^((1)/(6));
        
         if bk<-h^(6)*lambda5*L
        
            z5=h*(-lambda5*L)+z5;
            z4=h*(-lambda4*L^(5/6)*abs(z0-f_k)^(1/6)*sign(z0-f_k)+z5)+z4; 
            z3=h*(-lambda3*L^(4/6)*abs(z0-f_k)^(2/6)*sign(z0-f_k)+z4)+z3; 
            z2=h*(-lambda2*L^(3/6)*abs(z0-f_k)^(3/6)*sign(z0-f_k)+z3)+z2; 
            z1=h*(-lambda1*L^(2/6)*abs(z0-f_k)^(4/6)*sign(z0-f_k)+z2)+z1; 
            z0=h*(-lambda0*L^(1/6)*abs(z0-f_k)^(5/6)*sign(z0-f_k)+z1)+z0;  
        
         elseif bk>h^(6)*lambda5*L
            z5=h*(lambda5*L)+z5;
            z4=h*(-lambda4*L^(5/6)*abs(z0-f_k)^(1/6)*sign(z0-f_k)+z5)+z4; 
            z3=h*(-lambda3*L^(4/6)*abs(z0-f_k)^(2/6)*sign(z0-f_k)+z4)+z3; 
            z2=h*(-lambda2*L^(3/6)*abs(z0-f_k)^(3/6)*sign(z0-f_k)+z3)+z2; 
            z1=h*(-lambda1*L^(2/6)*abs(z0-f_k)^(4/6)*sign(z0-f_k)+z2)+z1; 
            z0=h*(-lambda0*L^(1/6)*abs(z0-f_k)^(5/6)*sign(z0-f_k)+z1)+z0;               
             
         else
             
            z5=bk/h+z5;
            z4=h*(-lambda4*L^(5/6)*abs(z0-f_k)^(1/6)*sign(z0-f_k)+z5)+z4; 
            z3=h*(-lambda3*L^(4/6)*abs(z0-f_k)^(2/6)*sign(z0-f_k)+z4)+z3; 
            z2=h*(-lambda2*L^(3/6)*abs(z0-f_k)^(3/6)*sign(z0-f_k)+z3)+z2; 
            z1=h*(-lambda1*L^(2/6)*abs(z0-f_k)^(4/6)*sign(z0-f_k)+z2)+z1; 
            z0=h*(-lambda0*L^(1/6)*abs(z0-f_k)^(5/6)*sign(z0-f_k)+z1)+z0;              
         end

	elseif differentiator_order==4
    lambda0=5;lambda1=10.03;lambda2=9.3;lambda3=4.57;lambda4=1.1;
        bk=f_k  -z0-h^1*z1-h^2*z2-h^3*z3-h^4*z4 ...
            +h*lambda0*L^((1)/(5))*sign(z0-f_k)*abs(z0-f_k)^((4)/(5))...
            +h^2*lambda1*L^((2)/(5))*sign(z0-f_k)*abs(z0-f_k)^((3)/(5))...
            +h^3*lambda2*L^((3)/(5))*sign(z0-f_k)*abs(z0-f_k)^((2)/(5))...
            +h^4*lambda3*L^((4)/(5))*sign(z0-f_k)*abs(z0-f_k)^((1)/(5));
        
        if bk<-h^(5)*lambda4*L
            z4=h*(-lambda4*L)+z4;
            z3=h*(-lambda3*L^(4/5)*abs(z0-f_k)^(1/5)*sign(z0-f_k)+z4)+z3;
            z2=h*(-lambda2*L^(3/5)*abs(z0-f_k)^(2/5)*sign(z0-f_k)+z3)+z2;
            z1=h*(-lambda1*L^(2/5)*abs(z0-f_k)^(3/5)*sign(z0-f_k)+z2)+z1;
            z0=h*(-lambda0*L^(1/5)*abs(z0-f_k)^(4/5)*sign(z0-f_k)+z1)+z0;            
             
        elseif bk>h^(5)*lambda4*L
            z4=h*(lambda4*L)+z4;
            z3=h*(-lambda3*L^(4/5)*abs(z0-f_k)^(1/5)*sign(z0-f_k)+z4)+z3;
            z2=h*(-lambda2*L^(3/5)*abs(z0-f_k)^(2/5)*sign(z0-f_k)+z3)+z2;
            z1=h*(-lambda1*L^(2/5)*abs(z0-f_k)^(3/5)*sign(z0-f_k)+z2)+z1;
            z0=h*(-lambda0*L^(1/5)*abs(z0-f_k)^(4/5)*sign(z0-f_k)+z1)+z0;  
        else
            z4=bk/h+z4;
            z3=h*(-lambda3*L^(4/5)*abs(z0-f_k)^(1/5)*sign(z0-f_k)+z4)+z3;
            z2=h*(-lambda2*L^(3/5)*abs(z0-f_k)^(2/5)*sign(z0-f_k)+z3)+z2;
            z1=h*(-lambda1*L^(2/5)*abs(z0-f_k)^(3/5)*sign(z0-f_k)+z2)+z1;
            z0=h*(-lambda0*L^(1/5)*abs(z0-f_k)^(4/5)*sign(z0-f_k)+z1)+z0;     
               
        end

	elseif differentiator_order==3
        lambda0=3;lambda1=4.16;lambda2=3.06;lambda3=1.1;        
            bk=f_k  -z0-h^1*z1-h^2*z2-h^3*z3...
            +h*lambda0*L^((1)/(4))*sign(z0-f_k)*abs(z0-f_k)^((3)/(4))...
            +h^2*lambda1*L^((2)/(4))*sign(z0-f_k)*abs(z0-f_k)^((2)/(4))...
            +h^3*lambda2*L^((3)/(4))*sign(z0-f_k)*abs(z0-f_k)^((1)/(4));
        
        if bk<-h^(4)*lambda3*L
            z3=h*(-lambda3*L)+z3;
            z2=h*(-lambda2*L^(3/4)*abs(z0-f_k)^(1/4)*sign(z0-f_k)+z3)+z2;
            z1=h*(-lambda1*L^(2/4)*abs(z0-f_k)^(2/4)*sign(z0-f_k)+z2)+z1;
            z0=h*(-lambda0*L^(1/4)*abs(z0-f_k)^(3/4)*sign(z0-f_k)+z1)+z0;            
             
        elseif bk>h^(4)*lambda3*L
            z3=h*(lambda3*L)+z3;
            z2=h*(-lambda2*L^(3/4)*abs(z0-f_k)^(1/4)*sign(z0-f_k)+z3)+z2;
            z1=h*(-lambda1*L^(2/4)*abs(z0-f_k)^(2/4)*sign(z0-f_k)+z2)+z1;
            z0=h*(-lambda0*L^(1/4)*abs(z0-f_k)^(3/4)*sign(z0-f_k)+z1)+z0;  
        else
            z3=bk/h+z3;
            z2=h*(-lambda2*L^(3/4)*abs(z0-f_k)^(1/4)*sign(z0-f_k)+z3)+z2;
            z1=h*(-lambda1*L^(2/4)*abs(z0-f_k)^(2/4)*sign(z0-f_k)+z2)+z1;
            z0=h*(-lambda0*L^(1/4)*abs(z0-f_k)^(3/4)*sign(z0-f_k)+z1)+z0;  
               
        end
        
	elseif differentiator_order==2
        lambda0=2;lambda1=2.12;lambda2=1.1;     
        bk=f_k  -z0-h^1*z1-h^2*z2...
            +h*lambda0*L^((1)/(3))*sign(z0-f_k)*abs(z0-f_k)^((2)/(3))...
            +h^2*lambda1*L^((2)/(3))*sign(z0-f_k)*abs(z0-f_k)^((1)/(3));
        
        if bk<-h^(3)*lambda2*L
            z2=h*(-lambda2*L)+z2;
            z1=h*(-lambda1*L^(2/3)*abs(z0-f_k)^(1/3)*sign(z0-f_k)+z2)+z1;
            z0=h*(-lambda0*L^(1/3)*abs(z0-f_k)^(2/3)*sign(z0-f_k)+z1)+z0;            
             
        elseif bk>h^(3)*lambda2*L
            z2=h*(lambda2*L)+z2;
            z1=h*(-lambda1*L^(2/3)*abs(z0-f_k)^(1/3)*sign(z0-f_k)+z2)+z1;
            z0=h*(-lambda0*L^(1/3)*abs(z0-f_k)^(2/3)*sign(z0-f_k)+z1)+z0;  
        else
            z2=bk/h+z2;
            z1=h*(-lambda1*L^(2/3)*abs(z0-f_k)^(1/3)*sign(z0-f_k)+z2)+z1;
            z0=h*(-lambda0*L^(1/3)*abs(z0-f_k)^(2/3)*sign(z0-f_k)+z1)+z0;     
               
        end    
    
    else
        disp("Error in E-AO-STD")
	end
        

    CT=CT+toc;
    
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
    
	calculation_time=CT;    
end

