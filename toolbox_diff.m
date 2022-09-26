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

function [u,calculation_time] = toolbox_diff(f,r,Ts, n,reset_flag,x0) %#codegen
% Calculates the derivatives up to n-th order of a given signal f
%   f... signal vector
%   n... maximum order of derivative
%   Ts... discretization time
%   r... robustness factor
% 
% outputs:
% e ... error signal (help for tuning r), should converge towards 0
% x ... containing the derivatives (n elements)

persistent Phi xp CT;
tic

if isempty(xp) || reset_flag==1
    xp=[f; zeros(10,1)];
    Phi = expm(diag(Ts*ones(10,1),1));
end

if isempty(CT) || reset_flag==1
    CT=0;
end

e=f-xp(1);

s = -r*(abs(e))^(-1/(n+1));
s0 = s;

z = exp(Ts*s);
z0 = exp(Ts*s0);

s0=s;
z0=z;

switch n
    case 1
        lambda =  ...
        [
                    2 - z0 - z
         ((z - 1)*(z0 - 1))/Ts
                             0
                             0
                             0
                             0
                             0
                             0
                             0
                             0
                             0
        ];
    case 2
        lambda =  ...
        [
                                   3 - z0 - 2*z
         ((z - 1)*(z + 3*z0 + z*z0 - 5))/(2*Ts)
                     -((z - 1)^2*(z0 - 1))/Ts^2
                                              0
                                              0
                                              0
                                              0
                                              0
                                              0
                                              0
                                              0
                                              
        ];
    case 3
        lambda =  ...
        [
                                                          4 - z0 - 3*z
         ((z - 1)*(7*z + 11*z0 + 5*z*z0 + 2*z^2*z0 + z^2 - 26))/(6*Ts)
                                   -((z - 1)^2*(2*z0 + z*z0 - 3))/Ts^2
                                             ((z - 1)^3*(z0 - 1))/Ts^3
                                                                     0
                                                                     0
                                                                     0
                                                                     0
                                                                     0
                                                                     0
                                                                     0
        ];
    case 4
        lambda = ...
        [
                                                                                5 - z0 - 4*z
         ((z - 1)*(23*z + 25*z0 + 13*z*z0 + 7*z^2*z0 + 3*z^3*z0 + 5*z^2 + z^3 - 77))/(12*Ts)
                       -((z - 1)^2*(35*z0 - 2*z + 26*z*z0 + 11*z^2*z0 + z^2 - 71))/(12*Ts^2)
                                               -((z - 1)^3*(z - 5*z0 - 3*z*z0 + 7))/(2*Ts^3)
                                                                  -((z - 1)^4*(z0 - 1))/Ts^4
                                                                                           0
                                                                                           0
                                                                                           0
                                                                                           0
                                                                                           0
                                                                                           0
        ];
    case 5
        lambda = ...
        [
                                                                                                             6 - z0 - 5*z
         ((z - 1)*(163*z + 137*z0 + 77*z*z0 + 47*z^2*z0 + 27*z^3*z0 + 12*z^4*z0 + 43*z^2 + 13*z^3 + 3*z^4 - 522))/(60*Ts)
                               -((z - 1)^2*(45*z0 - 7*z + 40*z*z0 + 25*z^2*z0 + 10*z^3*z0 + 2*z^2 + z^3 - 116))/(12*Ts^2)
                                                      -((z - 1)^3*(8*z - 17*z0 - 16*z*z0 - 7*z^2*z0 + z^2 + 31))/(4*Ts^3)
                                                                                 ((z - 1)^4*(z - 3*z0 - 2*z*z0 + 4))/Ts^4
                                                                                                ((z - 1)^5*(z0 - 1))/Ts^5
                                                                                                                        0
                                                                                                                        0
                                                                                                                        0
                                                                                                                        0
                                                                                                                        0
        ];
    case 6
        lambda = ...
        [
                                                                                                                                 7 - z0 - 6*z
         ((z - 1)*(213*z + 147*z0 + 87*z*z0 + 57*z^2*z0 + 37*z^3*z0 + 22*z^4*z0 + 10*z^5*z0 + 63*z^2 + 23*z^3 + 8*z^4 + 2*z^5 - 669))/(60*Ts)
                 -((z - 1)^2*(812*z0 - 232*z + 802*z*z0 + 597*z^2*z0 + 352*z^3*z0 + 137*z^4*z0 + 33*z^2 + 38*z^3 + 13*z^4 - 2552))/(180*Ts^2)
                                                   -((z - 1)^3*(39*z - 49*z0 - 57*z*z0 - 39*z^2*z0 - 15*z^3*z0 + 9*z^2 + z^3 + 111))/(8*Ts^3)
                                                                       ((z - 1)^4*(26*z - 35*z0 - 38*z*z0 - 17*z^2*z0 + 5*z^2 + 59))/(6*Ts^4)
                                                                                              -((z - 1)^5*(3*z - 7*z0 - 5*z*z0 + 9))/(2*Ts^5)
                                                                                                                   -((z - 1)^6*(z0 - 1))/Ts^6
                                                                                                                                            0
                                                                                                                                            0
                                                                                                                                            0
                                                                                                                                            0
        ];
    case 7
        lambda = ...
        [
                                                                                                                                                                    8 - z0 - 7*z
         ((z - 1)*(1851*z + 1089*z0 + 669*z*z0 + 459*z^2*z0 + 319*z^3*z0 + 214*z^4*z0 + 130*z^5*z0 + 60*z^6*z0 + 591*z^2 + 241*z^3 + 101*z^4 + 38*z^5 + 10*z^6 - 5772))/(420*Ts)
                              -((z - 1)^2*(938*z0 - 414*z + 994*z*z0 + 819*z^2*z0 + 574*z^3*z0 + 329*z^4*z0 + 126*z^5*z0 + 16*z^2 + 61*z^3 + 36*z^4 + 11*z^5 - 3490))/(180*Ts^2)
                                                 -((z - 1)^3*(1127*z - 967*z0 - 1277*z*z0 - 1077*z^2*z0 - 647*z^3*z0 - 232*z^4*z0 + 357*z^2 + 77*z^3 + 7*z^4 + 2632))/(120*Ts^3)
                                                                                    ((z - 1)^4*(68*z - 56*z0 - 77*z*z0 - 56*z^2*z0 - 21*z^3*z0 + 23*z^2 + 4*z^3 + 115))/(6*Ts^4)
                                                                                                        -((z - 1)^5*(43*z - 46*z0 - 55*z*z0 - 25*z^2*z0 + 10*z^2 + 73))/(6*Ts^5)
                                                                                                                                      ((z - 1)^6*(2*z - 4*z0 - 3*z*z0 + 5))/Ts^6
                                                                                                                                                       ((z - 1)^7*(z0 - 1))/Ts^7
                                                                                                                                                                               0
                                                                                                                                                                               0
                                                                                                                                                                               0
        ];
    case 8
        lambda = ...
        [
                                                                                                                                                                                                9 - z0 - 8*z
         ((z - 1)*(4437*z + 2283*z0 + 1443*z*z0 + 1023*z^2*z0 + 743*z^3*z0 + 533*z^4*z0 + 365*z^5*z0 + 225*z^6*z0 + 105*z^7*z0 + 1497*z^2 + 657*z^3 + 307*z^4 + 139*z^5 + 55*z^6 + 15*z^7 - 13827))/(840*Ts)
            -((z - 1)^2*(29531*z0 - 18254*z + 32926*z*z0 + 29013*z^2*z0 + 22468*z^3*z0 + 15293*z^4*z0 + 8622*z^5*z0 + 3267*z^6*z0 - 733*z^2 + 2172*z^3 + 1787*z^4 + 898*z^5 + 261*z^6 - 127251))/(5040*Ts^2)
                                                  -((z - 1)^3*(3777*z - 2403*z0 - 3457*z*z0 - 3302*z^2*z0 - 2442*z^3*z0 - 1367*z^4*z0 - 469*z^5*z0 + 1462*z^2 + 442*z^3 + 87*z^4 + 5*z^5 + 7667))/(240*Ts^3)
                                                                        ((z - 1)^4*(5572*z - 3207*z0 - 5092*z*z0 - 4682*z^2*z0 - 2852*z^3*z0 - 967*z^4*z0 + 2522*z^2 + 772*z^3 + 127*z^4 + 7807))/(240*Ts^4)
                                                                                                            -((z - 1)^5*(122*z - 81*z0 - 125*z*z0 - 95*z^2*z0 - 35*z^3*z0 + 50*z^2 + 10*z^3 + 154))/(6*Ts^5)
                                                                                                                                     ((z - 1)^6*(42*z - 39*z0 - 50*z*z0 - 23*z^2*z0 + 11*z^2 + 59))/(4*Ts^6)
                                                                                                                                                            -((z - 1)^7*(5*z - 9*z0 - 7*z*z0 + 11))/(2*Ts^7)
                                                                                                                                                                                  -((z - 1)^8*(z0 - 1))/Ts^8
                                                                                                                                                                                                           0
                                                                                                                                                                                                           0
        ];
    case 9
        lambda = ...
        [
                                                                                                                                                                                                                              10 - z0 - 9*z
         ((z - 1)*(15551*z + 7129*z0 + 4609*z*z0 + 3349*z^2*z0 + 2509*z^3*z0 + 1879*z^4*z0 + 1375*z^5*z0 + 955*z^6*z0 + 595*z^7*z0 + 280*z^8*z0 + 5471*z^2 + 2531*z^3 + 1271*z^4 + 641*z^5 + 305*z^6 + 125*z^7 + 35*z^8 - 48610))/(2520*Ts)
                -((z - 1)^2*(32575*z0 - 26477*z + 37754*z*z0 + 34905*z^2*z0 + 28864*z^3*z0 + 21689*z^4*z0 + 14514*z^5*z0 + 8095*z^6*z0 + 3044*z^7*z0 - 2712*z^2 + 2321*z^3 + 2566*z^4 + 1677*z^5 + 788*z^6 + 223*z^7 - 159826))/(5040*Ts^2)
                              ((z - 1)^3*(180920*z0 - 363543*z + 276981*z*z0 + 287607*z^2*z0 + 240602*z^3*z0 + 165702*z^4*z0 + 88737*z^5*z0 + 29531*z^6*z0 - 161922*z^2 - 60422*z^3 - 17337*z^4 - 2931*z^5 + 16*z^6 - 663941))/(15120*Ts^3)
                                                                            ((z - 1)^4*(9853*z - 4275*z0 - 7488*z*z0 - 7938*z^2*z0 - 6108*z^3*z0 - 3363*z^4*z0 - 1068*z^5*z0 + 5368*z^2 + 2198*z^3 + 638*z^4 + 101*z^5 + 12082))/(240*Ts^4)
                                                                                                    -((z - 1)^5*(6428*z - 3013*z0 - 5444*z*z0 - 5334*z^2*z0 - 3284*z^3*z0 - 1069*z^4*z0 + 3534*z^2 + 1244*z^3 + 229*z^4 + 6709))/(144*Ts^5)
                                                                                                                                            ((z - 1)^6*(129*z - 75*z0 - 126*z*z0 - 99*z^2*z0 - 36*z^3*z0 + 60*z^2 + 13*z^3 + 134))/(4*Ts^6)
                                                                                                                                                              -((z - 1)^7*(172*z - 145*z0 - 196*z*z0 - 91*z^2*z0 + 49*z^2 + 211))/(12*Ts^7)
                                                                                                                                                                                                 ((z - 1)^8*(3*z - 5*z0 - 4*z*z0 + 6))/Ts^8
                                                                                                                                                                                                                  ((z - 1)^9*(z0 - 1))/Ts^9
                                                                                                                                                                                                                                          0
        ];
    case 10
        lambda = ...
        [
                                                                                                                                                                                                                                                               11 - z0 - 10*z
                    ((z - 1)*(17819*z + 7381*z0 + 4861*z*z0 + 3601*z^2*z0 + 2761*z^3*z0 + 2131*z^4*z0 + 1627*z^5*z0 + 1207*z^6*z0 + 847*z^7*z0 + 532*z^8*z0 + 252*z^9*z0 + 6479*z^2 + 3119*z^3 + 1649*z^4 + 893*z^5 + 473*z^6 + 233*z^7 + 98*z^8 + 28*z^9 - 55991))/(2520*Ts)
         -((z - 1)^2*(177133*z0 - 181196*z + 211686*z*z0 + 202949*z^2*z0 + 175852*z^3*z0 + 140985*z^4*z0 + 104102*z^5*z0 + 68899*z^6*z0 + 38136*z^7*z0 + 14258*z^8*z0 - 27739*z^2 + 10278*z^3 + 16165*z^4 + 12728*z^5 + 7611*z^6 + 3454*z^7 + 962*z^8 - 976263))/(25200*Ts^2)
                               ((z - 1)^3*(420475*z0 - 1040321*z + 675075*z*z0 + 744585*z^2*z0 + 676405*z^3*z0 + 526605*z^4*z0 + 346845*z^5*z0 + 180175*z^6*z0 + 58635*z^7*z0 - 514467*z^2 - 222035*z^3 - 80075*z^4 - 21303*z^5 - 2669*z^6 + 427*z^7 - 1748357))/(30240*Ts^3)
                                                         ((z - 1)^4*(994506*z - 341693*z0 - 643092*z*z0 - 750990*z^2*z0 - 665660*z^3*z0 - 462765*z^4*z0 - 238632*z^5*z0 - 72368*z^6*z0 + 617430*z^2 + 304040*z^3 + 118155*z^4 + 33126*z^5 + 5084*z^6 + 1102859))/(15120*Ts^4)
                                                                                                       -((z - 1)^5*(24135*z - 8591*z0 - 17305*z*z0 - 19830*z^2*z0 - 15730*z^3*z0 - 8555*z^4*z0 - 2565*z^5*z0 + 16010*z^2 + 7550*z^3 + 2445*z^4 + 427*z^5 + 22009))/(288*Ts^5)
                                                                                                                                  ((z - 1)^6*(18188*z - 7513*z0 - 14948*z*z0 - 15378*z^2*z0 - 9548*z^3*z0 - 3013*z^4*z0 + 11418*z^2 + 4388*z^3 + 853*z^4 + 15553))/(240*Ts^6)
                                                                                                                                                                    -((z - 1)^7*(1139*z - 605*z0 - 1085*z*z0 - 875*z^2*z0 - 315*z^3*z0 + 581*z^2 + 133*z^3 + 1027))/(24*Ts^7)
                                                                                                                                                                                                      ((z - 1)^8*(56*z - 44*z0 - 62*z*z0 - 29*z^2*z0 + 17*z^2 + 62))/(3*Ts^8)
                                                                                                                                                                                                                            -((z - 1)^9*(7*z - 11*z0 - 9*z*z0 + 13))/(2*Ts^9)
                                                                                                                                                                                                                                                 -((z - 1)^10*(z0 - 1))/Ts^10
        ];
    otherwise
        lambda = zeros(11,1);
end
x=xp(2:end);

xp = Phi*xp + lambda*e;

CT=toc+CT;
calculation_time=CT;
u=x(n);

end