clc
clear all
close all

%% INITIALIZATION

load('exp_data.mat')

%% SPARSEPOP MIN

thetaMIN=zeros(5,1);
    
for k=1:5
    
    %POLYNOMIAL CONSTRAIN

    SUPPpoly=zeros(1,55);
    SUPPpoly(k)=1;
    
    objPoly.typeCone=1;
    objPoly.dimVar=55;
    objPoly.degree=1;
    objPoly.noTerms=1;
    objPoly.supports=SUPPpoly;
    objPoly.coef=[1];
    
    % INEQUALITY CONSTRAIN

    SUPPineq = zeros(9,55);
    SUPPineq(:,1:5)=   [0 0 0 0 0
                        1 0 0 0 0 
                        0 1 0 0 0
                        0 1 0 0 0
                        1 0 0 0 0
                        0 0 0 0 0
                        0 0 1 0 0
                        0 0 0 1 0
                        0 0 0 0 1];
    
    for i= 3:50
    
        SUPPineq(4,3+i) = 1;
        SUPPineq(5,4+i) = 1;
        SUPPineq(6,5+i) = 1;
    
        ineqPolySys{i-2}.typeCone=-1;
        ineqPolySys{i-2}.dimVar=55;
        ineqPolySys{i-2}.degree=2;
        ineqPolySys{i-2}.noTerms=9;
        ineqPolySys{i-2}.supports=SUPPineq;
        ineqPolySys{i-2}.coef=[y_tilde(i); y_tilde(i-1); y_tilde(i-2); -1; -1; -1; -u(i); -u(i-1); -u(i-2)];
    
        SUPPineq(4,3+i) = 0;
        SUPPineq(5,4+i) = 0;
        SUPPineq(6,5+i) = 0;
    end;
    
    %BOUNDED VALUE
    ubd=zeros(1,55);
    ubd(1:5)=1e10;
    ubd(6:55)=5;
    
    lbd=zeros(1,55);
    lbd(1:5)=-1e10;
    lbd(6:55)=-5;
    
    %PARAM
    param.relaxOrder=1;
    param.POPsolver="active-set";
    
    %RUN
    [param, SDPobjValue, POP, elapsedTime, SDPsolverInfo, SDPinfo]=sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);

    thetaMIN(k)=POP.objValueL;
    SUPPpoly(k)=0;

end

%% SPARSEPOP MAX

thetaMAX=zeros(5,1);

for k=1:5

    %POLYNOMIAL CONSTRAIN

    SUPPpoly=zeros(1,55);
    SUPPpoly(k)=1;

    objPoly.typeCone=1;
    objPoly.dimVar=55;
    objPoly.degree=1;
    objPoly.noTerms=1;
    objPoly.supports=SUPPpoly;
    objPoly.coef=[-1];

    % INEQUALITY CONSTRAIN

    SUPPineq = zeros(9,55);
    SUPPineq(:,1:5)=   [0 0 0 0 0
                        1 0 0 0 0 
                        0 1 0 0 0
                        0 1 0 0 0
                        1 0 0 0 0
                        0 0 0 0 0
                        0 0 1 0 0
                        0 0 0 1 0
                        0 0 0 0 1];

    for i= 3:50

        SUPPineq(4,3+i) = 1;
        SUPPineq(5,4+i) = 1;
        SUPPineq(6,5+i) = 1;

        ineqPolySys{i-2}.typeCone=-1;
        ineqPolySys{i-2}.dimVar=55;
        ineqPolySys{i-2}.degree=2;
        ineqPolySys{i-2}.noTerms=9;
        ineqPolySys{i-2}.supports=SUPPineq;
        ineqPolySys{i-2}.coef=[y_tilde(i); y_tilde(i-1); y_tilde(i-2); -1; -1; -1; -u(i); -u(i-1); -u(i-2)];

        SUPPineq(4,3+i) = 0;
        SUPPineq(5,4+i) = 0;
        SUPPineq(6,5+i) = 0;
    end;

    %BOUNDED VALUE
    ubd=zeros(1,55);
    ubd(1:5)=1e10;
    ubd(6:55)=5;

    lbd=zeros(1,55);
    lbd(1:5)=-1e10;
    lbd(6:55)=-5;

    %PARAM
    param.relaxOrder=1;
    param.POPsolver="active-set";

    %RUN
    [param, SDPobjValue, POP, elapsedTime, SDPsolverInfo, SDPinfo]=sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);
    thetaMAX(k)=-POP.objValueL;
    SUPPpoly(k)=0;
end

%% PUIs mean

thetaC=zeros(5,1);

for j=1:5;

    thetaC(j,1)=(thetaMIN(j)+thetaMAX(j))/2;

end
thetaC