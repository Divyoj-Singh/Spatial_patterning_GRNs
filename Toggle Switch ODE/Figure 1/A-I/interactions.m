%% Code which contains the equations:
function out = interactions
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,gammaA,gammaB,gA,gB,A0B,A0A,lB0A,lA0A,nB0A,nA0A,B0A,B0B,lA0B,lB0B,nA0B,nB0B)
%% The equations are entered here:
% the order of parameters given to fun_eval should be conserved throughout
% the codes. kmrgd contains the species, so kmrgd(1) = A. The equations are
% seprated via ";" and there is no ; after the last equation. Instead you
% should put it after closing the square bracket to supress diplaying the
% values on screen. For example, for ORDINARY TOGGLE SWITCH:
% dydt=[gA*hill(kmrgd(2),B0A,lB0A,nB0A)-gammaA*kmrgd(1);
%        gB*hill(kmrgd(1),A0B,lA0B,nA0B)-gammaB*kmrgd(2)];
%FOR SA_TOGGLE SWITCH :
% dydt=[gA*hill(kmrgd(2),B0A,lB0A,nB0A)*hill(kmrgd(1),A0A,lA0A,nA0A)-gammaA*kmrgd(1);
%        gB*hill(kmrgd(1),A0B,lA0B,nA0B)*hill(kmrgd(2),B0B,lB0B,nB0B)-gammaB*kmrgd(2)];
%For TOGGLE TRIAD WITHOUT ACTIVATION :
%   dydt=[ g_A*hill(kmrgd(2),B0A,lambda_BtoA,nBtoA)*hill(kmrgd(3),C0A,lambda_CtoA,nCtoA) - gamma_A*kmrgd(1); 
%          g_B*hill(kmrgd(1),A0B,lambda_AtoB,nAtoB)*hill(kmrgd(3),C0B,lambda_CtoB,nCtoB) - gamma_B*kmrgd(2); 
%          g_C*hill(kmrgd(1),A0C,lambda_AtoC,nAtoC)*hill(kmrgd(2),B0C,lambda_BtoC,nBtoC) - gamma_C*kmrgd(3)];

dydt=[gA*hill(kmrgd(2),B0A,lB0A,nB0A)-gammaA*kmrgd(1);
       gB*hill(kmrgd(1),A0B,lA0B,nA0B)-gammaB*kmrgd(2)];

    
% --------------------------------------------------------------------------
function [tspan,y0,options] = init
%% Let this be as it is:
handles = feval(interactions);
y0=[0,0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

%% This part is also not touched (didn't change the results as of now)
% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,s,gammaA,gammaB,gA,gB,B0,lB0A,nB0A,A0,lA0B,nA0B)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,s,gammaA,gammaB,gA,gB,B0,lB0A,nB0A,A0,lA0B,nA0B)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,s,gammaA,gammaB,gA,gB,B0,lB0A,nB0A,A0,lA0B,nA0B)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,s,gammaA,gammaB,gA,gB,B0,lB0A,nB0A,A0,lA0B,nA0B)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,s,gammaA,gammaB,gA,gB,B0,lB0A,nB0A,A0,lA0B,nA0B)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,s,gammaA,gammaB,gA,gB,B0,lB0A,nB0A,A0,lA0B,nA0B)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,s,gammaA,gammaB,gA,gB,B0,lB0A,nB0A,A0,lA0B,nA0B)
