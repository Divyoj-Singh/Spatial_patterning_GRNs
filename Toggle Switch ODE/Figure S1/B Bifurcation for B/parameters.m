%% Code which contains the parameters:
function [x,v,s,h,f] = parameters
curdir = pwd;
init;
cd(curdir);

%% Control the fine tuning of the bifurcation here:
% you can increase or decrease the step size, total number of points etc
% here:
opt = contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',20000);
opt=contset(opt,'MinStepsize',0.001);
opt=contset(opt,'MaxStepsize',0.01);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'Backward',1);


%% Paramaters:
% Mention all the parameters here.
% Important thing to note here is that matcont will start drawing
% bifurcation from the value you mention here.
gammaA=0.5;  gammaB=2;  %Degradation rate
gA=10;          gB=10;  %Production rate


% B0=5;  lB0A=0.2;   nB0A=10; 

A0B=2; B0A=2; A0A=80 ; B0B=80; %Threshold
lB0A=0.1; lA0B=0.1 ; lA0A=5; lB0B=5;%Fold change
nB0A=5; nA0B=5; nA0A=4; nB0B=4; %Hill's coefficient
 
% A0=5;  lA0B=0.2;   nA0B=10;


%

ap = 2; 
% describes the index of parameter for which the bifurcation is drawn using 
% the init_EP_EP function (second last line of the code). 
% Currently, ap=1, thus bifurcation parameter is gammaA (degradation rate of A)

handles = feval(@interactions);
tspan = 0:1:1000;

% initial condition for the species concentration:
x_start = [5,200];
%FOR ORDINARY TOGGLE SWITCH
% Just copy the parameters: gammaA,gammaB,gA,gB,B0,lB0A,nB0A,A0,lA0B,nA0B
% (FOR SA_TOGGLE SWITCH ) use : gammaA,gammaB,gA,gB,A0B,A0A,lB0A,lA0A,nB0A,nA0A,B0A,B0B,lA0B,lB0B,nA0B,nB0B

% in the conserved order to functions below as well.

%calculating steady state for given initial condition 
[t,x_time] = ode15s(@(t,kmrgd)handles{2}(t, kmrgd, gammaA,gammaB,gA,gB,A0B,A0A,lB0A,lA0A,nB0A,nA0A,B0A,B0B,lA0B,lB0B,nA0B,nB0B),tspan,x_start);
x_init = x_time(end,:)'

%drawing bifurcation using a continuation method
[x0,v0] = init_EP_EP(@interactions,x_init,[gammaA,gammaB,gA,gB,A0B,A0A,lB0A,lA0A,nB0A,nA0A,B0A,B0B,lA0B,lB0B,nA0B,nB0B],ap);
[x,v,s,h,f] = cont(@equilibrium, x0, v0,opt);

end


