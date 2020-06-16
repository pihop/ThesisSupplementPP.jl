%% MODEL
% Definition of symbolic variables:
syms x01 x00 x10 x11 z1 z2 %  species
syms y01 y00 y10 y11 m1 m2 %  species 
syms rs rm rm2 % kinetic rates
syms r0 % initial count
syms time % time

% Defininition of general fileds:
System.time = time;
System.compartments = {'cell'};
System.volumes = [1];

% Defininition of state fileds:
System.state.variable = [z1;z2;x01;x00;x10;x11;];
System.state.name = {'z_{1}';'z_{2}';'x_{01}';'x_{00}';'x_{10}';'x_{11}';};
System.state.compartment = {'cell';'cell';'cell';'cell';'cell'; 'cell'};
System.state.type     = {'stochastic';'stochastic';'moment';'moment';'moment';'moment';};
System.state.xmin     = [0;0;0;0;0;0];
System.state.xmax     = [1;1;1000;1000;1000;1000;];
System.state.mu0      = [1;0;0;r0;0;0;];
System.state.C0       = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1)+x(2)) == 1);

% Definition of parameters field:
System.parameter.variable = [rm;rm2;rs;r0];
% Definition of constant parameters filed:

% Specifying the scale of propensities and parameters
System.scaleIndicator = 'microscopic'; % options are 'microscopic' and 'macroscopic'

% Define propensities:
% Mode 1
% (R1)
System.reaction(1).educt      = [x01,z1];
System.reaction(1).product    = [x00,z1];
%System.reaction(1).propensity = rm*x01*z1;
System.reaction(1).propensity = rm*x01*z1;

% (R2)
System.reaction(2).educt      = [x00,z1];
System.reaction(2).product    = [x01,z1];
%System.reaction(2).propensity = 0.5*rm*x00*z1;
System.reaction(2).propensity = rm2*x00*z1;

% (R3)
System.reaction(3).educt      = [x00,z1];
System.reaction(3).product    = [x10,z1];
%System.reaction(3).propensity = 0.5*rm*x00*z1;
System.reaction(3).propensity = rm2*x00*z1;

% (R4)
System.reaction(4).educt      = [x10,z1];
System.reaction(4).product    = [x11,z1];
%System.reaction(4).propensity = 0.5*rm*x10*z1;
System.reaction(4).propensity = rm2*x10*z1;

% (R5)
System.reaction(5).educt      = [x10,z1];
System.reaction(5).product    = [x00,z1];
%System.reaction(5).propensity = 0.5*rm*x10*z1;
System.reaction(5).propensity = rm2*x10*z1;

System.reaction(10).educt      = [x11,z1];
System.reaction(10).product    = [x10,z1];
%System.reaction(5).propensity = 0.5*rm*x10*z1;
System.reaction(10).propensity = rm*x11*z1;


% (R6)
System.reaction(6).educt      = [x11,z1];
System.reaction(6).product    = [x11,z2];
%System.reaction(6).propensity = rs*x11*z1;
System.reaction(6).propensity = rs*x11*z1;

% Mode 2
% (R1)
System.reaction(7).educt      = [x01,z2];
System.reaction(7).product    = [x00,z2];
%System.reaction(7).propensity = rm*x01*z2;
System.reaction(7).propensity = rm*x01*z2;

% (R3)
System.reaction(8).educt      = [x00,z2];
System.reaction(8).product    = [x10,z2];
%System.reaction(8).propensity = rm*x00*z2;
System.reaction(8).propensity = rm*x00*z2;

% (R4)
System.reaction(9).educt      = [x10,z2];
System.reaction(9).product    = [x11,z2];
%System.reaction(9).propensity = rm*x10*z2;
System.reaction(9).propensity = rm*x10*z2;

System.output.variable = [ y01; y00; y10; y11;];
System.output.function = [ x01; x00; x10; x11;];
