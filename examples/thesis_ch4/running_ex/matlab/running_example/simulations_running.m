model_path = fileparts(which('simulations_running.m'));
addpath(genpath(model_path));

clear
clc
close all
warning('off','all')
modelDefName = 'modelDef_running';
t = linspace(0,10,1000);
theta = [1.0;0.5;0.2;100.0]

method = 'CMEC_2_ZC_2_a';
modelName = 'running_MCM';
tic
System_running = genSimFileIDA(modelName, modelDefName, method)
toc

tic
for i = 1:100
    System_running.sol = simulate_running_MCM(t,theta);
end
toc

save('../running_MCM.mat', 'System_running');
