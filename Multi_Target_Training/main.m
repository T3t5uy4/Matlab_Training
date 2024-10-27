clc
clear
close all

% Used to add the current working directory and all its subdirectories to MATLAB's search path
addpath(genpath(pwd));

% Set initial training parameters
% initial population size
searchAgentsNum = 30;
numOfRecord = 50;
% population dimension
dim = 30;
% Maximum number of iterations, recommended maximum training is 300,000
maxFes = 3000;

% algorithm name
algorithmName = 'MOPSO';

% Fold is recommended to be set to 30
fold = 3;

% Prepare output preprocessing
% Get the current year, month, day, time and minute
dateStr = datestr(now, 'yyyy-mm-dd');
timeStr = datestr(now, 'HH_MM');
% Set output directory name and file name
algNameStr = strjoin(algorithmName, '_');
dirName = ['output/', dateStr, '_', algNameStr, '/', timeStr];
mkdir(dirName);
fileName = [dirName, '/', timeStr];
