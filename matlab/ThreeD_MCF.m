function []=ThreeD_MCF()
%=======================================================================
%   A New 3-D Mininum Cost Function Phase Unwrapping Alogrithm Based on
%   Closure Phase
%
%   This is the main function of 3-D MCF PU algorithm.This function should 
%   be called after editing the 'InputParameters.m' file located in your 
%   work directory.
%   
%   Written by Fei Liu (liufei.whu@gmail.com)
%   Version 1.0
%   Release date: Oct 30, 2019
%=======================================================================

%   Set the input parameters
InputParameters;

%   Load and check input data
disp('Loading input data and obtaining closure phases...');
[ph,weight,closure,datatype]=LoadData(INPUTFILE,DATATYPE,WRAPPEDPHASES, ...
    WEIGHTCOEFFICIENTS,INTERFEROGRAMS,COORDINATE);
ph=single(ph/(2*pi));

%   Calculate the sum of unwrapped phase gradients (UPGs) in closure phase
disp('Calculating the sum of UPGs in closure phase...');
[dk_s]=CalSumUPGs(ph,closure,datatype);

%   Calculate UPGs and acquire residuals in each interferogram
disp('Calculating the UPGs and acquiring residuals in each interferogram...');
[res,dk]=CalRes(ph,datatype);

%	Save these temporary data for PU
disp('Saving temporary data...');
save('Tmpdata','res','dk','dk_s','ph','weight','datatype');
clear;

%   Unwrap all the phases jointly
disp('Unwrapping all the interferograms jointly...');
JointUnwrap(); 
disp('Finished.');

end
