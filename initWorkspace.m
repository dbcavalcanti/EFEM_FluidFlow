function initWorkspace()
    close all;
    % Clear the classes to avoid using a non-updated version
    clear Element RegularElement  EnrichedElement;
    clear EnrichedElement_SOS EnrichedElement_KOS EnrichedElement_KSON;
    clear Fracture Fracture_ConstantJump Fracture_LinearJump;
    clear IntPoint;
    clear Shape Shape_CST Shape_LST Shape_ISOQ4 Shape_ISOQ8;
    clear Model IntPoint Result;
    clear Anl Anl_Linear Anl_Nonlinear;
    % Clear the workspace and the command window
    clear; clc;
    %Use all folders and subfolders
    addpath(genpath('./'));
end