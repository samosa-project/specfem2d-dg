% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

% clear all;
% close all;
clc;

data=importdata("~/Downloads/595500_macrobloc302");

disp(strcat('x_min=', sprintf("%11.3e", min(data(:,2))),', x_max=',sprintf("%11.3e", max(data(:,2)))));
disp(strcat('z_min=', sprintf("%11.3e", min(data(:,3))),', z_max=',sprintf("%11.3e", max(data(:,3)))));