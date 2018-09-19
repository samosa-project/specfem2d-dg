% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

clear all;
% close all;
% clc;
format compact;
set(0, 'DefaultLineLineWidth', 2); set(0, 'DefaultLineMarkerSize', 8);
set(0, 'defaultTextFontSize', 12); set(0, 'defaultAxesFontSize', 12);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

% FILE='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200/OUTPUT_FILES_597316/Database00201';
disp(['1) Find the ID of CPU which crashed. Typically, look up the error message, and remember the number of the task which is mentionned. 2) Get the database file corresponding to that CPU.']);
disp(['For instance, if the error message is ''srun: error: eoscomp0: task 12: Floating point exception'', the CPU ID is 12, and you should get the Database file number 12.']);
FILE=input('Database file to scan? > ','s');

X=scan_database_file(FILE);

disp(strcat('x_min=', sprintf("%11.3e", min(X(:,1))),', x_max=',sprintf("%11.3e", max(X(:,1)))));
disp(strcat('z_min=', sprintf("%11.3e", min(X(:,2))),', z_max=',sprintf("%11.3e", max(X(:,2)))));

figure();
rectangle('Position', [min(X(:,1)), min(X(:,2)), max(X(:,1))-min(X(:,1)), max(X(:,2))-min(X(:,2))]);