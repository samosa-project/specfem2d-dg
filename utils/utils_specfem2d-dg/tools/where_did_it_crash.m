% Author:        Léo Martire.
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

clear all;
% close all;
% clc;
format compact;

% FILE='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/ON_EOS_STRATO_SAVE/stratoexplo_66_june_1200/OUTPUT_FILES_597316/Database00201';
disp(['[',mfilename,'] 1) Find the ID of CPU which crashed. Typically, look up the error message, and remember the number of the task which is mentionned.']);
disp(['[',mfilename,'] 2) Get the database file corresponding to that CPU.']);
disp(['[',mfilename,'] For instance, if the error message is ''srun: error: eoscomp0: task 12: Floating point exception'', the CPU ID is 12, and you should get the Database file number 12.']);
FILE=input(['[',mfilename,'] Path to database file to scan? > '],'s');

if(not(exist(FILE,'file')))
  error(['[, ERROR] File ''',FILE,''' does not exist.']);
end

[FOLDER,NAME,~] = fileparts(FILE);
CPU = regexp(NAME,'[0-9]+','match'); CPU = str2num(CPU{1});
boundary_only = 1;
plot_partitions(1, FOLDER, boundary_only, CPU);

% X=scan_database_file(FILE);
% disp(strcat('x_min=', sprintf("%11.3e", min(X(:,1))),', x_max=',sprintf("%11.3e", max(X(:,1)))));
% disp(strcat('z_min=', sprintf("%11.3e", min(X(:,2))),', z_max=',sprintf("%11.3e", max(X(:,2)))));
% figure();
% rectangle('Position', [min(X(:,1)), min(X(:,2)), max(X(:,1))-min(X(:,1)), max(X(:,2))-min(X(:,2))]);