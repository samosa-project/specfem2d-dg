% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         Based on /doc/note_on_Qkappa_versus_Qp.pdf.

cp=input('  v_p? > ');
cs=input('  v_s? > ');
qp=input('  q_p? > ');
qs=input('  q_s? > ');

rcsq=(cs/cp)^2;

qk=((1-rcsq)^(-1) * qp^(-1) - rcsq*qs^(-1))^(-1);
qm=qs;

disp(['q_\kappa = ', sprintf('%.16f',qk)]);
disp(['q_\mu = ', sprintf('%.16f',qm)]);