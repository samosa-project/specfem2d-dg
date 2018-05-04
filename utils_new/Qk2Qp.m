% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         Based on /doc/note_on_Qkappa_versus_Qp.pdf.

clear all;
close all;
% clc;

swi=-1;
while(not(ismember(swi,[0,1])))
  swi=input('  QpQs to QkappaQmu (0) or QkappaQmu to QpQs (1)? > ');
end

cp=input('  v_p? > ');
cs=input('  v_s? > ');
rcsq=(cs./cp).^2;

if(swi==0)
  qp=input('  q_p? > ');
  qs=input('  q_s? > ');
  qk=(1-rcsq).*(qp.^(-1) - rcsq.*qs.^(-1)).^(-1);
  qm=qs;
  disp(['q_\kappa = ', sprintf('%.16f ',qk)]);
  disp(['q_\mu    = ', sprintf('%.16f ',qm)]);
else
  qk=input('  q_\kappa? > ');
  qm=input('  q_\mu? > ');
  qp=((1-rcsq).*qk.^(-1) + rcsq.*qm.^(-1)).^(-1);
  qs=qm;
  disp(['q_p = ', sprintf('%.16f ',qp)]);
  disp(['q_s = ', sprintf('%.16f ',qs)]);
end