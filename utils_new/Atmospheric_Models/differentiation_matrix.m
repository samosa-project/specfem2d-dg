% Author:        Quentin Brissaud, Léo Martire.
% Mail:          quentinb@gps.caltech.edu
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         As of 2018/09/18, only used by 'modify_atmos_model.m' and 'Richardson_number.m'.

function [D] = differentiation_matrix(Z, conditionning)
  % Ref.: http://acoustique.ec-lyon.fr/publi/haniquecockenpot_thesis.pdf.
  dz = Z(2) - Z(1);
  nz = length(Z);
  D = sparse(zeros(nz));
  signe = 1;
  i = 1;
  D(i, i) = - 2.391602219538;
  D(i, i + signe * 1) = 5.832490322294; D(i, i + signe * 2) = - 7.650218001182; D(i, i + signe * 3) = 7.907810563576; D(i, i + signe * 4) = - 5.922599052629; D(i, i + signe * 5) = 3.071037015445;
  D(i, i + signe * 6) = - 1.014956769726; D(i, i + signe * 7) = 0.170022256519; D(i, i + signe * 8) = 0.002819958377; D(i, i + signe * 9) = - 0.004791009708; D(i, i + signe * 10) = - 0.000013063429;
  i = 2;
  D(i, i - signe * 1) = - 0.180022054228;
  D(i, i) = - 1.237550583044;
  D(i, i + signe * 1) = 2.484731692990; D(i, i + signe * 2) = - 1.810320814061; D(i, i + signe * 3) = 1.112990048440; D(i, i + signe * 4) = - 0.481086916514; D(i, i + signe * 5) = 0.126598690230;
  D(i, i + signe * 6) = - 0.015510730165; D(i, i + signe * 7) = 0.000021609059; D(i, i + signe * 8) = 0.000156447571; D(i, i + signe * 9) = - 0.000007390277;
  i = 3;
  D(i, i - signe * 2) = 0.057982271137; D(i, i - signe * 1) = - 0.536135360383;
  D(i, i) = - 0.264089548967;
  D(i, i + signe * 1) = 0.917445877606; D(i, i + signe * 2) = - 0.169688364841; D(i, i + signe * 3) = - 0.029716326170; D(i, i + signe * 4) = 0.029681617641; D(i, i + signe * 5) = - 0.005222483773;
  D(i, i + signe * 6) = - 0.000118806260; D(i, i + signe * 7) = - 0.000118806260; D(i, i + signe * 8) = - 0.000020069730;
  i = 4;
  D(i, i - signe * 3) = - 0.013277273810; D(i, i - signe * 2) = 0.115976072920; D(i, i - signe * 1) = - 0.617479187931;
  D(i, i) = - 0.274113948206;
  D(i, i + signe * 1) = 1.086208764655; D(i, i + signe * 2) = - 0.402951626982; D(i, i + signe * 3) = 0.131066986242; D(i, i + signe * 4) = - 0.028154858354; D(i, i + signe * 5) = 0.002596328316;
  D(i, i + signe * 6) = 0.000128743150; D(i, i + signe * 7) = - 0.0;
  i = 5;
  D(i, i - signe * 4) = 0.016756572303; D(i, i - signe * 3) = - 0.117478455239; D(i, i - signe * 2) = 0.411034935097; D(i, i - signe * 1) = - 1.130286765151;
  D(i, i) = 0.341435872100;
  D(i, i + signe * 1) = 0.556396830543; D(i, i + signe * 2) = - 0.082525734207; D(i, i + signe * 3) = 0.003565834658; D(i, i + signe * 4) = 0.001173034777; D(i, i + signe * 5) = - 0.000071772671;
  D(i, i + signe * 6) = - 0.000000352273;
  for i = 6:nz - 5
    D(i, i - 1) = 0.872756993962667; D(i, i - 2) = - 0.286511173973333; D(i, i - 3) = 0.090320001280000; D(i, i - 4) = - 0.020779405824000; D(i, i - 5) = 0.002484594688000;
    D(i, :) = - D(i, :);
    D(i, i + 1) = 0.872756993962667; D(i, i + 2) = - 0.286511173973333; D(i, i + 3) = 0.090320001280000; D(i, i + 4) = - 0.020779405824000; D(i, i + 5) = 0.002484594688000;
  end
  signe = - 1;
  i = nz;
  D(i, i) = - 2.391602219538;
  D(i, i + signe * 1) = 5.832490322294; D(i, i + signe * 2) = - 7.650218001182; D(i, i + signe * 3) = 7.907810563576; D(i, i + signe * 4) = - 5.922599052629; D(i, i + signe * 5) = 3.071037015445;
  D(i, i + signe * 6) = - 1.014956769726; D(i, i + signe * 7) = 0.170022256519; D(i, i + signe * 8) = 0.002819958377; D(i, i + signe * 9) = - 0.004791009708; D(i, i + signe * 10) = - 0.000013063429;
  i = nz - 1;
  D(i, i - signe * 1) = - 0.180022054228;
  D(i, i) = - 1.237550583044;
  D(i, i + signe * 1) = 2.484731692990; D(i, i + signe * 2) = - 1.810320814061; D(i, i + signe * 3) = 1.112990048440; D(i, i + signe * 4) = - 0.481086916514; D(i, i + signe * 5) = 0.126598690230;
  D(i, i + signe * 6) = - 0.015510730165; D(i, i + signe * 7) = 0.000021609059; D(i, i + signe * 8) = 0.000156447571; D(i, i + signe * 9) = - 0.000007390277;
  i = nz - 2;
  D(i, i - signe * 2) = 0.057982271137; D(i, i - signe * 1) = - 0.536135360383;
  D(i, i) = - 0.264089548967;
  D(i, i + signe * 1) = 0.917445877606; D(i, i + signe * 2) = - 0.169688364841; D(i, i + signe * 3) = - 0.029716326170; D(i, i + signe * 4) = 0.029681617641; D(i, i + signe * 5) = - 0.005222483773;
  D(i, i + signe * 6) = - 0.000118806260; D(i, i + signe * 7) = - 0.000118806260; D(i, i + signe * 8) = - 0.000020069730;
  i = nz - 3;
  D(i, i - signe * 3) = - 0.013277273810; D(i, i - signe * 2) = 0.115976072920; D(i, i - signe * 1) = - 0.617479187931;
  D(i, i) = - 0.274113948206;
  D(i, i + signe * 1) = 1.086208764655; D(i, i + signe * 2) = - 0.402951626982; D(i, i + signe * 3) = 0.131066986242; D(i, i + signe * 4) = - 0.028154858354; D(i, i + signe * 5) = 0.002596328316;
  D(i, i + signe * 6) = 0.000128743150; D(i, i + signe * 7) = - 0.0;
  i = nz - 4;
  D(i, i - signe * 4) = 0.016756572303; D(i, i - signe * 3) = - 0.117478455239; D(i, i - signe * 2) = 0.411034935097; D(i, i - signe * 1) = - 1.130286765151;
  D(i, i) = 0.341435872100;
  D(i, i + signe * 1) = 0.556396830543; D(i, i + signe * 2) = - 0.082525734207; D(i, i + signe * 3) = 0.003565834658; D(i, i + signe * 4) = 0.001173034777; D(i, i + signe * 5) = - 0.000071772671;
  D(i, i + signe * 6) = - 0.000000352273;
  D(end-4:end, 1:end) = - D(end-4:end, 1:end);
  D = D / dz;
  
  if(conditionning==1)
    D(1, 1) = 1.0; D(1, 2:end) = 0.0; % Conditioning.
  end
  % D(end, end) = 1.0; D(end, 1:end-1) = 0.0; % Conditioning.
end