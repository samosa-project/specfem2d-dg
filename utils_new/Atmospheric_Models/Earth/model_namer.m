% Author:        Léo Martire.
% Description:   TODO.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function comprehensive_name = model_namer(altMin, altMax, nLayers, lat, lon, yyears, ddays, ssecs, wProjAng)
  comprehensive_name=strcat(num2str(2000+yyears),"_",num2str(ddays),"_",sprintf("%.d",ssecs),"_",sprintf("%.5f",lat),"_",sprintf("%.5f",lon),"_",num2str(altMin),"_",num2str(altMax),"_",num2str(nLayers),"_",sprintf("%.5f",wProjAng));
  comprehensive_name=char(comprehensive_name);
end

