% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   Produces a variety of filters.
% Last modified: See file metadata.
% Usage:         N/A.
% Notes:         N/A.

function [freq,hlp,hhp,hbp] = custom_Filter(Fsample,Nsample,fclp,fchp,dlp,dhp,order)
  freq=linspace(-Fsample/2, Fsample/2, Nsample);
  switch(order)
    case(1)
      hhp = 1i*freq/fchp./(1 +  1i*freq/fchp);
      hlp = 1./(1 +  1i*freq/fclp);
    case(2)
      hhp = -(freq/fchp).^2 ./ (1 +  2*dhp*1i*freq/fchp - (freq/fchp).^2);
      hlp = 1./(1 +  2*dlp*1i*freq/fclp - (freq/fclp).^2);
    case(4)
      hhp = ((freq/fchp).^2 ./ (1 +  2*dhp*1i*freq/fchp - (freq/fchp).^2)).^2;
      hlp = 1./(1 +  2*dlp*1i*freq/fclp - (freq/fclp).^2).^2;
    otherwise
      error(['[',mfilename,', ERROR] Filter order unknown.']);
  end
  hbp=hlp.*hhp;
end

