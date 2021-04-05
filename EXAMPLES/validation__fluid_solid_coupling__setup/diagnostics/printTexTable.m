function [] = printTexTable(ic_deg, cases, saveR, saveT, effective_R, effective_T, err_R, err_T, stftag, ftstag)
  nsamples = 1;
  ncoef = 3;
  
  siunitxavailable = 0;
  
  setup_find_all_cases; % produces istfortho, istfslant, iftsslant, iftsortho
  
  stfptp = '\stfptp';
  stfprp = '\stfprp';
  stfprs = '\stfprs';
  ftsptp = '\ftsptp';
  ftspts = '\ftspts';
  ftsprp = '\ftsprp';
  
  if(siunitxavailable)
    numcmd = '',numcmd,'';
  else
    numcmd = '';
  end
  
  slantedtxt = [numcmd,'{',num2str(ic_deg),'}'];
  
  mr = '\multirow';
  
  str = '';
  str = [str, '\begin{tabular}{c|c|c|ccc}']; str=andl(str);

  str = [str, '    \hline']; str=andl(str);
  str = [str, '    ',mr,'{2}{*}{Case} & ',mr,'{2}{*}{$\theta_{\mathrm{i}}$ [\si{\degree}]} & ',mr,'{2}{*}{Coefficient} & Theoretical & Synthetic & Relative\\']; str=andl(str);
  str = [str, '    & & & Value & Value & Error [\si{\percent}]\\']; str=andl(str);
  str = [str, '    \hline\hline']; str=andl(str);

  % SOLID TO FLUID, ORTHO
  curi = istfortho;
  str = [str, '    ',mr,'{',num2str(ncoef*nsamples),'}{*}{',stftag,'} & ',mr,'{',num2str(ncoef*nsamples),'}{*}{0} & ',mr,'{',num2str(nsamples),'}{*}{',stfptp,'} & ',mr,'{',num2str(nsamples),'}{*}{',numcmd,'{',pnum(siunitxavailable,saveT{curi}),'}} & ',numcmd,'{',pnum(siunitxavailable,effective_T{curi}),'} & ',pnumpct(numcmd,err_T{curi}),' \\']; str=andl(str);
  % str = [str, '    & & & & ',numcmd,'{9.9963e-7} & 0.012 \\']; str=andl(str);
  % str = [str, '    & & & & ',numcmd,'{9.9984e-7} & 0.009 \\']; str=andl(str);
  str = [str, '    & & ',mr,'{',num2str(nsamples),'}{*}{',stfprp,'} & ',mr,'{',num2str(nsamples),'}{*}{',numcmd,'{',pnum(siunitxavailable,saveR{curi}(1)),'}} & ',numcmd,'{',pnum(siunitxavailable,effective_R{curi}(1)),'} & ',pnumpct(numcmd,err_R{curi}(1)),' \\']; str=andl(str);
  % EVENTUAL OTHER MEASUREMENT POINT
  % EVENTUAL OTHER MEASUREMENT POINT
  str = [str, '    & & ',mr,'{',num2str(nsamples),'}{*}{',stfprs,'} & ',mr,'{',num2str(nsamples),'}{*}{',numcmd,'{',pnum(siunitxavailable,saveR{curi}(2)),'}} & ',numcmd,'{',pnum(siunitxavailable,effective_R{curi}(2)),'} & ',pnumpct(numcmd,err_R{curi}(2)),' \\']; str=andl(str);
  % EVENTUAL OTHER MEASUREMENT POINT
  % EVENTUAL OTHER MEASUREMENT POINT

  str = [str, '    \hline']; str=andl(str);

  % SOLID TO FLUID, SLANTED
  curi = istfslant;
  str = [str, '    ',mr,'{',num2str(ncoef*nsamples),'}{*}{',stftag,'} & ',mr,'{',num2str(ncoef*nsamples),'}{*}{',slantedtxt,'} & ',mr,'{',num2str(nsamples),'}{*}{',stfptp,'} & ',mr,'{',num2str(nsamples),'}{*}{',numcmd,'{',pnum(siunitxavailable,saveT{curi}),'}} & ',numcmd,'{',pnum(siunitxavailable,effective_T{curi}),'} & ',pnumpct(numcmd,err_T{curi}),' \\']; str=andl(str);
  % EVENTUAL OTHER MEASUREMENT POINT
  % EVENTUAL OTHER MEASUREMENT POINT
  str = [str, '    & & ',mr,'{',num2str(nsamples),'}{*}{',stfprp,'} & ',mr,'{',num2str(nsamples),'}{*}{',numcmd,'{',pnum(siunitxavailable,saveR{curi}(1)),'}} & ',numcmd,'{',pnum(siunitxavailable,effective_R{curi}(1)),'} & ',pnumpct(numcmd,err_R{curi}(1)),' \\']; str=andl(str);
  % EVENTUAL OTHER MEASUREMENT POINT
  % EVENTUAL OTHER MEASUREMENT POINT
  str = [str, '    & & ',mr,'{',num2str(nsamples),'}{*}{',stfprs,'} & ',mr,'{',num2str(nsamples),'}{*}{',numcmd,'{',pnum(siunitxavailable,saveR{curi}(2)),'}} & ',numcmd,'{',pnum(siunitxavailable,effective_R{curi}(2)),'} & ',pnumpct(numcmd,err_R{curi}(2)),' \\']; str=andl(str);
  % EVENTUAL OTHER MEASUREMENT POINT
  % EVENTUAL OTHER MEASUREMENT POINT

  str = [str, '    \hline']; str=andl(str);

  % FLUID TO SOLID, ORTHO
  curi = iftsortho;
  str = [str, '    ',mr,'{',num2str(ncoef*nsamples),'}{*}{',ftstag,'} & ',mr,'{',num2str(ncoef*nsamples),'}{*}{0} & ',mr,'{',num2str(nsamples),'}{*}{',ftsptp,'} & ',mr,'{',num2str(nsamples),'}{*}{',numcmd,'{',pnum(siunitxavailable,saveT{curi}(1)),'}} & ',numcmd,'{',pnum(siunitxavailable,effective_T{curi}(1)),'} & ',pnumpct(numcmd,err_T{curi}(1)),' \\']; str=andl(str);
  % EVENTUAL OTHER MEASUREMENT POINT
  % EVENTUAL OTHER MEASUREMENT POINT
  str = [str, '    & & ',mr,'{',num2str(nsamples),'}{*}{',ftspts,'} & ',mr,'{',num2str(nsamples),'}{*}{',numcmd,'{',pnum(siunitxavailable,saveT{curi}(2)),'}} & ',numcmd,'{',pnum(siunitxavailable,effective_T{curi}(2)),'} & ',pnumpct(numcmd,err_T{curi}(2)),' \\']; str=andl(str);
  % EVENTUAL OTHER MEASUREMENT POINT
  % EVENTUAL OTHER MEASUREMENT POINT
  str = [str, '    & & ',mr,'{',num2str(nsamples),'}{*}{',ftsprp,'} & ',mr,'{',num2str(nsamples),'}{*}{',numcmd,'{',pnum(siunitxavailable,saveR{curi}),'}} & ',numcmd,'{',pnum(siunitxavailable,effective_R{curi}),'} & ',pnumpct(numcmd,err_R{curi}),' \\']; str=andl(str);
  % EVENTUAL OTHER MEASUREMENT POINT
  % EVENTUAL OTHER MEASUREMENT POINT

  str = [str, '    \hline']; str=andl(str);

  % FLUID TO SOLID, SLANTED
  curi = iftsslant;
  str = [str, '    ',mr,'{',num2str(ncoef*nsamples),'}{*}{',ftstag,'} & ',mr,'{',num2str(ncoef*nsamples),'}{*}{',slantedtxt,'} & ',mr,'{',num2str(nsamples),'}{*}{',ftsptp,'} & ',mr,'{',num2str(nsamples),'}{*}{',numcmd,'{',pnum(siunitxavailable,saveT{curi}(1)),'}} & ',numcmd,'{',pnum(siunitxavailable,effective_T{curi}(1)),'} & ',pnumpct(numcmd,err_T{curi}(1)),' \\']; str=andl(str);
  % EVENTUAL OTHER MEASUREMENT POINT
  % EVENTUAL OTHER MEASUREMENT POINT
  str = [str, '    & & ',mr,'{',num2str(nsamples),'}{*}{',ftspts,'} & ',mr,'{',num2str(nsamples),'}{*}{',numcmd,'{',pnum(siunitxavailable,saveT{curi}(2)),'}} & ',numcmd,'{',pnum(siunitxavailable,effective_T{curi}(2)),'} & ',pnumpct(numcmd,err_T{curi}(2)),' \\']; str=andl(str);
  % EVENTUAL OTHER MEASUREMENT POINT
  % EVENTUAL OTHER MEASUREMENT POINT
  str = [str, '    & & ',mr,'{',num2str(nsamples),'}{*}{',ftsprp,'} & ',mr,'{',num2str(nsamples),'}{*}{',numcmd,'{',pnum(siunitxavailable,saveR{curi}),'}} & ',numcmd,'{',pnum(siunitxavailable,effective_R{curi}),'} & ',pnumpct(numcmd,err_R{curi}),' \\']; str=andl(str);
  % EVENTUAL OTHER MEASUREMENT POINT
  % EVENTUAL OTHER MEASUREMENT POINT

  str = [str, '    \hline']; str=andl(str);


  str = [str, '\end{tabular}']; str=andl(str);

  clipboard('copy', str);
  
  disp(['[',mfilename,'] LaTeX tabular copied to clipboard. CTRL+V somewhere to see.']);
end

function str = pnum(siunitxavailable,x)
  sigdig = 4;
  
  if(log10(x)<-2)
    if(siunitxavailable)
      str = sprintf(['%.',num2str(sigdig),'e'], x);
    else
      str = scientific_latex_notation(x, sigdig, 1);
    end
  else
    str = sprintf(['%.',num2str(sigdig),'f'], x);
  end
end
function str = pnumpct(numcmd,x)
  sigdig = 3;
  if(isnan(x))
    str = '--';
  else
    str = ['',numcmd,'{',sprintf(['%.',num2str(sigdig),'f'], x*100), '}'];
  end
end

function str=andl(str)
  str = [str, sprintf('\n')];
end