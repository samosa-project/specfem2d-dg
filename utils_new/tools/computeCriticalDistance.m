function [criticaldistance, angleRad, angleDeg, SOUNDSPEED, vpvs] = computeCriticalDistance(atmfile,parfile,sourcefile)
  [~,~,~,SOUNDSPEED]=extract_atmos_model(atmfile,3,0,-1); SOUNDSPEED=SOUNDSPEED(1); % get C at ground
  vpvs = readExampleFiles_extractParfileModels(parfile); % get parfile models;
  vpvs = vpvs(2,:); disp(['[',mfilename,'] From parfile, assume first ground model is numbered 2.']);
  vpvs = vpvs(4:5); % extract vp and vs from this model
  [angleRad, angleDeg] = criticalAirToGroundTransmissionAngle(SOUNDSPEED, vpvs(1), vpvs(2)); % get critical angle.
  zsource = readExampleFiles_extractParam(sourcefile, 'zs', 'float'); % get source altitude
  % PYTHAGORUS, USE TRIGONOMETRY! (find the distance |x| corresponding to a line launched at angle angleRad from the vertical from the source)
  criticaldistance = zsource * tan(angleRad);
  % IT'S VERY EFFECTIVE
  disp(['[',mfilename,'] Critical distance (over which no atmospheric wave can be transmitted to the ground), with']);
  disp([blanks(length(mfilename)+2),'     c  = ',sprintf('%6.1f',SOUNDSPEED),' [m/s],']);
  disp([blanks(length(mfilename)+2),'     vp = ',sprintf('%6.1f',vpvs(1)),' [m/s],']);
  disp([blanks(length(mfilename)+2),'     vs = ',sprintf('%6.1f',vpvs(2)),' [m/s],']);
  disp([blanks(length(mfilename)+2),' is d^c = ',sprintf('%6.1f',criticaldistance),' [m] (critical incident angle being ',num2str(angleDeg),'° from vertical).']);
end