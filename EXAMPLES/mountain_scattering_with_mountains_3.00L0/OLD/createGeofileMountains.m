function [] = createGeofileMountains(outputFile, meshAlgorithm, xbound, zbound, nLayersGround, groundLayersHeightProportion, nPeaks, mwrp, mwlp, distRight, peaksAltitude, valleysAltitude, dxAir, dxGrdInts)
  if(not(exist('outputFile', 'var')))
    error('kek');
    outputFile = 'C:\Users\yanis\Desktop\ISAF\codesSpecfem\extMesh.geo';
  end
  if(not(exist('nLayersGround', 'var')))
    nLayersGround = 3;
  end
  if(not(exist('nPeaks', 'var')))
    nPeaks = 2; % Number of peaks of the mountain
  end
  if(not(exist('mwrp', 'var')))
    mwrp = [7,10];
  end
  if(not(exist('mwlp', 'var')))
    mwlp = [10,15];
  end
  if(not(exist('peaksAltitude', 'var')))
    peaksAltitude = [10,20];
  end
  if(not(exist('xbound','var')))
    xbound = [0,100]; 
  end
  if(not(exist('zbound','var')))
    zbound = [-50, 0, 70]; 
  end
  if(not(exist('groundLayersHeightProportion','var')))
  %     groundLayersHeightProportion = 1./nLayersGround*ones(1,nLayersGround);
    groundLayersHeightProportion = [2/5,2/5,1/5];
  end
  if(not(exist('valleysAltitude','var')))
    valleysAltitude = [5,0]; 
  end
  
  fid = fopen(outputFile, 'w');

  % Maillage
  xmin = xbound(1);
  xmax = xbound(2);
  zmin = zbound(1);
  zmax = zbound(3);
  zint = zbound(2);

  % Mountains
  % Height and width of each peaks
  mountain_1_x = xmax-distRight-0.5*mwrp(1);
  
  % Display parameters.
  showAxes = 0;
  displayLine = 1;
  displayLineNumber = 1;
  displayMeshLineNumber = 0;
  displayMeshPoints = 1;
  displayPointsNumber = 1;
  displayMeshPointsNumber = 0;
  displaySurfaceEdge = 1;
  displaySurfaceFaces = 1;
  displaySurfaceNumber = 1;
  displayMeshSurfaceNumber = 0;
  
  
  meshRandomFactor = 1e-4;
  
  nPointsTotal = (nLayersGround+1)*2 + nPeaks*3 - (nPeaks - 1) + 2;
  
  % Write display results
  fprintf(fid, '// Display. //---------------//\n');
  fprintf(fid, 'General.Axes = %d; // Show axes.\n', showAxes);
  fprintf(fid, 'Mesh.Lines = %d; // Display lines.\n', displayLine);
  fprintf(fid, 'Geometry.LineNumbers = %d; // Display line numbers.\n', displayLineNumber);
  fprintf(fid, 'Mesh.LineNumbers = %d; // Display line numbers.\n', displayMeshLineNumber);
  fprintf(fid, 'Mesh.Points = %d; // Display points.\n', displayMeshPoints);
  fprintf(fid, 'Geometry.PointNumbers = %d; // Display point numbers.\n', displayPointsNumber);
  fprintf(fid, 'Mesh.PointNumbers = %d; // Display point numbers\n', displayMeshPointsNumber);
  fprintf(fid, 'Mesh.SurfaceEdges = %d; // Display edges.\n', displaySurfaceEdge);
  fprintf(fid, 'Mesh.SurfaceFaces = %d; // Display surfaces.\n', displaySurfaceFaces);
  fprintf(fid, 'Geometry.SurfaceNumbers = %d; // Display surface numbers.\n', displaySurfaceNumber);
  fprintf(fid, 'Mesh.SurfaceNumbers = %d; // Display surface numbers.\n', displayMeshSurfaceNumber);
  fprintf(fid, '// Elements. //--------------//\n');
  fprintf(fid, 'Mesh.SubdivisionAlgorithm = 1; // Meshing algorithm: all quads.\n');
  fprintf(fid, 'Mesh.RecombineAll = 1; // Recombine all triangular elements.\n');
  fprintf(fid, '// Meshing algorithm. //-----//\n');
  fprintf(fid, 'Mesh.Algorithm = %d; // 1: MeshAdapt. 5: Delaunay for quads. 8: Delaunay (experimental). 9: structured (packing of parallelograms, experimental).\n', meshAlgorithm);
  fprintf(fid, 'Mesh.ElementOrder = 1; // Element order.\n');
  fprintf(fid, 'Mesh.RandomFactor = %f; // Perturbate every points positions in order to avoid 3 aligned points. Default is at 1e-9, maximum is 1e-3.\n', meshRandomFactor);
  fprintf(fid, '// Geometry. //--------------//\n');
  fprintf(fid, '// Physical parameters.\n');
  fprintf(fid, 'dxAirBot = %f;\n', dxAir(1));
  fprintf(fid, 'dxAirTop = %f;\n', dxAir(2));
  for i = 1:nLayersGround
    fprintf(fid, 'dxGrdInts%d = %f;\n', i, dxGrdInts(i));
  end
  fprintf(fid, '// Define geometry keypoints.\n');
  airBotPts = []; airTopPts = []; % points list for affecting dx
  grdPts = cell(nLayersGround,1);
  % Ground layers points
  curPt = 1;
  zsum = 0;
  for l = 0:nLayersGround
    if l>0
      zsum = zsum + (zint-zmin)*groundLayersHeightProportion(l);
    end
    fprintf(fid, 'Point(%d) = {%f, %f, 0}; // corner\n', curPt, xmin, zmin + zsum);
    if(zmin+zsum==zint)
      airBotPts=[airBotPts, curPt]; % at left interface, add current point to air points list
    else
      grdPts{l+1}=[grdPts{l+1}, curPt];
    end
    curPt = curPt+1;
    fprintf(fid, 'Point(%d) = {%f, %f, 0}; // corner\n', curPt, xmax, zmin + zsum);
    if(zmin+zsum==zint)
      airBotPts=[airBotPts, curPt]; % at right interface, add current point to air points list
    else
      grdPts{l+1}=[grdPts{l+1}, curPt];
    end
    curPt = curPt+1;
  end
  % Mountain points.
  fprintf(fid, 'Point(%d) = {%f, %f, 0}; // mountain point \n', curPt, mountain_1_x + 0.5*mwrp(1) , zint); airBotPts=[airBotPts, curPt]; curPt = curPt + 1;
  for i = 1:nPeaks
    fprintf(fid, 'Point(%d) = {%f, %f, 0}; // mountain point \n', curPt, mountain_1_x - mwrp(i)*(i-1) , zint + peaksAltitude(i)); airBotPts=[airBotPts, curPt]; curPt = curPt + 1;
    fprintf(fid, 'Point(%d) = {%f, %f, 0}; // mountain point \n', curPt, mountain_1_x - 0.5*mwlp(i) - mwrp(i)*(i-1) , zint+valleysAltitude(i)); airBotPts=[airBotPts, curPt]; curPt = curPt + 1;
  end
  % Finish box.
  fprintf(fid, 'Point(%d) = {%f, %f, 0};\n', curPt, xmin, zmax); airTopPts=[airTopPts, curPt]; curPt = curPt + 1;
  fprintf(fid, 'Point(%d) = {%f, %f, 0};\n', curPt, xmax, zmax); airTopPts=[airTopPts, curPt]; curPt = curPt + 1;
  
  % Lines.
  fprintf(fid, '// Define lines joining points. Warning: prefer counter-clockwise direction to prevent inverted elements. Warning: do not duplicate lines if two or more materials are needed.\n');
  curLine = 1;
  % Each rectangle.
  fprintf(fid, 'Line(%d) = {%d, %d};\n', curLine, 1 , 2);
  curLine = curLine + 1;
  for i = 1: nLayersGround-1
    fprintf(fid, 'Line(%d) = {%d, %d};\n', curLine, 2*i , 2*(i+1)); curLine = curLine + 1;
    fprintf(fid, 'Line(%d) = {%d, %d};\n', curLine, 2*(i+1) , 2*(i+1) - 1); curLine = curLine + 1;
    fprintf(fid, 'Line(%d) = {%d, %d};\n', curLine, 2*(i+1)-1 , 2*(i-1)+1); curLine = curLine + 1;
  end
  fprintf(fid, 'Line(%d) = {%d, %d};\n', curLine, 2*nLayersGround , 2*(nLayersGround+1)); curLine = curLine + 1;
  fprintf(fid, 'Line(%d) = {%d, %d};\n', curLine, 2*(nLayersGround+1) , 2*(nLayersGround+1)+1); curLine = curLine + 1;
  % Mountains.
  ref = 2*(nLayersGround+1)+1;
  for i = 1:nPeaks
    fprintf(fid, 'Line(%d) = {%d, %d};\n', curLine, ref + 2*(i-1) , ref + 2*(i-1)+1); curLine = curLine + 1;
    fprintf(fid, 'Line(%d) = {%d, %d};\n', curLine, ref + 2*(i-1)+1 , ref + 2*i); curLine = curLine + 1;
  end
  % Finish the big rectanle.
  fprintf(fid,'Line(%d) = {%d, %d};\n', curLine, ref + 2*nPeaks , 2*nLayersGround+1);
  curLine = curLine + 1;
  % Air.
  fprintf(fid,'Line(%d) = {%d, %d};\n', curLine, 2*(nLayersGround+1) , nPointsTotal); curLine = curLine + 1;
  fprintf(fid,'Line(%d) = {%d, %d};\n', curLine, nPointsTotal , nPointsTotal-1); curLine = curLine + 1;
  fprintf(fid,'Line(%d) = {%d, %d};\n', curLine, nPointsTotal-1 , 2*(nLayersGround+1) -1); curLine = curLine + 1;
  fprintf(fid,'Line(%d) = {%d, %d};\n', curLine, 2*(nLayersGround+1) -1 , 2*(nLayersGround) -1); curLine = curLine + 1;
  
  % Surfaces.
  curLineLoop = 1;
  % Air line loop.
  lineBeforeRight = 3*(nLayersGround-1)+2;
  lineFirstBeforeMountain = lineBeforeRight + 1;
  lineFirstAfterMountain = lineFirstBeforeMountain + 2*nPeaks + 1;
  lineRightAir = lineFirstAfterMountain + 1;
  lineTopAir = lineRightAir + 1;
  lineLeftAir = lineTopAir + 1;
  fprintf(fid, 'Line Loop(%d) = {%d,%d,%d,', curLineLoop, lineRightAir, lineTopAir, lineLeftAir);
  curLineLoop = curLineLoop + 1;
  refPeaks = lineFirstAfterMountain;
  for i = 0: 2*nPeaks
    fprintf(fid,'-%d,',refPeaks - i);
  end
  fprintf(fid,'-%d};\n',refPeaks - 2*nPeaks-1);
  % Ground.
  for i = 1:nLayersGround-1
%     fprintf(fid, 'Line Loop(%d) = {', curLineLoop);
    fprintf(fid, 'Line Loop(%d) = {', 2*nLayersGround-curLineLoop);
    curLineLoop = curLineLoop + 1;
    if i ==  1
      First = 1;
      Second = 2;
      Third = 3;
      Fourth = 4;
      fprintf(fid, '%d, %d, %d, %d};\n', First, Second, Third, Fourth);
    else
      First = 3*(i-1);
      Second = 2+3*(i-1);
      Third = Second + 1;
      Fourth = Third + 1;
      fprintf(fid, '-%d, %d, %d, %d};\n', First, Second, Third, Fourth);
    end
  end
%   fprintf(fid, 'Line Loop(%d', curLineLoop);
  fprintf(fid, 'Line Loop(%d', 2*nLayersGround-curLineLoop);
  fprintf(fid, ') = {-%d, %d, %d, ', (nLayersGround-1)*3, (nLayersGround-1)*3+2, (nLayersGround-1)*3+3);
  for i = 1:nPeaks
    ref = (nLayersGround-1)*3+3;
    fprintf(fid, '%d, %d, ', ref+2*i-1, ref+2*i);
  end
  fprintf(fid, '%d, %d}; \n', ref + nPeaks*2+1, ref + nPeaks*2+1 + 4);
  for i = 1:curLineLoop
    fprintf(fid, 'Plane Surface(%d) = {%d};\n', i , i);
  end
  
  % Elements' size.
  fprintf(fid, 'Mesh.CharacteristicLengthFactor = 1; // Scaling of all elements'' sizes. \n');
  dxTxt = sprintf('%d, ', airBotPts); fprintf(fid, ['Characteristic Length {',dxTxt(1:end-2),'} = 2*dxAirBot; \n']); % air bottom
  dxTxt = sprintf('%d, ', airTopPts); fprintf(fid, ['Characteristic Length {',dxTxt(1:end-2),'} = 2*dxAirTop; \n']); % air top
  for l=0:(nLayersGround-1)
    dxTxt = sprintf('%d, ', grdPts{l+1}); fprintf(fid, ['Characteristic Length {',dxTxt(1:end-2),'} = 2*dxGrdInts',sprintf('%d',l+1),'; \n']);
  end
  
  % Physical surfaces and lines.
  fprintf(fid, 'Physical Line("Top") = {%d};\n', lineTopAir);
  fprintf(fid, 'Physical Line("Left") = {%d, %d', lineLeftAir, lineLeftAir +1);
  if nLayersGround ==  1
    fprintf(fid, '};\n'); 
  else
    for i = 1:nLayersGround-1
       fprintf(fid, ',%d', 3*(nLayersGround-i)+1 );
    end
    fprintf(fid,'};\n');
  end
  fprintf(fid, 'Physical Line("Bottom") = {1};\n');
  fprintf(fid, 'Physical Line("Right") = {%d', 2);
  if nLayersGround > 1
    for i = 2: nLayersGround
      fprintf(fid, ',%d',  2+3*(i-1));
    end
  end
  fprintf(fid, ', %d}; \n', lineRightAir);
  for i = 1:curLineLoop
    fprintf(fid, 'Physical Surface("M%d") = {%d};\n',i,i); 
  end
  fprintf(fid, 'Recombine Surface {1');
  if curLineLoop>1
    for i = 1:curLineLoop
      fprintf(fid, ',%d',i);
    end
  end
  fprintf(fid, '};\n');
  fclose(fid);
end