function index=subsampleDataArea(lon,lat,grdlen,minum)
  %% index=subsampleDataArea(lon,lat,grdlen,minum)
  %
  % Returns the indices of randomly picked data points based on
  % Number of points in the equal-area cell.
  %
  % INPUT:
  %
  % lon     longitudinal coordinates of the points [degrees]
  % lat     latitudinal coordinates of the points [degrees]
  % grdlen  grid sides on equator in degrees
  % minum   min number of points
  %
  % OUTPUT:
  %
  % index   logical index with "true" for selected points, "false" for thrown out
  %
  % Last modified by plattner-at-alumni.ethz.ch, 10/22/2017

  addpath('./Frederik')	

  rplanet=1% It's just to show the unit cell area

  %% Create equal area grid
  epslat=grdlen/10;

  c11=[min(lon),max(lat)+epslat];
  cmn=[max(lon),min(lat)-grdlen];

  [latgrid,dlongrid,refarea,nmr]=authalic(c11,cmn,grdlen,grdlen,rplanet);

  %% For each point, find out in which cell
  [celnr,rownr,colnr]=acor2ind(latgrid,dlongrid,nmr,c11,[lon(:) lat(:)]);

  %% Now count how many points in each cell  
  hgrd=-0.5:max(celnr);
  cellcount=histc(celnr,hgrd);	
  
  %% Now give each point the probability minum/ncell
  p = minum./max(1,cellcount(celnr+1));
  
  %% Create uniformly randomly selected values between 0 and 1
  pcomp=unifrnd(0,1,length(p),1);
  
  %% We keep the points where the probability p is 
  %% greater than the randomly picked value pcomp
  index= (p>pcomp);
  