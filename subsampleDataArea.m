function index=subsampleDataArea(lon,lat,dlon,dlat,minum)
  %% index=subsampleDataArea(lon,lat,dlon,dlat,minum)
  %
  % Returns the indices of randomly picked data points based on
  % Number of points in the equal-area cell.
  %
  % INPUT:
  %
  % lon     longitudinal coordinates of the points [degrees]
  % lat     latitudinal coordinates of the points [degrees]
  % dlon    longitudinal grid side length on equator [degrees]
  % dlat    latitudinal grid side length
  % minum   min number of points
  %
  % OUTPUT:
  %
  % index   logical index with "true" for selected points, "false" for thrown out
  %
  % Last modified by plattner-at-alumni.ethz.ch, 10/23/2017

  %% Create equal area grid
  epslat=dlat/10;

  c11=[min(lon),max(lat)+epslat];
  cmn=[max(lon),min(lat)-dlat];

  [latgrid,dlongrid,refarea,nmr]=authalic(c11,cmn,dlat,dlon);

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
  
