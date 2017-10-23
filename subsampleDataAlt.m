function index=subsampleDataAlt(lon,lat,alt,alongtr,rplanet)
  %%  index=subsampleDataAlt(lon,lat,alt,alongtr,rplanet)
  %
  % Subsample based on altitude (the higher up, the lower we need the sampling rate). 
  %
  % INPUT:
  %
  % lon        vector containing longitudinal positions [in degrees]
  % lat        vector containing latitudinal positions [in degrees]
  % alt        vector containing altitude above planet surface [in km]
  % alongtr    along-track sampling interval [in degrees]
  % rplanet    planet radius [in km]
  %
  % OUTPUT:
  %
  % index      logical indices True (1) for points to keep
  %            (modeling data) and False (0) for the testing data
  %
  % Developed with Catherine Johnson at UBC
  % Last modified by plattner-at-alumni.ethz.ch, 10/23/2017
  
  %% Constants
  alongdist=(2*pi*rplanet)/(alongtr*360); % km 
  
  %% Calculate probability to be selected based on altitude
  % Idea: we can resolve roughly the same length scale on the
  % surface as our height above the surface, and we need    
  p=alongdist./(alt/2);

  %% Create a vector of uniformly randomly distributed numbers
  % to select which points live
  pcomp=unifrnd(0,1,length(p),1);

  index=( p>pcomp );
  