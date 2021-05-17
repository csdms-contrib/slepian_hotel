function index=subsampleDataTrackgap(lon,lat,alt,alongtr,betweentr,rplanet)
  %%  index=subsampleDataTrackgap(lon,lat,alt,alongtr,betweentr,rplanet)
  %
  % Subsample based on along-track vs between-track
  % distance. We assume that the tracks are roughly equally spaced
  % and are running north-south. 
  %
  % INPUT:
  %
  % lon        vector containing longitudinal positions [in degrees]
  % lat        vector containing latitudinal positions [in degrees]
  % alt        vector containing altitude above planet surface [in km]
  % alongtr    along-track sampling interval [in degrees]
  % betweentr  between-track sampling interval [in degrees]
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
  betweendisteq=(2*pi*rplanet)/(betweentr*360);  
  

  %% Calculate probability to be selected based on between-track filtering
  p = alongdist./(betweendisteq*cos(lat*pi/180));

  
  %% Create a vector of uniformly randomly distributed numbers
  % to select which points live
  pcomp=rand(length(p),1);

  index=( p>pcomp );
  




