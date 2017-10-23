function index=subsampleDataTrackgapAlt(lon,lat,alt,alongtr,betweentr,rplanet)
  %%  index=subsampleDataTrackgapAlt(lon,lat,alt,alongtr,betweentr,rplanet)
  %
  % Subsample based on altitude and along-track vs between-track
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
  % Last modified by plattner-at-alumni.ethz.ch, 10/22/2017
  
  %% Constants
  alongdist=(2*pi*r)/(alongtr*360); % km
  betweendisteq=(2*pi*r)/(betweentr*360);  
  
  %% Calculate probability to be selected based on altitude
  % Idea: we can resolve roughly the same length scale on the
  % surface as our height above the surface, and we need    
  p1=alongdist./(alt/2);

  %% Calculate probability to be selected based on between-track filtering

  p2 = alongdist./(betweendisteq*cos(lat*pi/180));

  %% We take the lower of the two probabilities
  p=min(p1,p2);
  
  %% Create a vector of uniformly randomly distributed numbers
  % to select which points live
  pcomp=unifrnd(0,1,length(p1),1);

  index=( p>pcomp );
  
  data=data(index,:);




