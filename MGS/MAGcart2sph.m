function [lon,cola,r,data,Brms,Bstat,Bdyn,index]=MAGcart2sph(...
    posX,posY,posZ,magX,magY,magZ,RMSBvalsX,RMSBvalsY,RMSBvalsZ,...
    BstatX,BstatY,BstatZ,BdynX,BdynY,BdynZ,sphericalcut)
% function [lon,cola,r,data,Brms,Bstat,Bdyn]=MAGcart2sph(...
%    posX,posY,posZ,magX,magY,magZ,RMSBvalsX,RMSBvalsY,RMSBvalsZ,...
%    BstatX,BstatY,BstatZ,BdynX,BdynY,BdynZ,sphericalcut)
%
% Transforms the data directly loaded from the .STS files into longitude,
% latitude, radial position and the magnetic fields into radial direction,
% colatitudinal direction, longitudinal direction, etc.
%
% INPUT:
%
% posX, posY, posZ  location in coordinates of the .STS file
% magX, magY, magZ  magnetic field components in coordinate system defined
%                   in the  .STS file (e.g. Planetocentric or SunState)
%
% If available, provide the following:
% RMSBvalsX,Y,Z     Some error measure for B field (I think from binning)
% BstatX,Y,Z        Some measure for B static noise (I think modeled)
% BdynX,Y,Z         Some measure for B dynamic noise (I think modeled)
%
% If you would like to only retain data within a spherical cap:
% sphericalcut      [clon clat cTH]: center point (lat, lon) and HALF 
%                   opening angle TH (see Plattner and Simons, 2014, ACHA)
%                   for spherical cap whithin which to retain the
%                   data (everything else gets thrown out). All values in
%                   degrees 0<=clon<=360, -90<=clat<=90, 0<TH<180
%
% OUPUT:
%
% lon       Longitudinal position (radians) on the sastellite in Mars 
%           planetocentric coordinates (between 0 and 2pi)
% cola      Colatitudinal position (radians) on the sastellite in Mars 
%           planetocentric coordinates (between 0 and pi)
% r         Radial position [km] on the sastellite in Mars planetocentric
%           coordinates.
% data      radial component data{1}, latitudinal component data{2}, and
%           longitudinal component data{3} of the magnetic measurement in
%           planeticentric coordinate system. Values in nanotesla nT.
% Brms      radial Brms{1}, latitudinal Brms{2} and longitudinal Brms{3}
%           values. See loadpds for detail.
% Bstat     radial Bstat{1}, latitudinal Bstat{2} and longitudinal Bstat{3}
%           values. See loadpds for detail.
% Bdyn      radial Bdyn{1}, latitudinal Bdyn{2} and longitudinal Bdyn{3}
%           values. See loadpds for detail.
%
% Example:
%
% [lon,cola,r,data,Brms,Bstat,Bdyn,index]=MAGcart2sph(...
%    posX,posY,posZ,magX,magY,magZ,RMSBvalsX,RMSBvalsY,RMSBvalsZ,...
%    BstatX,BstatY,BstatZ,BdynX,BdynY,BdynZ,[0 -90 14]);
% sap=sap(index);sam=sam(index);sao=sao(index);date=date(index);
% save('transformedData','lon','cola','r','data',...
%       'Brms','Bstat','Bdyn','sap','sam','sao','date');
%
% Last modified by plattner-at-alumni.ethz.ch, 07/04/2016

defval('sphericalcut',[]);

defval('RMSBvalsX',[]);
defval('BstatX',[]);
defval('BdynX',[])

% Step 1: Transform location from cartesian to spherical.
[lon,lat,r] = cart2sph(posX,posY,posZ);
% The output of this function is in latitudes. We need colatitudes to
% proceed:
cola=pi/2-lat;
% The longitudes are between -pi and pi. Let's put them between 0 and 2pi.
% That's easy, just use modulo
lon=mod(lon,2*pi);

% We don't need posX,posY,posZ anymore. Delete them (they take up memory)
clear posX; clear posY; clear posZ;

% Also delete the latitudes to avoid confusion:
clear lat;

% Prepare index in case you don't want to cut: All of them.
index=1:length(lon);
% In case we only want data from within a spherical cap: cut it now:
if ~isempty(sphericalcut)
    clon=sphericalcut(1);
    clat=sphericalcut(2);
    cTH =sphericalcut(3);
    
    % Of course we need the center point in colatitude:
    ccola=90-clat;
    clear clat;       
    
    % To find the indices of the points we want: simply rotate all the
    % points to the North Pole and then keep the ones with colatitudes
    % smaller than the opening angle.
    
    [thetap,phip]=rottp(cola,lon,clon*pi/180,ccola*pi/180,0);
    % Now find out which ones are close enough to the north pole
    index=find(thetap<=cTH*pi/180);
    

    % Only retain the points that are within the circle:
    magX=magX(index);
    magY=magY(index);
    magZ=magZ(index);

    lon=lon(index);
    cola=cola(index);
    r=r(index);

	if ~isempty(RMSBvalsX)
    	RMSBvalsX=RMSBvalsX(index);
    	RMSBvalsY=RMSBvalsY(index);
    	RMSBvalsZ=RMSBvalsZ(index);
	end
	
	if ~isempty(BstatX)	
    	BstatX=BstatX(index);
    	BstatY=BstatY(index);
    	BstatZ=BstatZ(index);
	end

	if ~isempty(BdynX)
    	BdynX=BdynX(index);
    	BdynY=BdynY(index);
    	BdynZ=BdynZ(index);
	end
end

% Now transform cartesian magnetic data into spherical format
[data{3},data{2},data{1}]=dcart2dsph(lon,cola,magX,magY,magZ);

% Now transform simulated B values and binning RMS into spherical
% coordinates
if ~isempty(BdynX)
	[Bdyn{3},Bdyn{2},Bdyn{1}]=dcart2dsph(lon,cola,BdynX,BdynY,BdynZ);
else
	Bdyn=[];
end
if ~isempty(BstatX)	
	[Bstat{3},Bstat{2},Bstat{1}]=dcart2dsph(lon,cola,BstatX,BstatY,BstatZ);
else
	Bstat=[];
end
if ~isempty(RMSBvalsX)
	[Brms{3},Brms{2},Brms{1}]=dcart2dsph(lon,cola,RMSBvalsX,RMSBvalsY,RMSBvalsZ);
else
	Brms=[];
end
