function [data,lon,lat,rad,index]=cut2cap(data,lon,lat,rad,cTH,clon,clat)
% [datacut,loncut,latcut,radcut,index]=cut2cap(data,lon,lat,rad,cTH,clon,clat)
% 
% Removes data points outside of a spherical cap or opening angle cTH
% centered at longitude clon and latitude clat.
%
% INPUT:
% 
% data      data{1} is radial component, data{2} is colatitudinal component
%           data{3} is longitudinal component
% lon, lat  longitudinal, latitudinal positions of the data points
%           (in degrees)
% rad       radial position of the data points
% cTH       spherical cap semi-opening angle
% clon      spherical cap center point longitude
% clat      spherical cap center point latitude
%
% OUTPUT
%
% The truncated data, lon, lat, and rad positions
% index     indices of points within the cap for the original vectors
%           (logical index, not list of indices) 
%
% Last modified by plattner-at-alumni.ethz.ch, 5/17/2017

% We need colatitude of center:
ccola=90-clat;
clear clat;

% To find the indices of the points we want: simply rotate all the
% points to the North Pole and then keep the ones with colatitudes
% smaller than the opening angle.

[thetap,~]=rottp((90-lat)*pi/180,lon*pi/180,clon*pi/180,ccola*pi/180,0);
% Now find out which ones are close enough to the north pole
% Use logical indexing for this
index=(thetap<=cTH*pi/180);

for cmp=1:3
    data{cmp}=data{cmp}(index);
end
lon=lon(index);
lat=lat(index);
rad=rad(index);
