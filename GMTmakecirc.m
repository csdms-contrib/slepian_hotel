function GMTmakecirc(TH,clon,ccola,fname)
% GMTmakecirc(TH,clon,ccola,fname)
%
% Writes out a circle on the sphere as a list of longitude and latitude
% values
%
% INPUT:
%
% TH         Cap opening angle in degrees (see Plattner & Simons 2014 ACHA)
% clon       Cap center longitude (degrees)
% ccola       Cap center colatitude (degrees)
% fname      Name of file to export to
%
% Last modified by plattner-at-alumni.ethz.ch, 05/16/2014

[lon,cola]=circonsphere(pi/180*clon,pi/180*ccola,2*TH,360/2/pi);


lon=180/pi*lon;
lat=90-180/pi*cola;
% Must have column vectors. Otherwise DESASTER
lon=lon(:);
lat=lat(:);

tab=table(lon,lat);
writetable(tab,fname,'Delimiter','\t');