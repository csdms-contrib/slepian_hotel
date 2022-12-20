function [G,V] = glmalphaPotSurf(TH,Lmax,rsat,rplanet)
% [G,V] = glmalphaPotSurf(TH,Lmax,rsat,rplanet)
%
% Slepian functions that are concentrated on the planet's surface but use data
% at satellite altitude. From potential field to potential field.
%  
% INPUT:
%
% TH       Angular extent of the spherical cap, in degrees OR
%          'england', 'eurasia',  'namerica', 'australia', 'greenland'
%          'africa', 'samerica', 'amazon', 'orinoco', 'antarctica', 'alloceans':
%          OR: [lon lat] an ordered list defining a closed curve [degrees]
%          OR: Angles of two spherical caps and we want the ring between 
%              them [TH1 TH2] 
% Lmax     Bandwidth (maximum angular degree), or passband (two degrees)
% rsat     radius of satellite altitudedata location
% rplanet  radius of planet
% srt      sorted output?
%
% OUTPUT:
%
% G        The unitary matrix of localization coefficients; note how
%          LOCALIZATION delivers these as LMCOSI arrays into PLM2XYZ
% V        The eigenvalues in this ordering (not automatically sorted)

% Last modified by plattner-at-alumni.ethz.ch, 07/29/2022  
 
  
TH = 'namerica';
Lmax=20;
rplanet = fralmanac('a_EGM96');
rsat = rplanet + 1000000;

% Calculate Slepian functions for A^-1 K A^-1
% For A = (rplanet/rsat)^(l+1) I can easily do this by just switching
% rplanet and rsat in the function call, which asks for rsat, then rplanet
[F,V] = glmalphapotup(TH,Lmax,rplanet,rsat);

% Get g from g = A^(-1)f
G = potup(F,rplanet,rsat,Lmax,0);
%G2 = potup(G,rplanet,rsat,Lmax,0);

figure(1)
clf
plotslep(F,15)

figure(2)
clf
plotslep(G,15)

keyboard