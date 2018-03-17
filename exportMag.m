function exportMag(coef,fname,rplanet)
% exportMag(coef,fname,rplanet)  
%
% Exports our potential field coefficients to the
% standard format. Don't forget to divide by sqrt(4*pi)
% if you obtained the coefficients from an inversion
% with LocalIntField or LocalIntExtField  
%
% INPUT:
%
% coef      potential field coefficients (normalized for
%           spherical-harmonics that have norm sqrt(4*pi)
%           over the unit sphere)  
% fname     file name for coefficients
% rplanet   planet radius
%
% Last modified by plattner-at-alumni.ethz.ch, 03/17/2018

  Lmax=sqrt(length(coef))-1;   
  [~,~,~,~,~,~,~,bigl]=addmon(Lmax);

  % Our coefficients are unit normalized. We need them
  % in 4 pi normalized
  coef=coef/sqrt(4*pi);
  
  coef=coef.*sqrt(2*bigl+1)*(-1)/rplanet;

  lmcosi=coef2lmcosi(coef);

  dlmwrite(fname,lmcosi,'\t')
