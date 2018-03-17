function coef=importMag(fname,rplanet,Lmin,Lmax)
% coef=importMag(fname,rplanet,Lmin,Lmax)
%
% Imports the crustal magnetic field coefficients from the
% standard format and converts them to our potential-field
% spherical-harmonic normalization.
%
% Delete all header lines beforehand
%  
% INPUT:
%
% fname       name of the data file
% rplanet     planet radius
% Lmin        minimum sperical-harmonic degree
%             given in the data file
% Lmax        maximum sperical-harmonic degree
%             given in the data file
%
% Last modified by plattner-at-alumni.ethz.ch, 03/17/2018

  [~,~,~,~,~,~,~,bigl]=addmon(Lmax);

  magdat=importdata(fname);

  lmcosi=magdat(:,1:4);

  if Lmin>0
    lmcosi=[zeros(Lmin^2,4);lmcosi];
  end
    
  coef=lmcosi2coef(lmcosi);

  coef=coef./sqrt(2*bigl+1)*rplanet*(-1);

  % The coefficients are now in 4pi normalized form.
  % We want them in unit-normailzed form
  coef=coef*sqrt(4*pi);

  
