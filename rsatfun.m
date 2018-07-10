function r=rsatfun(x)
% r=rsatfun(x)
%
% Example/Template for a function of latitudinally varying
% radial position forthe altitude-cognizant Slepian functions 
% calculated in kernelcppotupLatvar
%
% Original idea by Erwan Mazarico
%  
% INPUT:
%
% x       cos(colatitude) or sin(latitude)
%
% OUTPUT:
%
% y       radial position of the satellite
%  
% Last modified by plattner-at-alumni.ethz.ch, 7/9/2018

% Simple linear example:
r = 3390 + 100 + (1-x)*300;   
%r = 3390 + 100 - (1-x)*100; 
  
