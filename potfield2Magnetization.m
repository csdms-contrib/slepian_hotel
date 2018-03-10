function mcoef=potfield2Magnetization(coef,rp,dbot)
%%% mcoef=potfield2Magnetization(coef,rp,dbot)
%
% Calculates the coefficients for the magnetization similar to
% the work by Gubbins et al. (2011) but here we assume a
% magnetized layer between depth dtop and dbot. See magnetization.pdf
%
% Theoretically, it is possible to allow for dtop and dbot. For large l,
% however, numerical inaccuracies create problems. We therefore, for noe,
% only allow dbot
%
% INPUT:
%
% coef       spherical-harmonic coefficients describing the
%            potential field on the planet's surface
%            the physical unit of coef is typically [nT * km]
% rp         radius of planet [in km]
% dbot       depth to the bottom of magnetized layer (thickness)
%            [in km]  
%
% OUTPUT:
%
% mcoef      coefficients for Flm to describe the magnetization's
%            lateral variation
%
% Last modified by plattner-at-alumni.ethz.ch, 03/08/2018
%
%%% NOT YET SURE ABOUT THE UNIT TRANSITIONS HERE!!!
%%% See magnetization.pdf
  
mu0=4*pi*1e-7;
  
Lmax=sqrt(length(coef))-1;
[~,~,~,~,~,~,~,l]=addmon(Lmax);

%%% Turn coef into the right units.
%%% It's typically given in nT * km. Turn it into T*m
coef=coef*1e-9*1e3;
%%% Because coef was obtained through integration against Y*,
%%% we need to change the km into m with the units of Ylm:
coef=coef*1e-6;

%%% Turn rp into meters
rp=rp*1e3;
%%% Turn dbot into meters
dbot=dbot*1e3;



%%%fact=rp.^(l+1).*sqrt(l.*(2*l+1))./(mu0*( (rp-dtop).^l - (rp-dbot).^l ));
%%% Because of numeical problems we need to rearrange this
%%%fact=(rp.^l)./( (rp-dtop).^l - (rp-dbot).^l ) .*sqrt(l.*(2*l+1)) *rp/mu0;
%%% Even after rearrangement: With high l, this breaks. Therefore:
%%% Don't allow dtop anymore. Then you can do
fact= 1./( 1 - ((rp-dbot)/rp).^l  ) .*sqrt(l.*(2*l+1)) *rp/mu0;


mcoef=coef.*fact;

%%% The output coefficients are for the Flm, which have no l=0.
%%% Here, fact(1)=0 anyhow. So we just remove the l=0 coef.
mcoef=mcoef(2:end);


%%% The unit of F* in the expansion is 1/m^2,
%%% so we need to multiply with 1e-6 to go from 1/km^2
mcoef=mcoef*1e-6;


%%% I don't yet understand why, but my magnetization results were
%%% always sqrt(4*pi) off from other people's results. It is not because 
%%% of the Flm2xyz plotting routine normalization. Same thing with flm,
%%% which is unit-sphere normalized.
mcoef=mcoef/sqrt(4*pi);
