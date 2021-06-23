function varargout=glmalphapotuptoJp(TH,L,rnew,rold,phi,theta,omega,J)
% [Grot,V]=glmalphapotuptoJp(TH,L,rnew,rold,phi,theta,omega,J)
%
% glmalphapotup for upward transformed scalar Slepian functions.
% Loads glmalphapotup and rotates the first J functions only.
%  
% This function is designed for the potential field at satellite altitude  
%
% INPUT:
%
% TH        Angular extent of the spherical cap (degrees)
%           OR: Angles of two spherical caps and we want the ring between 
%               them [TH1 TH2] 
% L         Bandwidth (maximum angular degree) or passband (two degrees)
% rnew      radius for radial component (at satellite altitude)
% rold      radius for scalar potential (on surface)
% phi       Longitude of the center (degrees) [default: 0]
% theta     Colatitude of the center (degrees) [default: 0]
% omega     Anticlockwise azimuthal rotation (degrees) [default: 0]
% J         The number of eigenfunctions that are being asked (and saved)
%
% OUTPUT:
%
% G        The unitary matrix of localization coefficients
% V        The eigenvalues 
%
% Last modified by plattner-at-alumni.ethz.ch, 27/06/2018

toout=1; % Should the result be in addmout?

defval('phi',0);
defval('theta',0);
defval('omega',0);

if length(L)==1
    if length(TH)==1
        fname=fullfile(getenv('IFILES'),'GLMALPHAPOTUPTOJP',...
             sprintf('glmalphapotuptoJp-%g-%i-%g-%g-%i-%i-%i-%i.mat',...
                 TH,L,rnew,rold,phi,theta,omega,J));
    elseif length(TH)==2
        fname=fullfile(getenv('IFILES'),'GLMALPHAPOTUPTOJP',...
             sprintf('glmalphapotuptoJp-%g_%g-%i-%g-%g-%i-%i-%i-%i.mat',...
                 max(TH),min(TH),L,rnew,rold,phi,theta,omega,J));
    else
        error('Bad choice for TH')
    end
elseif length(L)==2
    if length(TH)==1
        fname=fullfile(getenv('IFILES'),'GLMALPHAPOTUPTOJP',...
             sprintf('glmalphapotuptoJp-%g-%i_%i-%g-%g-%i-%i-%i-%i.mat',...
                 TH,min(L),max(L),rnew,rold,phi,theta,omega,J));
    elseif length(TH)==2
        fname=fullfile(getenv('IFILES'),'GLMALPHAPOTUPTOJP',...
             sprintf('glmalphapotuptoJp-%g_%g-%i_%i-%g-%g-%i-%i-%i-%i.mat',...
                 max(TH),min(TH),min(L),max(L),rnew,rold,phi,theta,omega,J));
    else
        error('Bad choice for TH')
    end
end
    

if exist(fname,'file')==2 
    load(fname)
    disp(sprintf('Loading %s',fname))
else
% Maybe this can be sped up by using the individual m blocks and not one
% big matrix
[G,V]=glmalphapotup(TH,L,rnew,rold);
Grot=sparse((max(L)+1)^2,J);
[V,isrt]=sort(V,2);
V=fliplr(V);
G=G(:,fliplr(isrt));
G=out2on(G(:,1:J),max(L));
[dems,dels,mz,lmcosi,mzi,mzo,bigm,bigl,rinm,ronm]=addmon(max(L));
%h=waitbar(0,sprintf('Rotating the first %d Slepian functions',J));
disp('Calculating rotations in parallel mode')
%try
%parpool
%end
parfor j=1:J
    lmcosip=lmcosi;
    lmcosip(:,3:4)=reshape(insert(G(:,j),0,mzi),2,length(dems))';
    lmcosirot=plm2rot(lmcosip,omega,-theta,-phi);
    %tempry=kindeks(plm2rot(lmcosi,omega,-theta,-phi),3:4);
    g=reshape(lmcosirot(:,3:4),1,2*length(dems));
    g=g(mzo);
    Grot(:,j)=g(:);%tempry(ronm); 
%    waitbar(j/J,h);
end
%delete(gcp('nocreate'))
%delete(h)
V=V(1:J);
% Transform back to addmout 
if toout
    Grot=Grot(rinm,:);
end
try
  % If you are running Matlab
  save(fname,'Grot','V','-v7.3')
catch
  % If you are running octave
  save(fname,'Grot','V') 
end 
    
    
end
varns={Grot,V};
varargout=varns(1:nargout);
