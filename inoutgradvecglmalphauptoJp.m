function varargout=inoutgradvecglmalphauptoJp(TH,Lin,Lout,rsat,rpotin,rpotout,phi,theta,omega,J)
% [Grot,V]=inoutgradvecglmalphauptoJp(TH,Lin,Lout,rsat,rpotin,rpotout,phi,theta,omega,J)
%
% Rotates the best J Slepian functions from inoutgradvecglmalphaup from the
% poles to a center location of your choice and allows you to also rotate
% them around the center location.
%
% INPUT:
%
% TH        Angular extent of the spherical cap (degrees)
% Lin       Maximum spherical-harmonic degree for inner sources field
%           or passband [Lintmin Lintmax]
% Lout      Maximum spherical-harmonic degree for outer sources field
% rsat      average satellite radial position
% rpotin    focusing radius for inner source scalar potential
% rpotout   focusing radius for outer source scalar potential
% phi       Longitude of the center (degrees) [default: 0]
% theta     Colatitude of the center (degrees) [default: 0]
% omega     Anticlockwise azimuthal rotation (degrees) [default: 0]
% J         The number of eigenfunctions that are being asked (and saved)
%
% OUTPUT:
%
% Grot      The unitary matrix of localization coefficients. 
%           First (Lin+1)^2 coefficients are for Elm, last (Lout+1)^2-1
%           coefficients are for Flm.
% V         The eigenvalues (conditioning values) 
%
% Last modified by plattner-at-alumni.ethz.ch, 05/31/2016

toout=1; % Should the result be in addmout?

if length(Lin)==2
    Lin=[min(Lin) max(Lin)];
end

defval('phi',0);
defval('theta',0);
defval('omega',0);


% This here is just to get a unique name for saving. This allows us to do
% much faster calculations later by just loading
if length(Lin)==1
    if length(TH)==1
    fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHAUPTOJP',...
             sprintf('inoutgradvecglmalphauptoJp-%g-%i-%i-%g-%g-%g-%g-%g-%g-%i.mat',...
                 TH,Lin,Lout,rsat,rpotin,rpotout,phi,theta,omega,J));
    elseif length(TH)==2
    fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHAUPTOJP',...
             sprintf('gradvecglmalphauptoJp-%g-%g-%i-%i-%g-%g-%g-%g-%g-%g-%i.mat',...
                 max(TH),min(TH),Lin,Lout,rsat,rpotin,rpotout,phi,theta,omega,J));    
    else
        error('Bad choice for TH')
    end
else
    if length(TH)==1
    fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHAUPTOJP',...
             sprintf('inoutgradvecglmalphauptoJp-%g-%i_%i-%i-%g-%g-%g-%g-%g-%g-%i.mat',...
                 TH,min(Lin),max(Lin),Lout,rsat,rpotin,rpotout,phi,theta,omega,J));
    elseif length(TH)==2
    fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHAUPTOJP',...
             sprintf('gradvecglmalphauptoJp-%g-%g-%i_%i-%i-%g-%g-%g-%g-%g-%g-%i.mat',...
                 max(TH),min(TH),min(Lin),max(Lin),Lout,rsat,rpotin,rpotout,phi,theta,omega,J));    
    else
        error('Bad choice for TH')
    end
end

% If we have it we can load it
if exist(fname,'file')==2 
    load(fname)
    disp(sprintf('Loading %s',fname))
else % we need to calculate it
    
% Maybe this can be sped up by using the individual m blocks and not one
% big matrix

% So first we load the north polar caps
[H,S]=inoutgradvecglmalphaup(TH,Lin,Lout,rsat,rpotin,rpotout);

% Make sure that it is sparse and sorted. This doesn't take long, even for
% large H
H=sparse(H);
% And build the sparse matrix of rotated coefficients
Hrot=sparse((max(Lin)+1)^2+(Lout+1)^2-1,J);
% This is Frederik's classic way of sorting in desceding order
[S,isrt]=sort(S,2);
S=fliplr(S);
H=H(:,fliplr(isrt));

parfor j=1:J
    % First get the coefficients for Elm and Flm
    Ecoef=H(1:(max(Lin)+1)^2,j);
    Fcoef=H((max(Lin)+1)^2+1:end,j);
    % The turn both of them into the lmcosi format because that is what
    % plm2rot eats
    Elmcosi= coef2lmcosi(Ecoef,1);
    % Make sure to include the L=0 part for the Flm
    Flmcosi=[0 0 0 0;fcoef2flmcosi(Fcoef,1)];
    Elmcosirot=plm2rot(Elmcosi,omega,-theta,-phi);
    Flmcosirot=plm2rot(Flmcosi,omega,-theta,-phi);
    % And turn them back into coefficients. Here we include if it should be
    % in addmout or addmon format
    Ecoefrot=lmcosi2coef(Elmcosirot,toout);
    % Remove the L=0 component from Flmcosi
    Fcoefrot=flmcosi2fcoef(Flmcosirot(2:end,:),toout);
    % And put them both into the big matrix Hrot
    Hrot(:,j)=[Ecoefrot;Fcoefrot];
end


S=S(1:J);

try
	save(fname,'Hrot','S','-v7.3')
catch
  save(fname,'Hrot','S') 
end
	

end
varns={Hrot,S};
varargout=varns(1:nargout);
