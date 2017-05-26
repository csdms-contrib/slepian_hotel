function rG=rGvecInOut(coefs,Lin,theta,phi,rad,rplanet,router,onorout)
% rG=rGvecInOut(coefs,Lin,theta,phi,rad,rplanet,router,onorout)
%
% Evaluates any inner and outer combined gradient function 
% (e.g. gradient vector Slepian functions)
% given by the provided coefficients of Elm and Flm at the provided points. 
% This function is relatively efficient as it does not assemble the entire 
% Elm matrix but iterates through the degrees to save time and memory. 
% It does, however, assemble the entire Flm matrix because usually, 
% the outer sources have a much lower max degree
% 
% INPUT:
%
% coefs     [(Lin+1)^2+(Lout+1)^2-1] x J matrix whose columns are the 
%           (Slepian) coefficients for Elm and Flm
% Lin       Max degree for Elm coefficients
% theta     colatitude of the locations in radians [0<=cola<=pi]
% phi       longitude of the locations in radians [0<=lon<=2*pi]
% rad       radial location of the data points 
% rplanet   radius for which the inner source coefficients are defined
% router    radius for which the outer source coefficients are defined
% onorout   are the columns of "coefs" in ADDMON (0) or ADDMOUT (1) format?
%           default: 0
%
% OUTPUT:
%
% r         Matrix of evaluated Slepian functions, size J x 3*npoints
%
% See also elm, flm, rGvec
%
% Last modified by plattner-at-alumni.ethz.ch, 11/21/2015

defval('onorout',0)

% First calculate rG for inner source
rG=rGvec(coefs(1:(Lin+1)^2,:),theta,phi,rad,rplanet,onorout);


% Now outer rG
coefsOut=coefs((Lin+1)^2+1:end,:);
Lout=sqrt(size(coefsOut,1)+1)-1;
[~,~,~,~,~,~,~,Ls,rinm]=addmon(Lout);
% If coefs are in addmout, then we are good, if they are in addmon, then we
% need to switch them up
if ~onorout
    coefsOut=[zeros(1,size(coefsOut,2));coefsOut];
    coefsOut=coefsOut(rinm,:);
    coefsOut=coefsOut(2:end,:);
end

Flmcell=flm(Lout,theta,phi+pi,1);
Flm=[Flmcell{1}(2:end,:) Flmcell{2}(2:end,:) Flmcell{3}(2:end,:)];
% Now upward continue Flm from its outer radius location, to the data point
% locations
Ls=Ls(2:end);
% Turn rad into a matrix
radMat=repmat(rad(:)',length(Ls),3);
Lmat=repmat(Ls(:),1,3*length(rad));
B=1/router*sqrt(Lmat.*( 2*Lmat+1 )).*((radMat/router).^(Lmat-1));

%B=1/router*sqrt(repmat(Ls(:),1,3*length(rad)).*( 2*repmat(Ls(:),1,3*length(rad))+1 )).*((repmat(rad(:)',length(Ls),3)/router).^(repmat(Ls(:),1,3*length(rad))-1));
Flmup=B.*Flm;

rG = rG + coefsOut'*Flmup;





