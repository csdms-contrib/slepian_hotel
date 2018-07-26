function rG=rGvec(coefs,theta,phi,rad,rplanet,onorout)
% rG=rGvec(coefs,theta,phi,rad,rplanet,onorout)
%
% Evaluates any gradient functions (e.g. gradient vector Slepian functions)
% given by the provided coefficients of Elm at the provided points. 
% This function is quite efficient as it does not assemble the entire Elm
% matrix but iterates through the degrees to save time and memory  
% 
% INPUT:
%
% coefs     (L+1)^2 x J matrix whose columns are the (Slepian) coefficients 
% theta     colatitude of the locations in radians [0<=theta<=pi]
% phi       longitude of the locations in radians [0<=phi<=2*pi]
% rad       radial location of the data points 
% rplanet   radius for which the coefficients are defined
% onorout   are the columns of "coefs" in ADDMON (0) or ADDMOUT (1) format?
%           default: 0
%
% OUTPUT:
%
% r         Matrix of evaluated Slepian functions, size J x 3*npoints
%           Order: radial, colatitudinal, longitudinal
%
% See also elm, rGscal
%
% Last modified by plattner-at-alumni.ethz.ch, 04/07/2016
% Clarified output: plattner-at-alumni.ethz.ch, 14/07/2018

defval('onorout',0)

Lmax=sqrt(size(coefs,1))-1;

% We need the coefficients in addmout for this to work efficiently. If they
% are in addmon, transform them to addmout
if ~onorout
    [~,~,~,~,~,~,~,~,rinm]=addmon(Lmax);
    coefs=coefs(rinm,:);
end   

rG=zeros(size(coefs,2),3*length(theta));

% Make sure phi and rad are a row vectors
phi=phi(:)';
rad=rad(:)';

% Also, include the phase shift
phi=phi+pi;

divsinvals=1./sin(theta(:)');

% Start with L=0 because it is special. for L=0 Elm is equal to Ylm:
L=0;
m=0;
% Calculating the Xlm
X=xlm(L,abs(m),theta);
% HERE IS THE UPWARD CONTINUATION PART:
% Multiply each Xlm with the corresponding r factor:
X=X.*( (-L-1)/rplanet*(rad/rplanet).^(-L-2) );        
% Make the longitudinal phase: ones, sines or cosines, sqrt(2) or not 
% The negative m is the cosine
P=diag(sqrt(2-(m(:)==0)))*...
     cos(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi))); 
rG(:,1:length(theta))=coefs(L^2+1:(L+1)^2,:)'*(X.*P); 
 

% Iterate over the other Ls
for L=1:Lmax
    % Setting the ms
    m=-L:L;    
    % Calculating the Xlm and the dXlm
    [X,dX]=xdxlm(L,abs(m),theta);
    % Make the longitudinal phase: ones, sines or cosines, sqrt(2) or not 
    % The negative m is the cosine
    P=diag(sqrt(2-(m(:)==0)))*...
         cos(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi)));
    dP=diag(sqrt(2-(m(:)==0)))*...
	   diag(m)*(-sin(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi))));
    % radial factor:
    Rfac=repmat( 1/rplanet*(rad/rplanet).^(-L-2) ,length(m),1);    
    % Formula for renormalized E 
    % (see Plattner and Simons (2015) Handbook eq. 22):
    % Erad   = (-l-1)*X*P
    % Etheta = dX*P
    % Ephi   = 1/sin(theta)*X*dP
    % And of course don't forget the radial factors.
    % Do the same thing as fin rGscal: sum up in each degree:
    % Radial component
    rG(:,1:length(theta)) =  rG(:,1:length(theta)) + ...
        coefs(L^2+1:(L+1)^2,:)'*( (-L-1)*Rfac.*X.*P );
    % Colatitudinal component
    rG(:,length(theta)+1:2*length(theta)) = ...
        rG(:,length(theta)+1:2*length(theta)) + ...
        coefs(L^2+1:(L+1)^2,:)'*( Rfac.*dX.*P );
    % Longitudinal component
    rG(:,2*length(theta)+1:end) = ...
        rG(:,2*length(theta)+1:end) + ... 
         coefs(L^2+1:(L+1)^2,:)'* ...
         ( Rfac.*repmat(divsinvals,length(m),1).*X.*dP );
end

    

