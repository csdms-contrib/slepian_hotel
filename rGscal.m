function rG=rGscal(coefs,theta,phi,rad,rplanet,onorout)
% rG=rGscal(coefs,theta,phi,rad,rplanet,onorout)
% 
% Evaluates any radial derivative functions (e.g. scalar Slepian functions)
% given by the provided coefficients at the provided points. This function 
% is quite efficient as it does not assemble the entire Ylm matrix but 
% iterates through the degrees to save time and memory  
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
% rG        Matrix of evaluated Slepian functions, size J x npoints
%
% See also ylm
%
% Last modified by plattner-at-alumni.ethz.ch, 04/07/2015

defval('onorout',0)

Lmax=sqrt(size(coefs,1))-1;

% We need the coefficients in addmout for this to work efficiently. If they
% are in addmon, transform them to addmout
if ~onorout
    [~,~,~,~,~,~,~,~,rinm]=addmon(Lmax);
    coefs=coefs(rinm,:);
end   

rG=zeros(size(coefs,2),length(theta));

% Make sure phi and rad are a row vectors
phi=phi(:)';
rad=rad(:)';

% Also, include the phase shift
phi=phi+pi;


% Iterate over the Ls
for L=0:Lmax
    % Setting the ms
    m=-L:L;
    % Calculating the Xlm
    X=xlm(L,abs(m),theta);
    % HERE IS THE UPWARD CONTINUATION PART:
    % Multiply each Xlm with the corresponding r factor:
    X=X.*repmat( (-L-1)/rplanet*(rad/rplanet).^(-L-2) ,length(m),1);        
    % Make the longitudinal phase: ones, sines or cosines, sqrt(2) or not 
    % The negative m is the cosine
    P=diag(sqrt(2-(m(:)==0)))*...
	     cos(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi))); 
    % The Ylm are the product of the X and the P
    % Now instead of storing each submatrix YL=X.*P for each degree in the 
    % large matrix Y, we multiply it with the coefficients and sum them up  
    rG=rG+coefs(L^2+1:(L+1)^2,:)'*(X.*P);    
end
    
