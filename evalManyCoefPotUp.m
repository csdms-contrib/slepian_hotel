function Geval = evalManyCoefPotUp(coefs,theta,phi,onorout)
% Geval = evalManyCoefPotUp(coefs,theta,phi,onorout)
%
% Evaluates any set of coefficients (e.g. scalar Slepian functions)
% at the provided points after Continuation. This function 
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
% Geval        Matrix of evaluated Slepian functions, size J x npoints
%
% See also ylm
%
% Last modified by plattner-at-alumni.ethz.ch, 05/09/2018

defval('onorout',0)

Lmax=sqrt(size(coefs,1))-1;

% We need the coefficients in addmout for this to work efficiently. If they
% are in addmon, transform them to addmout
if ~onorout
    [~,~,~,~,~,~,~,~,rinm]=addmon(Lmax);
    coefs=coefs(rinm,:);
end

Geval = zeros(size(coefs,2),length(theta));

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
    X=X.*repmat((rad/rplanet).^(-L-1) ,length(m),1);        
    % Make the longitudinal phase: ones, sines or cosines, sqrt(2) or not 
    % The negative m is the cosine
    P=diag(sqrt(2-(m(:)==0)))*...
	     cos(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi))); 
    % The Ylm are the product of the X and the P
    % Now instead of storing each submatrix YL=X.*P for each degree in the 
    % large matrix Y, we multiply it with the coefficients and sum them up  
    Geval=Geval+coefs(L^2+1:(L+1)^2,:)'*(X.*P);    
end



