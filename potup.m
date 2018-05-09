function varargout = potup(coefin,rnew,rold,Lmax,inv,Lrange)
% coefout=potup(coefin,rsat,rEarth,Lmax,inv,Lrange)
%
% Upward-continues potential-field scalar spherical harmonic coefficients
% 
% Calculates cnew_{lm} = cold_{lm}*(rnew/rold)^(-l-1)
%
% INPUT:
%
% coefin    spherical harmonic coefficient vector or Matrix of point 
%           evaluations or matrix of coefficients in either addmout or 
%           addmon ordering for degrees 0 to Lmax. The lm cycling must be
%           in the columns. 
%           These are the coefficients for the scalar potential 
%           at altitude rold
% rsat      Satellite radius for coefficients of radial field 
% rEarth    Earth radius for coeficients of potential field 
% Lmax      maximum degree
% inv       up (0) or down (1)? (Handbook: A->0 , A^{-1}->1 )
%           up also means derivative, down means integral
% Lrange    all the L that are used if not all l,m combinations between 0 
%           and Lmax are included in the matrix (for example in sdwcapup)
%           WARNING: Use correct Ls. Otherwise you seriously mess things
%           up.
%
% OUTPUT:
%
% coefout   spherical harmonic coefficient vector for radial field at
%           altitude rnew. Same format as coefin
%
% Last modified by plattner-at-alumni.ethz.ch, 05/09/2018


defval('Lrange',[])
defval('inv',0)

if ~isempty(Lrange)
    %disp('WARNING: Using provided L range in scalupderivative')
    bigl=Lrange(:);
else
    [~,~,~,~,~,~,~,bigl]=addmon(Lmax);
end

if inv
    % Actually we build A^(-1)
    A=(rnew/rold).^(bigl+1);
else
    A=(rold/rnew).^(bigl+1);
end
coefout=nan(size(coefin));

if (size(coefin,1)==length(bigl))
 for i=1:size(coefin,2)     
    coefout(:,i)=A.*coefin(:,i);
 end  
else
 error('Wrong number of entries or not column vector')
end

varns={coefout};
varargout=varns(1:nargout);




