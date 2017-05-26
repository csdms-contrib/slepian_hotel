function coefout=outupderivative(coefin,rsat,rpot,Lmax,inv,Lrange)
% coefout=outupderivative(coefin,rsat,rpot,Lmax,inv,Lrange)
%
% Continuation from radial position rpot to radial position
% rsat and derivative for coefficients of the OUTER gradient vector
% spherical harmonics
%
% INPUT:
%
% coefin    spherical harmonic coefficient vector or Matrix of point 
%           evaluations or matrix of coefficients in either addmout or 
%           addmon ordering for degrees 1 to Lmax. The lm cycling must be
%           in the columns. 
%           For inv=0:           
%           These are the coefficients or point values for the scalar 
%           potential at altitude rpot. 
%           OR: for inv=1:
%           These are the coefficients or point values for the vector 
%           gradient at altitude rsat and we want the scalar potential 
%           coefficients at altitude rpot
% rsat      gradient vector radial position 
% rpot      scalar potential radial position 
% Lmax      maximum degree
% inv       derivative (up) (0) or "integral" (down) (1)? 
%           (Handbook: A->0 , A^{-1}->1 )
% Lrange    all the L that are used if not all l,m combinations between 0 
%           and Lmax are included in the matrix (for example in sdwcapup)
%           WARNING: Use correct Ls. Otherwise you seriously mess things
%           up.
% 
% OUTPUT:
%
% coefout   spherical harmonic coefficient vector for Flm field at
%           altitude rsat Same format as coefin.
%           OR, if inv=1:
%           scalar spherical-harmonic coefficient vector at altitude rpot
%
% See also vecupderivative
%
% Last modified by plattner-at-alumni.ethz.ch, 4/15/2015


defval('Lrange',[])
defval('inv',0)

if ~isempty(Lrange)
    bigl=Lrange(:);
else
    [~,~,~,~,~,~,~,bigl]=addmon(Lmax);
    % And remove l=0
    bigl=bigl(2:end);
end

if inv
    % Actually we build B^(-1)
    A=rpot*(rpot/rsat).^(bigl-1)./sqrt(bigl.*( 2*bigl+1 ));
else
    A=sqrt(bigl.*( 2*bigl+1 )).*(rsat/rpot).^(bigl-1)/rpot;
end
   
coefout=nan(size(coefin));
if (size(coefin,1)==length(bigl))
 for i=1:size(coefin,2)     
    coefout(:,i)=A.*coefin(:,i);
 end  
else
 error('Wrong number of entries or not column vector')
end
    


