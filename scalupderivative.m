function varargout=scalupderivative(coefin,rnew,rold,Lmax,inv,Lrange)
% coefout=scalupderivative(coefin,rsat,rEarth,Lmax,inv,Lrange)
%
% Upward continue calculate derivative of scalar spherical harmonic
% expansion: 
% calculate cnew_{lm} = cold_{lm}*(-l-1)*rold^(-1)*(rnew/rold)^(-l-2)
%
% INPUT:
%
% coefin    spherical harmonic coefficient vector or Matrix of point 
%           evaluations or matrix of coefficients in either addmout or 
%           addmon ordering for degrees 0 to Lmax. The lm cycling must be
%           in the columns. 
%           These are the coefficients or point values for the scalar 
%           potential at altitude rold
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
% See also: RAD2POTINVERSION
%
% Last modified by plattner-at-alumni.ethz.ch, 04/28/2015 

defval('Lrange',[])
defval('inv',0)

if ~ischar(coefin)
if ~isempty(Lrange)
    %disp('WARNING: Using provided L range in scalupderivative')
    bigl=Lrange(:);
else
    [~,~,~,~,~,~,~,bigl]=addmon(Lmax);
end
if inv
    % Actually we build A^(-1)
    A=rold*(rnew/rold).^(bigl+2)./(-bigl-1);
else
    A=(-bigl-1).*(rold/rnew).^(bigl+2)/rold;
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




elseif strcmp(coefin,'demo1')
% Upward continuation and radial derivative for a Slepian function
index=10;
dom='africa';
Lmax=40;
rad_E=6371;
a=300;
res=1;


% Get the Slepian coefficients
[G,V]=glmalpha(dom,Lmax);
if ~ischar(dom)
    [V,isrt]=sort(V,2);
    V=fliplr(V);
    G=G(:,fliplr(isrt));
end
% And transform the spherical harmonic sequence into the ADDMON format   
G=out2on(G,Lmax);
Sfcoef=G(:,index);
upradSfcoef=scalupderivative(Sfcoef,rad_E+a,rad_E,Lmax,0);

% Now transform to lmcosi for plotting
[dems,dels,mz,lmc,mzin]=addmon(Lmax);

coeflm=reshape(insert(Sfcoef,0,mzin),2,length(dems))';  
lmcosi=[dels dems coeflm];
upcoeflm=reshape(insert(upradSfcoef,0,mzin),2,length(dems))';  
uplmcosi=[dels dems upcoeflm];


% plot
subplot(1,2,1)
plotplm(lmcosi,[],[],2,res)
range=caxis;
newrange=[-max(abs(range)) max(abs(range))];
caxis(newrange)
kelicol(1);

subplot(1,2,2)
plotplm(uplmcosi,[],[],2,res)
range=caxis;
newrange2=[-max(abs(range)) max(abs(range))];
caxis(newrange2)
kelicol(1);

elseif strcmp(coefin,'demo2')
% Plot an co-co Slepian function
index=1;
dom='africa';
Lmax=40;
rad_E=6371;
a=300;
res=1;

% Get the Slepian coefficients
[G,V]=glmalphaup(dom,Lmax,rad_E+a,rad_E);
if ~ischar(dom)
    [V,isrt]=sort(V,2);
    V=fliplr(V);
    G=G(:,fliplr(isrt));
end
% And transform the spherical harmonic sequence into the ADDMON format   
G=out2on(G,Lmax);
Sfcoef=G(:,index);
upradSfcoef=scalupderivative(Sfcoef,rad_E+a,rad_E,Lmax,0);

% Now transform to lmcosi for plotting
[dems,dels,mz,lmc,mzin]=addmon(Lmax);

coeflm=reshape(insert(Sfcoef,0,mzin),2,length(dems))';  
lmcosi=[dels dems coeflm];
upcoeflm=reshape(insert(upradSfcoef,0,mzin),2,length(dems))';  
uplmcosi=[dels dems upcoeflm];


% plot
subplot(1,2,1)
plotplm(lmcosi,[],[],2,res)
range=caxis;
newrange=[-max(abs(range)) max(abs(range))];
caxis(newrange)
kelicol(1);

subplot(1,2,2)
plotplm(uplmcosi,[],[],2,res)
range=caxis;
newrange2=[-max(abs(range)) max(abs(range))];
caxis(newrange2)
kelicol(1);

end

