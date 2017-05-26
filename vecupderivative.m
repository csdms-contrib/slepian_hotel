function varargout=vecupderivative(coefin,rnew,rold,Lmax,inv,Lrange)
% coefout=vecupderivative(coefin,rsat,rEarth,Lmax,inv,Lrange)
%
% Upward continue calculate derivative of scalar spherical harmonic
% expansion then expressed in coefficients of Elm: 
% calculate cnew_{lm} = -cold_{lm}*sqrt((l+1)*(2*l+1))*rold^(-1)*(rnew/rold)^(-l-2)
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
% coefout   spherical harmonic coefficient vector for Elm field at
%           altitude rnew. Same format as coefin
%
% See also: RAD2POTINVERSION
%
% Last modified by plattner-at-alumni.ethz.ch, 04/28/2015 

defval('Lrange',[])
defval('inv',0)

if ~ischar(coefin)
if ~isempty(Lrange)
    %disp('WARNING: Using provided L range in vecupderivative')
    bigl=Lrange(:);
else
    [~,~,~,~,~,~,~,bigl]=addmon(Lmax);
end
if inv
    % Actually we build B^(-1)
    B=-rold*(rnew/rold).^(bigl+2)./sqrt((bigl+1).*(2*bigl+1));
else
    B=-sqrt((bigl+1).*(2*bigl+1)).*(rold/rnew).^(bigl+2)/rold;
end
coefout=nan(size(coefin));

if (size(coefin,1)==length(bigl))
 for i=1:size(coefin,2)     
    coefout(:,i)=B.*coefin(:,i);
 end  
else
 error('Wrong number of entries or not column vector')
end

varns={coefout};
varargout=varns(1:nargout);

elseif strcmp(coefin,'demo1')
% Upward continuation and gradient for a Slepian function
index=2;
dom='africa';%30
Lmax=30;
rad_E=6371;
a=3000;
res=1;
plottype=2;

% Get the Slepian coefficients
[H,S]=gradvecglmalpha(dom,Lmax);
if ~isstr(dom)
    [S,isrt]=sort(S,2);
    S=fliplr(S);
    H=H(:,fliplr(isrt));
end
% And transform the spherical harmonic sequence into the ADDMON format   
H=out2on(H,Lmax);
GVSfcoef=H(:,index);
upGVSfcoef=vecupderivative(GVSfcoef,rad_E+a,rad_E,Lmax,0);

% Now evaluate on regular grid
lon=(0:res:360)*pi/180;%(-180:res:180)*pi/180;
lat=(1:res:180-1)*pi/180;
lonp=lon;
latp=lat;
[lon,lat]=meshgrid(lonp,latp);
lo=lon(:);
la=lat(:);
[E,theta,phi]=elm(Lmax,la,lo-pi,1);
E{1}=out2on(E{1},Lmax);
E{2}=out2on(E{2},Lmax);
E{3}=out2on(E{3},Lmax);
fvals{1}=reshape(E{1}'*GVSfcoef,length(latp),length(lonp));% rad
fvals{2}=reshape(E{2}'*GVSfcoef,length(latp),length(lonp));% th
fvals{3}=reshape(E{3}'*GVSfcoef,length(latp),length(lonp));% ph



fvalsup{1}=reshape(E{1}'*upGVSfcoef,length(latp),length(lonp));% rad
fvalsup{2}=reshape(E{2}'*upGVSfcoef,length(latp),length(lonp));% th
fvalsup{3}=reshape(E{3}'*upGVSfcoef,length(latp),length(lonp));% ph



% plot:
% Eplot=reshape(E{1}(10,:),length(latp),length(lonp));
% plotplm(Eplot,lonp,pi/2-latp,2)

% Plot the upward continued on top of the Earth surface functions
subplot(2,3,1)
plotplm(fvalsup{1},lonp,pi/2-latp,plottype)
kelicol(1)
ax=max(abs(caxis));
caxis([-ax ax]);
subplot(2,3,2)
plotplm(fvalsup{2},lonp,pi/2-latp,plottype)
kelicol(1)
ax=max(abs(caxis));
caxis([-ax ax]);
subplot(2,3,3)
plotplm(fvalsup{3},lonp,pi/2-latp,plottype)
kelicol(1)
ax=max(abs(caxis));
caxis([-ax ax]);

subplot(2,3,4)
plotplm(fvals{1},lonp,pi/2-latp,plottype)
kelicol(1)
ax=max(abs(caxis));
caxis([-ax ax]);
subplot(2,3,5)
plotplm(fvals{2},lonp,pi/2-latp,plottype)
kelicol(1)
ax=max(abs(caxis));
caxis([-ax ax]);
subplot(2,3,6)
plotplm(fvals{3},lonp,pi/2-latp,plottype)
kelicol(1)
ax=max(abs(caxis));
caxis([-ax ax]);


end



