function varargout=gradvecglmalphauptoJp(TH,L,rnew,rold,phi,theta,omega,J)
% [Grot,V]=gradvecglmalphauptoJp(TH,L,rnew,rold,phi,theta,omega,J)
%
% Gradient vectorial version if glmalphauptoJp. Loads gradvecglmalphaup and rotates the first
% J Slepian functions only. 
%
% INPUT:
%
% TH        Angular extent of the spherical cap (degrees)
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
% Grot     The unitary matrix of localization coefficients
% V        The eigenvalues 
%
% Last modified by plattner-at-alumni.ethz.ch, 10/14/2015

toout=1; % Should the result be in addmout?

defval('phi',0);
defval('theta',0);
defval('omega',0);

if length(L)==1
    if length(TH)==1
    fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUPTOJP',...
             sprintf('gradvecglmalphauptoJp-%g-%i-%g-%g-%g-%g-%g-%i.mat',...
                 TH,L,rnew,rold,phi,theta,omega,J));
    elseif length(TH)==2
    fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUPTOJP',...
             sprintf('gradvecglmalphauptoJp-%g-%g-%i-%g-%g-%g-%g-%g-%i.mat',...
                 max(TH),min(TH),L,rnew,rold,phi,theta,omega,J));    
    else
        error('Bad choice for TH')
    end
else
    if length(TH)==1
    fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUPTOJP',...
             sprintf('gradvecglmalphauptoJp-%g-%i-%i-%g-%g-%g-%g-%g-%i.mat',...
                 TH,min(L),max(L),rnew,rold,phi,theta,omega,J));
    elseif length(TH)==2
    fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUPTOJP',...
             sprintf('gradvecglmalphauptoJp-%g-%g-%i-%i-%g-%g-%g-%g-%g-%i.mat',...
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
[H,S]=gradvecglmalphaup(TH,L,rnew,rold);
H=sparse(H);
Hrot=sparse((max(L)+1)^2,J);
[S,isrt]=sort(S,2);
S=fliplr(S);
H=H(:,fliplr(isrt));
%H=out2on(full(H(:,1:J)),L);
H=out2on(H(:,1:J),max(L));
[dems,dels,mz,lmcosi,mzi,mzo,bigm,bigl,rinm,ronm]=addmon(max(L));
%hh=waitbar(0,sprintf('Rotating the first %d Slepian functions',J));
%disp('Calculating rotations in parallel mode')
%try
%matlabpool open
%end

% Avoid overwriting of dlmb matrix in parallel mode:
%dlmb(max(L));

%parfor j=1:J
% For large problems the parallel version breaks down. Have not yet figured
% out why...
% For J moderately large this does not take too long 
for j=1:J 
    tstart=tic;
    lmcosip=lmcosi;
    lmcosip(:,3:4)=reshape(insert(H(:,j),0,mzi),2,length(dems))';
    lmcosirot=plm2rot(lmcosip,omega,-theta,-phi);
    h=reshape(lmcosirot(:,3:4),1,2*length(dems));
    h=h(mzo);
    Hrot(:,j)=h(:); 
    tend=toc(tstart);
    fprintf('done with %d out of %d. It took %d seconds\n',j,J,tend);
end

    
%delete(hh)
S=S(1:J);
% Transform back to addmout 
if toout
    Hrot=Hrot(rinm,:);
end

if exist('octave_config_info')
   % If you are running octave
   save(fname,'Hrot','S')
else 
   % If you are running Matlab
   save(fname,'Hrot','S','-v7.3')
end  

end
varns={Hrot,S};
varargout=varns(1:nargout);
