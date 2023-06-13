function varargout=kernelcppotupNoise(Lmax,dom,rnew,rold,pars,ngl,rotb,Ninv,Ndetails)
% K=kernelcppotup(Lmax,dom,rnew,rold,pars,ngl,rotb)
%
% This function is designed for the potential field at satellite altitude  
% 
% Continuation-cognizant scalar Slepian kernel for the upward continuation
% and derivative described in scalupderivative from altitude rold to rnew.
% This is done by loading kernelcp and then upward continuing the kernel
% and save it.
%
% INPUT:
%
% Lmax       Maximum angular degree (bandwidth)
% dom        'patch'   spherical patch [default], with specs in 'pars'
%                      NOTE: better use GRUNBAUM / PLM2ROT in this case
%            'sqpatch' square patch with [thN thS phW phE] in 'pars'
%            'england', 'eurasia',  'namerica', 'australia', 'greenland'
%            'africa', 'samerica', 'amazon', 'orinoco', 'gpsnamerica',
%            'antarctica', 'alloceans', with specs in 'pars'
%            OR: [lon lat] an ordered list defining a closed curve [degrees] 
% rnew       radius for radial component (at satellite altitude)
% rold       radius for scalar potential (on surface)
% pars       [th0,ph0,thR] for 'patch'
%                 th0  Colatitude of the cap center, in radians
%                 ph0  Longitude of the cap center, in radians
%                 thR  Radius of the cap, in radians
%            N  splining smoothness for geographical regions [default: 10]
%            OR: the string with the name you want the result saved as
% ngl        The degree of the Gauss-Legendre integration [default: 200]            
% rotb       0 That's it, you got it [default: 0]
%            1 For, e.g., 'antarctica', if you were given rotated coordinates
%            to make the integration procedure work, this option makes
%            sure that the kernel matrix reflects this. If not, you have
%            to apply counterrotation after diagonalizing in LOCALIZATION.
% Ninv       An inverse noise covariance matrix (i.e. a weight) that you 
%              want to minimize against in the eigen decomposition
%
% OUTPUT:
%
% Klmlmp     The localization kernel whose eigenfunctions we want,
%            indexed as: degree  [0  1  1  1  2  2  2  2  2]
%                        order   [0  0 -1  1  0 -1  1 -2  2]
%            The function LOCALIZATION later reindexes this in LMCOSI
%            fashion. Note: you can use ADDMOUT and ADDMON to modify, and
%            see, e.g. PLOTSLEP and KLMLMP2ROT for some implementations
% XY         The outlines of the region into which you are localizing
% K1         An intermediate result useful when rotb=1, see KLMLMP2ROT
% K          An verification result useful when rotb=1, see KLMLMP2ROT
%
% See also SCALUPDERIVATIVE, GLMALPHAUP
%
% Last modified by charig-at-arizona.edu, 05/18/2023
% Last modified by plattner-at-alumni.ethz.ch, 06/27/2018

defval('ngl',200)
defval('rotb',0)
defval('pars',[])
defval('Ninv',diag(ones(1,(max(Lmax)+1)^2)))

% Generic path name that I like
filoc=fullfile(getenv('IFILES'),'KERNELCPPOTUPNOISE');
if length(Lmax)==1
    Lwrite=[0,Lmax];
else
    Lwrite=Lmax;
end

% % Make a hash from the matrix
% try
%    h2=hash(Ninv(1:10,1:10),'sha1');
% catch
%    h2=hash(Ninv,'sha1');
% end


if isstr(dom)
    switch dom
      % If the domain is a square patch
     case 'sqpatch'
      fnpl=sprintf('%s/%s-%i-%i-%i-%i-%i-%i-%g-%g.mat',filoc,dom,min(Lwrite),max(Lwrite),...
           round(pars(1)*180/pi),round(pars(2)*180/pi),...
           round(pars(3)*180/pi),round(pars(4)*180/pi),...
           rnew,rold);
      % If the domain is a spherical patch
     case 'patch'
      fnpl=sprintf('%s/%s-%i-%i-%i-%i-%i-%g-%g.mat',filoc,dom,min(Lwrite),max(Lwrite),...
           round(pars(1)*180/pi),round(pars(2)*180/pi),...
           round(pars(3)*180/pi),rnew,rold);
      % If the domain is a named region or a closed contour
     otherwise
      fnpl=sprintf('%s/WREG-%s-%i-%i-%g-%g-%s.mat',filoc,dom,min(Lwrite),max(Lwrite),...
          rnew,rold,Ndetails);
      % For some of the special regions it makes sense to distinguish
      % It it gets rotb=1 here, it doesn't in LOCALIZATION
       if strcmp(dom,'antarctica') && rotb==1 
     fnpl=sprintf('%s/WREG-%s-%i-%i-%i-%g-%g.mat',filoc,dom,min(Lwrite),max(Lwrite),rotb,...
         rnew,rold);
       end
    end
 elseif iscell(dom)
    % We have a cell, and we expect the format to be {'greenland' 1.0}
    % where we have a region and want to add a buffer region around it.
    % However, if dom{2} turns out to be zero, we should ignore it.
    if dom{2}==0; h=dom{1}; else h=[dom{1} num2str(dom{2})]; end
    buf=dom{2};
    % However, if dom{2} turns out to be zero, we should ignore it.
    fnpl=sprintf('%s/WREG-%s-%i-%i-%g-%g-%s.mat',filoc,h,min(Lwrite),max(Lwrite),...
          rnew,rold,h2);
else
    % If, instead of a string, we have closed form coordinates, then make a
    % hash from the coordinates and use it as the filename.
    try
      h=hash(dom,'sha1');
    catch
      h=builtin('hash','sha1',dom);
    end
    fnpl=sprintf('%s/%s-%i-%i-%g-%g.mat',filoc,h,min(Lwrite),max(Lwrite),rnew,rold);  
end

% If the continued kernel already exists, load it.
if exist(fnpl,'file')==2 && ~isstr(ngl)      
    load(fnpl)
    disp(sprintf('%s loaded by KERNELCPPOTUPNOISE',fnpl))
else
    % Otherwise obtain the uncontinued kernel and continue it.
    try
        K=kernelcp(max(Lmax),dom,pars,ngl,rotb);
    catch
        K=kernelc(max(Lmax),dom,pars,ngl,rotb);
    end
    % Set all entries of K for l<min(Lmax) to zero
    if length(Lmax)==2
        maxind=(min(Lmax)+1)^2;
        K(1:maxind,:)=zeros(maxind,size(K,2));
        K(:,1:maxind)=zeros(size(K,2),maxind);
    end
    
    % Calculate D N^-1 D
    K = K*Ninv*K;

    % % Now calculate AKA' = (A(AK)')'
    % disp('Multiplication with A')
    % K=potup(K,rnew,rold,max(Lmax),0);% This is AK
    % K=K';% This is (AK)'
    % disp('Multiplication with Aprime')
    % K=potup(K,rnew,rold,max(Lmax),0);% This is A(AK)'
    % K=K';% This is (A(AK)')'=AKA';
    % In order to avoid numerical desymmetrification:
    K=(K+K')/2;
    try
      save(fnpl,'Lmax','dom','ngl','K','-v7.3')
    catch
      save(fnpl,'Lmax','dom','ngl','K')
    end
end

varns={K};
varargout=varns(1:nargout);
    
    
