function varargout=kernelcppotup(Lmax,dom,rnew,rold,pars,ngl,rotb)
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
% Last modified by plattner-at-alumni.ethz.ch, 06/27/2018

defval('ngl',200)
defval('rotb',0)
defval('pars',[])

% Generic path name that I like
filoc=fullfile(getenv('IFILES'),'KERNELCPPOTUP');
if length(Lmax)==1
    Lwrite=[0,Lmax];
else
    Lwrite=Lmax;
end
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
      fnpl=sprintf('%s/WREG-%s-%i-%i-%g-%g.mat',filoc,dom,min(Lwrite),max(Lwrite),...
          rnew,rold);
      % For some of the special regions it makes sense to distinguish
      % It it gets rotb=1 here, it doesn't in LOCALIZATION
       if strcmp(dom,'antarctica') && rotb==1 
     fnpl=sprintf('%s/WREG-%s-%i-%i-%i-%g-%g.mat',filoc,dom,min(Lwrite),max(Lwrite),rotb,...
         rnew,rold);
       end
    end
else
    % If, instead of a string, we have closed form coordinates, then make a
    % hash from the coordinates and use it as the filename.
    if exist('octave_config_info')
      h=builtin('hash','sha1',dom);
    else
      h=hash(dom,'sha1');
    end
    fnpl=sprintf('%s/%s-%i-%i-%g-%g.mat',filoc,h,min(Lwrite),max(Lwrite),rnew,rold);  
end

% If the continued kernel already exists, load it.
if exist(fnpl,'file')==2 && ~isstr(ngl)      
    load(fnpl)
    disp(sprintf('%s loaded by KERNELCPPOTUP',fnpl))
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
    
    % Now calculate AKA' = (A(AK)')'
    disp('Multiplication with A')
    K=potup(K,rnew,rold,max(Lmax),0);% This is AK
    K=K';% This is (AK)'
    disp('Multiplication with Aprime')
    K=potup(K,rnew,rold,max(Lmax),0);% This is A(AK)'
    K=K';% This is (A(AK)')'=AKA';
    % In order to avoid numerical desymmetrification:
    K=(K+K')/2;
    if exist('octave_config_info')% If you are running octave   
      save(fnpl,'Lmax','dom','ngl','K')
    else
      save(fnpl,'Lmax','dom','ngl','K','-v7.3')
    end
end

varns={K};
varargout=varns(1:nargout);
    
    
