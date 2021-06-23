function varargout=kernelepup(Lmax,dom,rnew,rold,pars,ngl,rotb,anti)
% K=kernelepup(Lmax,dom,rnew,rold,pars,ngl,rotb,anti)
%
% Continuation-cognizant scalar Slepian kernel for the upward continuation
% and derivative described in vecupderivative from altitude rold to rnew.
% This is done by loading kernelep and then upward continuing the kernel
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
% anti       Calculate for opposite of region (outside instead of inside)
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
% See also VECUPDERIVATIVE, GRADVECGLMALPHAUP, KERNELCPUP
%
% Last modified by plattner-at-alumni.ethz.ch, 07/14/2017

defval('ngl',200)
defval('rotb',0)
defval('pars',[])
defval('anti',0)

% Generic path name that I like
filoc=fullfile(getenv('IFILES'),'KERNELEPUP');
if isstr(dom)
    switch dom
      % If the domain is a square patch
     case 'sqpatch'
      fnpl=sprintf('%s/%s-%i-%i-%i-%i-%i-%g-%g-%i.mat',filoc,dom,Lmax,...
           round(pars(1)*180/pi),round(pars(2)*180/pi),...
           round(pars(3)*180/pi),round(pars(4)*180/pi),...
           rnew,rold,anti);
      % If the domain is a spherical patch
     case 'patch'
      fnpl=sprintf('%s/%s-%i-%i-%i-%i-%g-%g-%i.mat',filoc,dom,Lmax,...
           round(pars(1)*180/pi),round(pars(2)*180/pi),...
           round(pars(3)*180/pi),rnew,rold,anti);
      % If the domain is a named region or a closed contour
     otherwise
      fnpl=sprintf('%s/WREG-%s-%i-%g-%g-%i.mat',filoc,dom,Lmax,...
          rnew,rold,i);
      % For some of the special regions it makes sense to distinguish
      % It it gets rotb=1 here, it doesn't in LOCALIZATION
       if strcmp(dom,'antarctica') && rotb==1 
     fnpl=sprintf('%s/WREG-%s-%i-%i-%g-%g-%i.mat',filoc,dom,Lmax,rotb,...
         rnew,rold,i);
       end
    end
else
    % If, instead of a string, we have closed form coordinates, then make a
    % hash from the coordinates and use it as the filename.
    try
      h=hash(dom,'sha1');
    catch
      h=builtin('hash','sha1',dom);
    end
    fnpl=sprintf('%s/%s-%i-%g-%g.mat',filoc,h,Lmax,rnew,rold);  
end

% If the continued kernel already exists, load it.
if exist(fnpl,'file')==2 && ~isstr(ngl)      
    load(fnpl)
    disp(sprintf('%s loaded by KERNELEPUP',fnpl))
else
    % Otherwise obtain the uncontinued kernel and continue it.
    K=kernelep(Lmax,dom,pars,ngl,rotb);
    if anti
        K = eye(size(K))-K;
    end
    % Now calculate BKB' = (B(BK)')'
    disp('Multiplication with B')
    K=vecupderivative(K,rnew,rold,Lmax,0);% This is BK
    K=K';% This is (BK)'
    disp('Multiplication with Bprime')
    K=vecupderivative(K,rnew,rold,Lmax,0);% This is B(BK)'
    K=K';% This is (B(BK)')'=BKB';
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
    
    
    
    
    
    
    
    
