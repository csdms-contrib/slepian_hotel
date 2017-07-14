function K=kernelfpup(Lmax,dom,rsat,rpot,pars,ngl,rotb)
% K=KERNELFP(Lmax,dom,rsat,rpot,pars,ngl,rotb)
%
% Upward continued kernel for the OUTER gradient vector Slepian functions
%
% INPUT:
%
% Lmax       Maximum spherical-harmonic degree (bandwidth)
% dom        'patch'   spherical patch [default], with specs in 'pars'
%                      NOTE: better use GRUNBAUM / PLM2ROT in this case
%            'sqpatch' square patch with [thN thS phW phE] in 'pars'
%            'england', 'eurasia',  'namerica', 'australia', 'greenland'
%            'africa', 'samerica', 'amazon', 'orinoco', 'gpsnamerica',
%            'antarctica', 'alloceans', with specs in 'pars'
%            OR: [lon lat] an ordered list defining a closed curve [degrees] 
% rsat       radius satellite altitude
% rpot       radius where to evaluate scalar potential
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
% K          The localization kernel whose eigenfunctions we want,
%            indexed as: degree  [0  1  1  1  2  2  2  2  2]
%                        order   [0  0 -1  1  0 -1  1 -2  2]
%
% See also kernelepup
%
% Last modified by plattner-at-alumni.ethz.ch, 7/14/2017

defval('ngl',200)
defval('rotb',0)
defval('pars',[])


% Generic path name that I like
filoc=fullfile(getenv('IFILES'),'KERNELFPUP');

if ischar(dom)
    switch dom
      % If the domain is a square patch
     case 'sqpatch'
      fnpl=sprintf('%s/Fcoef_%s-%i-%i-%i-%i-%i-%g-%g.mat',filoc,dom,Lmax,...
           round(pars(1)*180/pi),round(pars(2)*180/pi),...
           round(pars(3)*180/pi),round(pars(4)*180/pi),...
           rsat,rpot);
      % If the domain is a spherical patch
     case 'patch'
      fnpl=sprintf('%s/Fcoef_%s-%i-%i-%i-%i-%g-%g.mat',filoc,dom,Lmax,...
           round(pars(1)*180/pi),round(pars(2)*180/pi),...
           round(pars(3)*180/pi),rsat,rpot);
      % If the domain is a named region or a closed contour
     otherwise
      fnpl=sprintf('%s/Fcoef_WREG-%s-%i-%g-%g.mat',filoc,dom,Lmax,...
          rsat,rpot);
      % For some of the special regions it makes sense to distinguish
      % It it gets rotb=1 here, it doesn't in LOCALIZATION
       if strcmp(dom,'antarctica') && rotb==1 
     fnpl=sprintf('%s/Fcoef_WREG-%s-%i-%i-%g-%g.mat',filoc,dom,Lmax,rotb,...
         rsat,rpot);
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
    fnpl=sprintf('%s/Fcoef_%s-%i-%g-%g.mat',filoc,h,Lmax,rsat,rpot);  
end

% If the continued kernel already exists, load it.
if exist(fnpl,'file')==2 && ~isstr(ngl)      
    load(fnpl)
    disp(sprintf('%s loaded by KERNELFPUP',fnpl))
else
    % Otherwise obtain the uncontinued kernel and continue it.
    K=kernelfp(Lmax,dom,pars,ngl,rotb);
    % Now calculate AKA' = (A(AK)')'
    disp('Multiplication with Acirc')
    K=outupderivative(K,rsat,rpot,Lmax,0);% This is AK
    K=K';% This is (AK)'
    disp('Multiplication with Acircprime')
    K=outupderivative(K,rsat,rpot,Lmax,0);% This is B(BK)'
    K=K';% This is (A(AK)')'=AKA';
    % In order to avoid numerical desymmetrification:
    K=(K+K')/2;
    if exist('octave_config_info')% If you are running octave
      save(fnpl,'Lmax','dom','ngl','K')
    else
      save(fnpl,'Lmax','dom','ngl','K','-v7.3')
    end
end
    

    
    
