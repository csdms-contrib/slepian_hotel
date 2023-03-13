function varargout=makelocalfield(data,rad,lat,lon,dom,Lmax,J,rplanet,avgsat,rotcoord,loadJ,savename,niter,weights,lambda,coef0,altcog)
% [coef,condi,dataweights,fnpl]=makelocalfield(data,rad,lat,lon,dom,Lmax,J,rplanet,avgsat,rotcoord,loadJ,savename,niter,weights,lambda,coef0)
%
% Calculates potential crustal field at radius rplanet from local satellite
% data within chosen region dom.
%
% Works for radial conponent OR vector data (figures out depending on the 
% input data).
%
% INPUT:
%
% Required:
%
% data      EITHER radial component data:
%              column vector of values for the points given as rad,cola,lon
%           OR vectorial data:
%              givpien as data{1}=rad component, data{2}=colat component,
%              data{3}=lon component, 
%              or data is a vector of length 3*length(rad) as [rad;cola;lon]
% rad       radial position of satellite (planet radius + altitude)
% cola,lon  colatitude/longitude positions of the data values (both as
%           column values), 0<=cola<=pi; 0<=lon<=2pi
% dom       Slepian concentration target region either given as string 
%           (one of the region names)
%           OR integer (polar cap opening angle in degrees)
%           OR two polar cap opening angles for the ring in between them
%           OR [lon lat] an ordered list defining a closed curve [degrees]
%           OR several regions to add up/subtract (ONLY FOR VECTORIAL DATA): 
%           struct with
%           dom.name    for name of the combined region
%           dom.parts   for the cell array of names of the parts, 
%                       or cap opening angles
%                       or [cap,lon,colat] for rotated caps
%           dom.sign    for adding or subtracting 
%           Example: dom.parts{1}='namerica'; dom.parts{2}='samerica';
%                    dom.sign=[1,1]; dom.name='americas';
%                    dom.name='weirdRing'
%                    dom.parts{1}=30; dom.parts{2}=[5,5,10]; dom.sign=[1,-1]
%                    subtracts the ring of cTH=5, clon=5, ccola=10 from the
%                    larger polar cap
% Lmax      Maximum spherical harmonic degree
% J         How many Slepian functions should be used to calculate the
%           solution? More means more sensitive to noise but higher spatial
%           resolution
% rplanet   Planet radius to which the solution should be calculated
%
% Recommended:
%
% avgsat    Average satellite radial position for which the Slepian functions
%           should be calculated. The crustal magnetic field is of course 
%           calculated taking the true radial positions into account.
%
% Optional, depending on what you want to do:
%
% rotcoord  This option is needed if you work with spherical caps or rings
%           that are not centered at the north pole. Provide center as
%           [longitude colatitude] IN DEGREES [0<=360 0<=180].
% loadJ     For polar caps/rings that are not centered at north pole: You
%           might only want to calculate and save the first "loadJ"
%           evaluated Slepian functions. Default: loadJ=J
% savename  To make future calculations faster: Povide a name to save the 
%           evaluated Slepian functions. Default: no saving
% niter     Numer of iterations for iteratively reweighted residual
% weights   Weights you want to use in least squares, it is a vector, 
%           same size as either number of data locations or 3*number of 
%           data locations if you want to weight each component: 
%           [wBr wBth wBph]
% lambda    Damping factor [default = 0]
% coef0     reference model for damping: 
%           min( ||A*coef - dat||^2 + lambda*||coef-coef0||^2 )
%           [default = zeros]
% DAMPING CURRENTLY ONLY FOR RADIAL-ONLY INVERSION
%
% OUTPUT:
%
% coef          spherical-harmonic coefficients of the potential field
%               in ADDMON format
% condi         Conditioning number of the square matrix of evaluated 
%               Slepian functions
% dataweights   Final weights used in the iteratively reweighted least 
%               squares solution
% fnpl          path to saved evaluated matrix file
%
% Last modified by plattner-at-alumni.ethz.ch,  04/23/2018

defval('avgsat',[])
defval('rotcoord',[0 0])
defval('loadJ',J)
defval('savename',[])
defval('niter',10)
defval('weights',[])
defval('lambda',0)
defval('coef0',[])
defval('altcog',1)

if length(data)==3
    data=[data{1};data{2};data{3}];
end

if isempty(avgsat)
    avgsat=mean(rad);
end

% If we are working with radial derivative data only:
if length(data)==length(rad)
    if length(Lmax)>1
        warning('Bandpass for only radial data not yet implemented')
    end
    % Polar caps/rings?
    if (~ischar(dom)) & (rotcoord(1)~=0 || rotcoord(2)~=0)
        [G,~]=glmalphapotuptoJp(dom,Lmax,avgsat,rplanet,...
            rotcoord(1),rotcoord(2),0,loadJ);        
    else
        if altcog
            [G,~]=glmalphapotup(dom,Lmax,avgsat,rplanet);
        else
            [G,V]=glmalpha(dom,Lmax);
            if ~ischar(dom)
                 [V,isrt]=sort(sum(V,1),'descend');
                 G=G(:,isrt);
            end
        end
    end
    
    % Now, if available, load the previously calculated evaluated 
    % Slepian functions up to J=loadJ. For this make the dataname:
    filoc=fullfile(getenv('IFILES'),'SLEPEVAL');
    
    if length(Lmax)==1
    if ischar(dom)
        fnpl=sprintf('%s/%s-rad-%s-L%i-%g-%g-Jmax%i.mat',filoc,...
        savename,dom,Lmax,avgsat,rplanet,loadJ);
    elseif iscell(dom)
        % We have a cell, and we expect the format to be {'greenland' 1.0}
        % where we have a region and want to add a buffer region around it.
        % However, if dom{2} turns out to be zero, we should ignore it.
        if dom{2}==0; h=dom{1}; else h=[dom{1} num2str(dom{2})]; end
        buf=dom{2};
        % However, if dom{2} turns out to be zero, we should ignore it.
        fnpl=sprintf('%s/%s-rad-%s-L%i-%g-%g-Jmax%i.mat',filoc,...
        savename,h,Lmax,avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)==1
        fnpl=sprintf('%s/%s-rad-TH%i-L%i-%g-%g-Jmax%i.mat',filoc,...
        savename,dom,Lmax,avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)==2
        fnpl=sprintf('%s/%s-rad-TH%i_%i-L%i-%g-%g-Jmax%i.mat',filoc,...
		     savename,dom(1),dom(2),Lmax,avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)>2
      try
	nam=hash(dom,'sha1');
      catch
  nam=builtin('hash','sha1',dom);
      end
      fnpl=sprintf('%s/%s-rad-%s-L%i-%g-%g-Jmax%i.mat',filoc,...
		   savename,nam,min(Lmax),avgsat,rplanet,loadJ);      
    end
    else
    if ischar(dom)
        fnpl=sprintf('%s/%s-rad-%s-L%i_%i-%g-%g-Jmax%i.mat',filoc,...
        savename,dom,min(Lmax),max(Lmax),avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)==1
        fnpl=sprintf('%s/%s-rad-TH%i-L%i_%i-%g-%g-Jmax%i.mat',filoc,...
        savename,dom,min(Lmax),max(Lmax),avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)==2
        fnpl=sprintf('%s/%s-rad-TH%i_%i-L%i_%i-%g-%g-Jmax%i.mat',filoc,...
		     savename,dom(1),dom(2),min(Lmax),max(Lmax),avgsat,rplanet,loadJ);
    elseif ~ischar(dom) && length(dom)>2
      try
	nam=hash(dom,'sha1');
      catch
  nam=builtin('hash','sha1',dom);
      end
       fnpl=sprintf('%s/%s-rad-%s-L%i_%i-%g-%g-Jmax%i.mat',filoc,...
		    savename,nam,min(Lmax),max(Lmax),avgsat,rplanet,loadJ);
    end
    end
    
    
    % .. and load it
    if exist(fnpl,'file')==2 && ~isempty(savename)
        load(fnpl)
        fprintf('%s loaded by LocalInnerField\n',fnpl)
    else
        % We need to calculate it  
        %MloadJ=rGscal(G(:,1:loadJ),cola,lon,rad,rplanet,1);
        MloadJ=potupxyz(G(:,1:loadJ),rplanet,lon,lat,rad,Lmax); 
        MloadJ = MloadJ';
        % Only addmout format is available in potupxyz

        % And save it if dataname is provided
        if ~isempty(savename) 
            % See if it's Octave or Matlab
            try
              % It's Matlab
              save(fnpl,'MloadJ','-v7.3')
            catch
              % It's octave
              save(fnpl,'MloadJ')
            end
        end
    end
    % Now fit the data with the evaluate Slepian functions using
    % iteratively reweighted least squares.
       
    % First: Only select the Slepian functions we want:
    M=MloadJ(1:J,:);
    
    if ~isempty(weights)
        % Here I am multiplying the matrix with the weights
        M=M.*repmat(weights(:)',J,1);
        % and here the right hand side.
        data=weights(:).*data;
    end    
        
    % Prepare damping if requested:
    lambda=abs(lambda);
    if isempty(coef0)
      slepcoef0 = zeros(J,1);
    else
      % transform the provided spherical-harmonic reference 
      % coefficients into Slepian reference coefficients
      slepcoef0 = (G(:,1:J))'*coef0;
    end

    % Start with a least squares step:
    MM=M*M' + lambda*eye(J);
    Md=M*data + lambda*slepcoef0;
    slepcoef=MM\Md;

    % And now iteratively reweighted residual calculation
    [slepcoef,dataweights]=itweighres(M,data,slepcoef,niter,lambda,slepcoef0);
    
    
    % Turn the coefficients into spherical-harmonic coefficients: 
    coef=G(:,1:J)*slepcoef;       
   
   
       
else
    error('Something is not right with the provided data')        
end

% Coefs are in ADDMOUT. Transform to ADDMON:
coef=out2on(coef,max(Lmax));

if nargout>2
    condi=cond(MM);
else
    condi=[];
end   


varns={coef,condi,dataweights,fnpl};
varargout=varns(1:nargout);

end

