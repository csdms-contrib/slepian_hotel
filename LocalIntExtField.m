function varargout=LocalIntExtField(data,rad,cola,lon,dom,LmaxIn,LmaxOut,J,rplanet,router,avgsat,rotcoord,loadJ,savename,niter,weights)
% [coef,condi,dataweights,fnpl]=LocalIntExtField(data,rad,cola,lon,dom,LmaxIn,LmaxOut,J,rplanet,router,avgsat,rotcoord,loadJ,savename,niter,weights)
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
%              given as data{1}=rad component, data{2}=colat component,
%              data{3}=lon component, 
%              or data is a vector of length 3*length(rad) as [rad;cola;lon]
% rad       radial position of satellite (planet radius + altitude)
% cola,lon  colatitude/longitude positions of the data values (both as
%           column values), 0<=cola<=pi; 0<=lon<=2pi
% dom       Slepian concentration target region either given as string 
%           (one of the region names)
%           OR integer (polar cap opening angle)
%           OR two polar cap opening angles for the ring in between them
%           OR [lon lat] an ordered list defining a closed curve [degrees]
% LmaxIn    Maximum spherical harmonic degree for the internal source part
%           or passband [Lintmin Lintmax]
% LmaxOut   Maximum spherical harmonic degree for the external source part
% J         How many Slepian functions should be used to calculate the
%           solution? More means more sensitive to noise but higher spatial
%           resolution
% rplanet   Planet radius to which the inner source solution should be 
%           calculated
% router    Outer radius for which the outer source solution should be
%           calculated
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
%
% OUTPUT:
%
% coef          spherical-harmonic coefficients of the inner and outer 
%               potential field 
%               (the first (LmaxIn+1)^2 are for the inner field)
%               in ADDMON format
% condi         Conditioning number of the square matrix of evaluated 
%               Slepian functions
% dataweights   Final weights used in the iteratively rewweighted least 
%               squares solution
% fnpl          path to saved evaluated matrix file
%
% Last modified by plattner-at-alumni.ethz.ch,  7/14/2017

defval('avgsat',[])
defval('rotcoord',[0 0])
defval('loadJ',J)
defval('savename',[])
defval('niter',10)
defval('weights',[])

if length(data)==3
  useit(1)=~isempty(data{1})
  useit(2)=~isempty(data{2});
  useit(3)=~isempty(data{3});
  data=[data{1};data{2};data{3}];
end

if isempty(avgsat)
    avgsat=mean(rad);
end

% If we are working with radial derivative data only:
%if length(data)==length(rad)

%    error('Inner and Outer solution only for vector data (so far)')
    
    % Now the full gradient case: This is almost the same:
%elseif length(data)==3*length(rad)
    % Vectorial case  
    % First load or calculate the Slepian function coefficients
    if (~ischar(dom)) & (rotcoord(1)~=0 || rotcoord(2)~=0)        
        [H,~]=inoutgradvecglmalphauptoJp(dom,LmaxIn,LmaxOut,avgsat,...
            rplanet,router,...        
           rotcoord(1),rotcoord(2),0,loadJ);       
    else
        if length(LmaxIn)==2 & ischar(dom)
            error('Bandpass only for spherical caps at this point')
        end
        [H,~]=inoutgradvecglmalphaup(dom,LmaxIn,LmaxOut,avgsat,...
            rplanet,router);
    end
    
    % Now, if available, load the previously calculated evaluated Slepian
    % functions up to J=loadJ
    filoc=fullfile(getenv('IFILES'),'INOUTSLEPEVAL');
    if length(LmaxIn)==1
    if ischar(dom)
        fnpl=sprintf('%s/%s-vec-%s-Lin%i-Lout%i-%g-%g-%g-Jmax%i.mat',filoc,...
        savename,dom,LmaxIn,LmaxOut,avgsat,rplanet,router,loadJ);
    elseif ~ischar(dom) && length(dom)==1
        fnpl=sprintf('%s/%s-vec-TH%f-Lin%i-Lout%i-%g-%g-%g-Jmax%i.mat',filoc,...
        savename,dom,LmaxIn,LmaxOut,avgsat,rplanet,router,loadJ);
    elseif ~ischar(dom) && length(dom)==2
        fnpl=sprintf('%s/%s-vec-TH%f_%f-Lin%i-Lout%i-%g-%g-%g-Jmax%i.mat',filoc,...
		     savename,dom(1),dom(2),LmaxIn,LmaxOut,avgsat,rplanet,router,loadJ);
    elseif ~ischar(dom) && length(dom)>2
      try
	      nam=hash(dom,'sha1');
      catch
        nam=builtin('hash','sha1',dom);
      end
      fnpl=sprintf('%s/%s-vec-%s-Lin%i-Lout%i-%g-%g-%g-Jmax%i.mat',filoc,...
        savename,nam,LmaxIn,LmaxOut,avgsat,rplanet,router,loadJ);
    end
    else
      if ischar(dom)
         fnpl=sprintf('%s/%s-vec-%s-L%i_%i-%g-%g-Jmax%i.mat',filoc,...
         savename,dom,min(Lmax),max(Lmax),avgsat,rplanet,loadJ);
      elseif ~ischar(dom) && length(dom)==1
        fnpl=sprintf('%s/%s-vec-TH%f-Lin%i_%i-Lout%i-%g-%g-%g-Jmax%i.mat',filoc,...
		     savename,dom,min(LmaxIn),max(LmaxIn),LmaxOut,avgsat,rplanet,router,loadJ);
      elseif ~ischar(dom) && length(dom)==2
        fnpl=sprintf('%s/%s-vec-TH%f_%f-Lin%i_%i-Lout%i-%g-%g-%g-Jmax%i.mat',filoc,...
		     savename,dom(1),dom(2),min(LmaxIn),max(LmaxIn),LmaxOut,avgsat,rplanet,router,loadJ);
      elseif ~ischar(dom) && length(dom)>2
	try
	  nam=hash(dom,'sha1');
  catch
    nam=builtin('hash','sha1',dom);
	end
	fnpl=sprintf('%s/%s-vec-%s-L%i_%i-%g-%g-Jmax%i.mat',filoc,...
         savename,nam,min(Lmax),max(Lmax),avgsat,rplanet,loadJ);
      end
      
    end
    
    if exist(fnpl,'file')==2 && ~isempty(savename)
        load(fnpl)
        fprintf('%s loaded by LocalInnerField\n',fnpl)

    else
      % We need to calculate it  
      MloadJ=rGvecInOut(H(:,1:loadJ),max(LmaxIn),cola,lon,rad,rplanet,router,1);
      % If we only have some of the components, only use that part of M:
      MloadJ=MloadJ(:,repelem(useit,length(rad)));      
      % And save it if dataname is provided
      if ~isempty(savename)
        try
          save(fnpl,'MloadJ','-v7.3')
        catch
          save(fnpl,'MloadJ')
        end
      end                
    end 
    
    % Now fit the data with the evaluate Slepian functions using
    % iteratively reweighted least squares.
       
    % First: Only select the Slepian functions we want:
    M=MloadJ(1:J,:);
    
    if ~isempty(weights)
        % if you want to weight each point, 
        % instead of each component of each point:
        %t1=tic;
        if length(weights)~=length(data)
            weights=[weights(:);weights(:);weights(:)];
        end 
        
        % Here I am reweighing the data points
        % Here the matrix
        M=M.*repmat(weights',J,1);                
        % Here the right hand side
        data=weights.*data;
    end        
        
    % Start with a least squares step:
    MM=M*M';
    Md=M*data;
    slepcoef=MM\Md;
    
    % And now iteratively reweighted residual calculation
    [slepcoef,dataweights]=itweighres(M,data,slepcoef,niter);    
    
     % Turn the coefficients into spherical-harmonic coefficients: 
    coef=H(:,1:J)*slepcoef;    
       
%else
%    error('Something is not right with the provided data')        
%end

% Coefs are in ADDMOUT. Transform to ADDMON:
coefIn=coef(1:(max(LmaxIn)+1)^2);
coefOut=[0;coef((max(LmaxIn)+1)^2+1:end)];
LmaxOut=sqrt(length(coefOut)+1)-1;

coefIn=out2on(coefIn,max(LmaxIn));
coefOut=out2on(coefOut,max(LmaxOut));
coef=[coefIn;coefOut(2:end)];

if nargout>2
    condi=cond(MM);
else
    condi=[];
end   


varns={coef,condi,dataweights,fnpl};
varargout=varns(1:nargout);

end



