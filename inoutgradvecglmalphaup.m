function varargout=inoutgradvecglmalphaup(TH,Lin,Lout,rsat,rpotin,rpotout,srt)
% [G,V]=inoutgradvecglmalphaup(TH,Lin,Lout,rsat,rpotin,rpotout,srt)
% 
% Construction of gradient vector Slepian functions for inner and outer
% sources. Inner and outer can each have their own max spherical-harmonic
% degree and individual reference radius.
%
% INPUT:
% 
% TH        Region name (e.g. 'africa' etc), OR opening angle of spherical 
%           cap, 
%           OR Angles of two spherical caps and we want the ring between 
%           them [TH1 TH2]
%           OR [lon lat] an ordered list defining a closed curve [degrees]  
% Lin       Maximum spherical-harmonic degree for inner sources
%           or passband [Lintmin Lintmax]
% Lout      Maximum spherical-harmonic degree for outer sources
% rsat      Satellite reference altitude for gradient
% rpotin    Reference radius for inner sources potential
% rpotout   Reference radius for outer sources potential
% srt       Should the Slepian functions be sorted? [default = 1 = yes]
%
% OUTPUT:
%
% G     Matrix containing vector Slepian function coefficients for the Elm
%       vector spherical harmonics IN ADDMOUT FORMAT
% V     Eigenvalues (conditioning values) UNSORTED FOR REGULAR REGIONS
%
% Last modified by plattner-at-alumni.ethz.ch, 07/14/2017

defval('anti',0)
defval('srt',1)

% Figure out if lowpass or bandpass
lp=length(Lin)==1;
bp=length(Lin)==2;
maxLin=max(Lin);
if bp
    Lin=[min(Lin) max(Lin)];
end

% The spherical harmonic dimension of the internal part
ldim=(Lin(2-lp)+1)^2-bp*Lin(1)^2;


if ~isstr(TH) && length(TH)==1 % POLAR CAPS
    defval('sord',1) % SINGLE OR DOUBLE CAP
    if lp
        fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHAUP',...
		   sprintf('inoutgradvecglmalphaup-%g-%i-%i-%g-%g-%g.mat',...
           TH,Lin,Lout,rsat,rpotin,rpotout));    
    elseif bp
        fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHAUP',...
		   sprintf('inoutgradvecglmalphaup-%g-%i_%i-%i-%g-%g-%g.mat',...
           TH,min(Lin),max(Lin),Lout,rsat,rpotin,rpotout)); 
    else
        error('The degree range is either one or two numbers')   
    end
       
elseif ~isstr(TH) && length(TH)==2 % Ring between polar caps
  defval('sord',1) % SINGLE OR DOUBLE CAP   
  if lp
    fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHAUP',...
		   sprintf('inoutgradvecglmalphaup-%g_%g-%i-%i-%g-%g-%g.mat',max(TH),...
			   min(TH),Lin,Lout,rsat,rpotin,rpotout));
  elseif bp
    fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHAUP',...
		   sprintf('inoutgradvecglmalphaup-%g_%g-%i_%i-%i-%g-%g-%g.mat',max(TH),...
			   min(TH),Lin(1),Lin(2),Lout,rsat,rpotin,rpotout));
  else
    error('The degree range is either one or two numbers')       
  end

         
else % GEOGRAPHICAL REGIONS and XY REGIONS
    defval('sord',10) % SPLINING SMOOTHNESS
    % This is in case we give the region as a list of lon lat coordinates
    if ischar(TH) % Geographic (keep the string)
      h=TH;
    else % Coordinates (make a hash)
      if exist('octave_config_info')
	h=builtin('hash','sha1',TH);
      else
	h=hash(TH,'sha1');
      end
    end  
    if lp
        fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHAUP',...
		     sprintf('inoutgradvecglmalphaup-%s-%i-%i-%g-%g-%g.mat',...
             h,Lin,Lout,rsat,rpotin,rpotout));
    elseif bp
        fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHAUP',...
		     sprintf('inoutgradvecglmalphaup-%s-%i_%i-%i-%g-%g-%g.mat',...
             h,Lin(1),Lin(2),Lout,rsat,rpotin,rpotout));
    else
        error('The degree range is either one or two numbers')       
    end
end

if exist(fname,'file')==2
  load(fname)
  disp(sprintf('Loading %s',fname))
else
    
  % Find row indices into G belonging to the orders
  [EM,EL,mz,blkm]=addmout(maxLin); 
    
  % Find increasing column index; that's how many belong to this order
  %alpha=cumsum([1 Lin(2-lp)-bp*Lin(1)+1 ...
  %		gamini(Lin(2-lp)-bp*(Lin(1)-1),bp*2*Lin(1)) ...
  %		gamini(Lin(2-lp)-bp*Lin(1):-1:1,2)]);  
  
  
  % For GEOGRAPHICAL REGIONS or XY REGIONS
  if isstr(TH) || length(TH)>2  
      if bp
        error('Bandpass geographical tapers are not ready yet')
      end
      
     % warning('Something is fishy. Test and fix!')
    % Load directly upward-continued kernel for inner sources
    KEE=kernelepup(Lin,TH,rsat,rpotin);
    % Load directly upward-continued kernel for outer sources
    KFF=kernelfpup(Lout,TH,rsat,rpotout);
    % Now the mixed part:
    % First load non-continued kernel:
    KEF=kernelmixefp(Lin,Lout,TH);
    % Now upward continue this to the right individual altitudes and take
    % corresponding derivatives. We do this as in kernelfpup or kernelepup:
    % AKB'=(B(AK)')'
    % First do AK: The Elm part    
    KEF=vecupderivative(KEF,rsat,rpotin,Lin,0);
    % Now (AK)'
    KEF=KEF';
    % Now B(AK)', that's the Flm part
    KEF=outupderivative(KEF,rsat,rpotout,Lout,0);
    % Now AKB' = (B(AK)')'
    KEF=KEF';
    % Now assemble the whole matrix K
    K=[KEE KEF;KEF' KFF];
    % Now remove asymmetries due to numerical noise:
    fprintf('Numerical asymmetry: %d\n',norm(K-K'))
    K=(K+K')/2;
    % Now eigenvalue decomposition
    [G,V]=eig(K);
    [V,isrt]=sort(sum(real(V),1));
    V=fliplr(V);
    G=G(:,fliplr(isrt));
    % Reorder eigenvector coefficient ordering from ADDMON to
    % ADDMOUT
    [~,~,~,~,~,~,~,~,Rin,~]=addmon(Lin);
    [~,~,~,~,~,~,~,~,Rout,~]=addmon(Lout);
    Gin=G(1:(Lin+1)^2,:);
    % To do addmon to addmout: include the Lout=0 row
    Gout=[nan(1,size(G,2));G((Lin+1)^2+1:end,:)];
    % Resort them individually
    Gin=Gin(Rin,:);    
    Gout=Gout(Rout,:);
    % Remove the nan line for L=0
    Gout=Gout(2:end,:);
    % And put them together
    G=[Gin;Gout];  
    % And save:
    if exist('octave_config_info')
        % If you are running octave
        save(fname,'Lin','Lout','TH','G','V')         
    else
        save(fname,'Lin','Lout','TH','G','V','-v7.3') 
    end
else
    % And now the polar caps:
    
    % Plattner 5/11/2016:
    % For now I'm using the strategy: Set up the entire max(Lin) part but
    % the first min(Lin)^2+1 entries will be zero.    
    
    deME=addmout(maxLin);
    deMF=addmout(Lout);
    
    % In the alpha vector we need to include the correct sizes of the
    % matrices
    mvec=0:max(maxLin,Lout);
    if length(Lin)==2
        sizE=maxLin+1-max(mvec,min(Lin));
    else
        sizE=max(maxLin+1-mvec, zeros(size(mvec)));
    end
    sizF=max(Lout+1-max(mvec,1), zeros(size(mvec)));    
    
    
    % The vectors sizE and sizF contain the sizes of the matrices for the E
    % components and the F components. Of course we will have them both.
    siztot=sizE+sizF;

    % The entries in alpha show the beginning (alpha(i)) and end 
    % (alpha(i+1))-1 of each block i. The first block, for m=0, only occurs
    % once. then for m=-1 m=1, etc we always have the positive and negative
    % m. Therefore we start with index 1, then the size of the m=0 block, 
    % then all other sizes twice. That's what gamini does. Then of course 
    % we need the cumulative index to know where in the big matrix we want
    % to put things.
    
    if length(Lin==2)
        % We need to + and - m, but m=0 is only once
        sizindices=[1 gamini(2:max(maxLin,Lout)+1,2)];
        % This could be done in a direct way without the for loop, but it's not
        % expensive to do this here..
        alpha=[1;nan(length(sizindices)-1,1)];
        for i=2:length(alpha)
            alpha(i)=alpha(i-1)+siztot(sizindices(i-1));
        end
        % The last one has the same number of Ls as the 
        % previous to last one because of m, -m 
        alpha=[alpha;alpha(end)+(alpha(end)-alpha(end-1))];
    else
        alpha=cumsum([1 siztot(1) gamini(siztot(2:end),2) ]);
    end    
    
    
    % The redistribution is not even as complicated as it seems. The alpha
    % intervalls need to be exactly the size of the number of eigenvectors.
    % The right location with the ms is taken care of by EM==m 
    % (might need to make this faster?)
    
    % To simplify things and avoid mistakes: Treat the Flm coefficients the
    % same way as the Elm coefficients and just remove the L=0 part if
    % necessary.
    
    % Initialize matrices  
    if length(Lin)==2
        GE=sparse((maxLin+1)^2,(maxLin+1)^2 - (min(Lin))^2 + (Lout+1)^2-1);
        % Treat Flm like Elm and remove the zero part later
        GF=sparse((Lout+1)^2,  (maxLin+1)^2 - (min(Lin))^2 + (Lout+1)^2-1);     
    else
        GE=sparse((maxLin+1)^2,(Lin+1)^2 + (Lout+1)^2-1);
        % Treat Flm like Elm and remove the zero part later
        GF=sparse((Lout+1)^2,  (Lin+1)^2 + (Lout+1)^2-1); 
    end
    V=zeros(1,(maxLin+1)^2 - (min(Lin)+1)^2 + (Lout+1)^2-1);
    disp('Calculating in parallel mode')
    try
        matlabpool open
    end
    % Now the same as in gradvecglmalphaup: Calculate the individual
    % solutions for the m and put them back in the right place
    parfor mm=1:max(maxLin,Lout)+1       
    %for mm=1:max(maxLin,Lout)+1   
            m=mm-1;            
            [Vpp,Cp]=inoutgradvecsdwcapup(TH,Lin,Lout,m,rsat,rpotin,rpotout);
            Vp{mm}=Vpp;
            sizE=max(maxLin+1-m,zeros(size(m)));
            sizF=max(Lout+1-max(m,1),zeros(size(m)));
            CE{mm}=Cp(1:sizE,:);            
            CF{mm}=Cp(sizE+1:end,:);
    end 
    
    % To make distribution a bit simpler, add the L=0 row to CF.
    % This is only necessary for m=0.
    CF{1}=[zeros(1,size(CF{1},2)); CF{1}];
    
    % Distribute this at the right point in the huge matrix
    for m=0:max(maxLin,Lout)   
        if m>0
            % Here you supply the negative orders
            GE(deME==-m,alpha(2*m):alpha(2*m+1)-1)=CE{m+1};  
            GF(deMF==-m,alpha(2*m):alpha(2*m+1)-1)=CF{m+1};  
            V(alpha(2*m):alpha(2*m+1)-1)=Vp{m+1};
        end
        % Duplicate for the positive order in case the region is axisymmetric
            GE(deME==m,alpha(2*m+1):alpha(2*m+2)-1)=CE{m+1};     
            GF(deMF==m,alpha(2*m+1):alpha(2*m+2)-1)=CF{m+1};  
            V(alpha(2*m+1):alpha(2*m+2)-1)=Vp{m+1};   
    end           

    
    GF=GF(2:end,:);
    G=[GE;GF];
    
    if srt
        [V,isrt]=sort(V,'descend');
        % Now remove the L=0 part from GF. It should be a zero row anyhow
        % Uncomment this if you want to test it
        % fprintf('Norm of L=0 in GF is %g\n',norm(GF(1,:)))        
        G=G(:,isrt);
    end
        
    if exist('octave_config_info')
    	% Octave
    	save(fname,'G','V')
    else
    	% Matlab
    	save(fname,'G','V','-v7.3')
    end
    	
    
end
end
    
% Provide output
varns={G,V};
varargout=varns(1:nargout);
         
