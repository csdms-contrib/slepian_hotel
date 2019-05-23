function varargout=gradvecglmalphaup(TH,L,rnew,rold,srt,anti)
% [G,V]=gradvecglmalphaup(TH,L,rnew,rold,srt,anti)
% 
% Construction of gradient vector Slepian functions
%
% INPUT:
% 
% TH    Region name (e.g. 'africa' etc),
%       OR opening angle of spherical cap,
%       OR Angles of two spherical caps and we want the ring between 
%       them [TH1 TH2]
%       OR [lon lat] an ordered list defining a closed curve [degrees]
%       OR several regions to add up/subtract: 
%          struct with
%          TH.name    for name of the combined region
%          TH.parts   for the cell array of names of the parts, 
%                     or cap opening angles
%                     or [cap,lon,colat] for rotated caps
%          TH.sign    for adding or subtracting 
%          Example: TH.parts{1}='namerica'; TH.parts{2}='samerica';
%                   TH.sign=[1,1]; TH.name='americas';
%                   TH.name='weirdRing'
%                   TH.parts{1}=30; TH.parts{2}=[5,5,10]; TH.sign=[1,-1]
%                   subtracts the ring of cTH=5, clon=5, ccola=10 from the
%                   larger polar cap
% L     Bandwidth (maximum angular degree), or passband (two degrees)
% rnew  Satellite altitude
% rold  planet radius
% srt   Should the Slepian functions be sorted? [default = 1 = yes]
% anti  The opposite of the region (everything outside)? [default = 0 = no]
%       Only for named regions or closed curves.
%
% OUTPUT:
%
% G     Matrix containing vector Slepian function coefficients for the Elm
%       vector spherical harmonics IN ADDMOUT FORMAT
% V     Eigenvalues (conditioning values) UNSORTED FOR REGULAR REGIONS
%
% Last modified by plattner-at-alumni.ethz.ch, 5/23/2019

defval('anti',0)
defval('srt',1)

% Have to check if struct is correctly set up
if isstruct(TH)
  if ~(ischar(TH.name)&iscell(TH.parts)&isvector(TH.sign) )
    error('Something wrong with the struct you used for combining regions')
    error('Need TH.name (string) and TH.parts (cell array) and TH.sign (vector of 1 and -1)')
  end
end

% Figure out if it's lowpass or bandpass
lp=length(L)==1;
bp=length(L)==2;
maxL=max(L);

% The spherical harmonic dimension
ldim=(L(2-lp)+1)^2-bp*L(1)^2;

% First check if already calculated
  if ~isstr(TH) && ~isstruct(TH) && length(TH)==1 % POLAR CAPS
    defval('sord',1) % SINGLE OR DOUBLE CAP
    if lp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUP',...
		     sprintf('gradvecglmalphaup-%g-%i-%g-%g.mat',TH,L,rnew,rold));
    elseif bp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUP',...
		     sprintf('gradvecglmalphablup-%g-%i-%i-%g-%g.mat',...
             TH,L(1),L(2),rnew,rold));
    else
      error('The degree range is either one or two numbers')       
    end
    
  elseif ~isstr(TH) && ~isstruct(TH) && length(TH)==2 % Ring between polar cap angles
    defval('sord',1) % SINGLE OR DOUBLE CAP
    if lp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUP',...
		     sprintf('gradvecglmalphaup-%g_%g-%i-%g-%g.mat',max(TH),...
             min(TH),L,rnew,rold));
    elseif bp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUP',...
		     sprintf('gradvecglmalphablup-%g_%g-%i-%i-%g-%g.mat',...
			     max(TH),min(TH),L(1),L(2),rnew,rold));
    else
      error('The degree range is either one or two numbers')       
    end

    % Initialize ordering matrices
    %MTAP=repmat(0,1,ldim);
    %IMTAP=repmat(0,1,ldim);        
  else % GEOGRAPHICAL REGIONS and XY REGIONS
    defval('sord',10) % SPLINING SMOOTHNESS
    % We'll put in a Shannon number based on the area only, not based on
    % an actual sum of the eigenvalues
    defval('J',ldim)
    % Note the next line, though we can change our minds
    %defval('J',ldim*spharea(TH))
    if isstr(TH) % Geographic (keep the string)
      h=TH;
    elseif isstruct(TH)
      h=TH.name;
    else % Coordinates (make a hash)
      if exist('octave_config_info')
	h=builtin('hash','sha1',TH);
      else
	h=hash(TH,'sha1');
      end
    end
    if lp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUP',...
		     sprintf('gradvecglmalphaup-%s-%i-%i-%g-%g-%i.mat',...
             h,L,J,rnew,rold,anti));
    elseif bp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUP',...
		     sprintf('gradvecglmalphablup-%s-%i-%i-%i-%g-%g-%i.mat',...
             h,L(1),L(2),J,rnew,rold,anti));
    else
     error('The degree range is either one or two numbers')       
    end
    defval('GM2AL',NaN) % If not, calculate order per taper
    defval('MTAP',NaN) % If not, calculate order per taper
    defval('IMTAP',NaN) % And rank ordering within that taper
    defval('xver',0) % For excessive verification of the geographical case
  end
  
if exist(fname,'file')==2
  load(fname)
  disp(sprintf('Loading %s',fname))
else
  % Initialize matrices
  % plattner-at-alumni.ethz.ch, 2/11/2015:
  % This is now moved to the individual cases because generic regions will
  % reqire full matrices and polar caps require sparse matrices  
  %G=repmat(0,(maxL+1)^2,ldim);
  %V=repmat(0,1,ldim);
  
   % Find row indices into G belonging to the orders
  [EM,EL,mz,blkm]=addmout(maxL);
  
  % Find increasing column index; that's how many belong to this order
  % alpha=cumsum([1 L+1 gamini(L:-1:1,2)]);
  % The middle bit is twice for every nonzero order missing
  % alpha=cumsum([1 L(2-lp)-bp*L(1)+1 ...
  %   		gamini(L(2-lp)-bp*(L(1)-1),bp*2*(L(1)-1)) ...
  %   		gamini(L(2-lp)-bp*(L(1)-1):-1:1,2)]);
  % This should be the same for L and [0 L]
  alpha=cumsum([1 L(2-lp)-bp*L(1)+1 ...
  		gamini(L(2-lp)-bp*(L(1)-1),bp*2*L(1)) ...
  		gamini(L(2-lp)-bp*L(1):-1:1,2)]);
      
  % For GEOGRAPHICAL REGIONS or XY REGIONS
  if isstr(TH) || isstruct(TH) || length(TH)>2
  % Initialize matrices
  %G=zeros((maxL+1)^2,ldim);
  %V=zeros(1,ldim);  
    if bp
      error('Bandpass geographical tapers are not ready yet')
    end

    if isstruct(TH)
      % Several named regions. We will add them up.
      Klmlmp=zeros((L+1)^2,(L+1)^2);
      for reg=1:length(TH.parts)
          % If the subregion is a named region
          if ischar(TH.parts{reg})
             Kreg=kernelepup(L,TH.parts{reg},rnew,rold);
          else
             % If the subregion is a polar cap
             if length(TH.parts{reg})==1
                 % North-polar cap
                 Kreg=kernelepupcap(L,TH.parts{reg},rnew,rold);
             else
                 % Cap that needs to be rotated
                 cTH=TH.parts{reg}(1);
                 rotlon=TH.parts{reg}(2);
                 rotcola=TH.parts{reg}(3);
                 Kreg=kernelepupcap(L,cTH,rnew,rold,[rotlon,rotcola]);
             end
          end
        Klmlmp=Klmlmp + TH.sign(reg)*Kreg;
        
      end
      
    else 
        % If it's not a struct, it's either a string or a list of
        % coordinates. In both cases, kernelepup takes care of it.
        Klmlmp=kernelepup(L,TH,rnew,rold,[],[],[],anti);
    end
    
    % if anti==1
    %   % Get the complimentary region
    %   Klmlmp=eye(size(Klmlmp))-Klmlmp;
    % end
    
    
    % Calculates the eigenfunctions/values for this localization problem
    [G,V]=eig(Klmlmp);
    [V,isrt]=sort(sum(real(V),1));
    V=fliplr(V);
    G=G(:,fliplr(isrt));
    
    [a,b,c,d,e,f,ems,els,R1,R2]=addmon(L);
    % This indexes the orders of G back as 0 -101 -2-1012 etc
    G=G(R1,:);
    % Check indexing
%     difer(els(R1)-EL,[],[],mesg)
%     difer(ems(R1)-EM,[],[],mesg)
    
    % Calculate Shannon number and compare this with the theory
    N=sum(V);
    G=G(:,1:J);
    V=V(1:J);
    
    if exist('octave_config_info')
    % If you are running octave
    	save(fname,'G','V','EL','EM','N')  
    else
    % If you are running Matlab
    	save(fname,'G','V','EL','EM','N','-v7.3')
    end 
    
    
  else
        % For AXISYMMETRIC REGIONS
        % Initialize matrices
        G=sparse((maxL+1)^2,ldim);%repmat(0,(maxL+1)^2,ldim);
        V=zeros(1,ldim);
        disp('Calculating in parallel mode')
        try
            parpool
        end

        parfor mm=1:maxL+1       
            m=mm-1;
            %[E,Vpp,Np,th,Cp]=gradvecsdwcapup(TH,L,m,0,-1,[],[],rnew,rold);
            [~,Vpp,~,~,Cp]=gradvecsdwcapup(TH,L,m,0,-1,[],[],rnew,rold);
            Vp{mm}=Vpp;
            C{mm}=Cp;
        end 
        % Distribute this at the right point in the huge matrix
        for m=0:maxL
        if m>0
            % Here you supply the negative orders
            G(EM==-m,alpha(2*m):alpha(2*m+1)-1)=C{m+1};
            V(alpha(2*m):alpha(2*m+1)-1)=Vp{m+1};
            %MTAP(alpha(2*m):alpha(2*m+1)-1)=-m;
            % It's all neatly ordered here, downgoing within every order
            %IMTAP(alpha(2*m):alpha(2*m+1)-1)=1:length(Vp{m+1});
        end
      % Duplicate for the positive order in case the region is axisymmetric
      G(EM==m,alpha(2*m+1):alpha(2*m+2)-1)=C{m+1};
      V(alpha(2*m+1):alpha(2*m+2)-1)=Vp{m+1};
      %MTAP(alpha(2*m+1):alpha(2*m+2)-1)=m;
      % It's all neatly ordered here, downgoing within every order
      %IMTAP(alpha(2*m+1):alpha(2*m+2)-1)=1:length(Vp{m+1});
        end
    
        if srt
            [V,isrt]=sort(V,'descend');
            G=G(:,isrt);
        end
        
        if exist('octave_config_info')
    	% If you are running octave
    		save(fname,'G','V','EL','EM')  
    	else
    	% If you are running Matlab
    		save(fname,'G','V','EL','EM','-v7.3')
    	end 
        
       %save(fname,'G','V','EL','EM')
  end
end

% Provide output
varns={G,V,EL,EM};
varargout=varns(1:nargout);
