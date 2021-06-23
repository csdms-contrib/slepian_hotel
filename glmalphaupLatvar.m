function varargout=glmalphaupLatvar(TH,L,rplanet,rsatfun,savename,srt)
% [G,V]=GLMALPHAUPLATVAR(TH,L,rplanet,rsatfun,savename,srt)
%
% This function is designed for the radial component of the gradient at
% satellite altitude  
%   
% UPWARD CONTINUED AND DERIVATIVE VERSION OF:
%
% Returns an (lm)X(alpha) matrix with spherical harmonic coefficients of
% the BANDLIMITED or PASSBAND Slepian functions of the SINGLE or DOUBLE
% polar cap, or of a geographical region of interest. Only in the
% geographical case are the eigenvalues automatically sorted; if not, the
% column dimension is always block-ordered by virtue of the
% construction. The matrix G is orthogonal, G'*G is the identity.
%
% Should put an option to save only the essentials up to a certain truncation
%
% INPUT:
%
% TH       Angular extent of the spherical cap, in degrees OR
%          'england', 'eurasia',  'namerica', 'australia', 'greenland'
%          'africa', 'samerica', 'amazon', 'orinoco', 'antarctica', 'alloceans':
%          OR: [lon lat] an ordered list defining a closed curve [degrees]
%          OR: Angles of two spherical caps and we want the ring between 
%              them [TH1 TH2] 
% L        Bandwidth (maximum angular degree), or passband (two degrees)
% rplanet    planetary radial position
% rsatfun    function handle for satellite radial position depending
%            on sin(latitude) or cos(colatitude)
% savename   unique name you want to use for this function handle to
%            save and restore the calculated kernels              
% srt      sorted output?
%
% OUTPUT:
%
% G        The unitary matrix of localization coefficients; note how
%          LOCALIZATION delivers these as LMCOSI arrays into PLM2XYZ
% V        The eigenvalues in this ordering (not automatically sorted)
%
% Note that using ADDMOUT you can get this back to block-diagonal form
% G=glmalpha; [a,b,c,bl]=addmout(18); imagesc(G(bl,:))
% difer(G'*G-eye(size(G)))
%
% SEE ALSO:
%
% GLMALPHAPTO, ADDMOUT, ADDMON, KERNELC, LOCALIZATION, GALPHA, DLMLMP
%
% With contributions from charig-at-princeton.edu, 10/10/2011
% Last modified by fjsimons-at-alum.mit.edu, 10/11/2011
% Last modified by plattner-at-alumni.ethz.ch, 5/26/2017 

% Should be able to update this to retain the rank order per m as well as
% the global ordering. Does this work for the whole-sphere? In that case,
% should really want G to be the identity - all though of course,
% anything works, too. You don't get necessarily the spherical harmonics
% back...

defval('TH',30)
defval('srt',1)
defval('L',18)
% This is only relevant for the axisymmetric cap
blox=0;
upco=0;
resc=0;
anti=0;
%defval('blox',0)
%defval('upco',0)
%defval('resc',0)
%defval('anti',0)

defval('mesg','GLMALPHA Check passed')
% Hold all messages
mesg=NaN;

% Figure out if it's lowpass or bandpass
lp=length(L)==1;
bp=length(L)==2;
maxL=max(L);

% The spherical harmonic dimension
ldim=(L(2-lp)+1)^2-bp*L(1)^2;

% Just get the file names here
if upco==0 && resc==0
  if ~isstr(TH) && length(TH)==1 % POLAR CAPS
    defval('sord',1) % SINGLE OR DOUBLE CAP
    if lp
      fname=fullfile(getenv('IFILES'),'GLMALPHAUPLATVAR',...
		     sprintf('glmalphaupLatvar-%g-%i-%g-%s-%i-%i.mat',TH,L,...
             rplanet,savename,sord,blox));
    elseif bp
      fname=fullfile(getenv('IFILES'),'GLMALPHAUPLATVAR',...
		     sprintf('glmalphablupLatvar-%g-%i-%i-%g-%s-%i-%i.mat',...
			     TH,L(1),L(2),rplanet,savename,sord,blox));
    else
      error('The degree range is either one or two numbers')       
    end

    % Initialize ordering matrices
    MTAP=repmat(0,1,ldim);
    IMTAP=repmat(0,1,ldim);
  elseif ~isstr(TH) && length(TH)==2 % Ring between polar cap angles
    defval('sord',1) % SINGLE OR DOUBLE CAP
    if lp
      fname=fullfile(getenv('IFILES'),'GLMALPHAUPLATVAR',...
		     sprintf('glmalphaupLatvar-%g_%g-%i-%g-%g-%i-%i.mat',max(TH),...
             min(TH),L,rplanet,savename,sord,blox));
    elseif bp
      fname=fullfile(getenv('IFILES'),'GLMALPHAUPLATVAR',...
		     sprintf('glmalphablupLatvar-%g_%g-%i-%i-%g-%g-%i-%i.mat',...
			     max(TH),min(TH),L(1),L(2),rplanet,savename,sord,blox));
    else
      error('The degree range is either one or two numbers')       
    end

    % Initialize ordering matrices
    MTAP=repmat(0,1,ldim);
    IMTAP=repmat(0,1,ldim);
    
  else % GEOGRAPHICAL REGIONS and XY REGIONS
    defval('sord',10) % SPLINING SMOOTHNESS
    % We'll put in a Shannon number based on the area only, not based on
    % an actual sum of the eigenvalues
    defval('J',ldim)
    % Note the next line, though we can change our minds
    defval('J',ldim*spharea(TH))
    if isstr(TH) % Geographic (keep the string)
      h=TH;
    else % Coordinates (make a hash)
      try
	h=hash(TH,'sha1');
      catch
	h=builtin('hash','sha1',TH);
      end
    end
    if lp
      fname=fullfile(getenv('IFILES'),'GLMALPHAUPLATVAR',...
		     sprintf('glmalphaupLatvar-%s-%i-%g-%s-%i.mat',h,L,rplanet,...
             savename,J));
    elseif bp
      fname=fullfile(getenv('IFILES'),'GLMALPHAUPLATVAR',...
		     sprintf('glmalphablup-%s-%i-%i-%g-%s-%i.mat',h,L(1),L(2),...
             rplanet,savename,J));
    else
     error('The degree range is either one or two numbers')       
    end
    defval('GM2AL',NaN) % If not, calculate order per taper
    defval('MTAP',NaN) % If not, calculate order per taper
    defval('IMTAP',NaN) % And rank ordering within that taper
    defval('xver',0) % For excessive verification of the geographical case
  end
else
  fname='neveravailable';
  defval('xver',0) % For excessive verification of the upco'd case
end

if anti==1
  % Update the file name to reflect the complimentarity of the region
  fname=sprintf('%s-1.mat',pref(fname)); 
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
  if isstr(TH) || length(TH)>2
    % Initialize matrices
    G=repmat(0,(maxL+1)^2,ldim);
    V=repmat(0,1,ldim);             
    if bp
      error('Bandpass geographical tapers are not ready yet')
    end
    % Calculates the localization kernel for this domain
    Klmlmp=kernelcpupLatvar(L,TH,rplanet,rsatfun,savename);
    
    if anti==1
      % Get the complimentary region
      Klmlmp=eye(size(Klmlmp))-Klmlmp;
    end

    % Calculates the eigenfunctions/values for this localization problem
    [G,V]=eig(Klmlmp);
    [V,isrt]=sort(sum(real(V),1));
    V=fliplr(V);
    G=G(:,fliplr(isrt));
    
    [a,b,c,d,e,f,ems,els,R1,R2]=addmon(L);
    % This indexes the orders of G back as 0 -101 -2-1012 etc
    G=G(R1,:);
    % Check indexing
    difer(els(R1)-EL,[],[],mesg)
    difer(ems(R1)-EM,[],[],mesg)
    
    % Calculate Shannon number and compare this with the theory
    N=sum(V);

    if lp
      % Is the Shannon number right? Need the area of the region
      difer(ldim*Klmlmp(1)-N,[],[],mesg)
    elseif bp
      difer(ldim*spharea(TH)-N,[],[],mesg)
    end


      
    % Here's it's going to be almost trivial to get the integral over the
    % region of the Slepian eigenfunctions, given that it is a linear
    % combination of a scaled row of Klmlmp
    %% sarea=G*Klmlmp(:,1);
    % Which you could verify in the space domain using galpha
    
    % Here I should save the actual eigenfunctions
    % defval('J',round(N))
    % Truncate to the smaller amount of eigenfunctions and -values
    G=G(:,1:J);
    V=V(1:J);
    try
      save(fname,'G','V','EL','EM','N','-v7.3')
    catch
      save(fname,'G','V','EL','EM','N')
    end
  else
    % For AXISYMMETRIC REGIONS
    % Initialize matrices
    G=sparse((maxL+1)^2,ldim);%repmat(0,(maxL+1)^2,ldim);
    V=repmat(0,1,ldim);
    if blox~=0 & blox~=1
      error('Specify valid block-sorting option ''blox''')
    end
    % For the SINGLE or DOUBLE POLAR CAPS
    disp('Calculating in parallel mode')
    %try
    %    parpool
    %end
%    for m=0:maxL
%        Vp{m+1}=[];C{m+1}=[];
%    end
    C=cell(maxL+1,1);
    Vp=cell(maxL+1,1);
%     for m=0:maxL
%         C{m+1}=nan(maxL-m);
%         V{m+1}
%     end
    %parfor m=0:maxL
     parfor mm=1:maxL+1
         m=mm-1;
      % Same results for +/- m; set nth=0 thus no sign correction!
%       if sord==1
% 	if lp
% 	  [E,Vg,th,C,T,Vp]=grunbaum(TH,L,m,0);
% 	elseif bp
	  % Note that the small-eigenvalue eigenfunctions might be
          % numerically degenerate and thus not as using Grunbaum - if
          % you need to compare, compare where the eigenvalues are "big"      
	  [E,Vpp,Np,th,Cp]=sdwcapupLatvar(TH,L,m,0,-1,[],[],rplanet,rsatfun,savename);
      Vp{mm}=Vpp;
      C{mm}=Cp;

% 	end
%       elseif sord==2
% 	if lp
% 	  [E,Vg,th,C,T,Vp]=grunbaum2(TH,L,m,0);
% 	elseif bp
% 	  error('Bandpass double-cap tapers not ready yet')
% 	end
%       else
% 	error('Specify single or double polar cap')
%       end
      
    
      if upco~=0
	if upco>0
	  % The upward continuation matrix
	  A=diag((1+upco).^[-(m:L)-1]);
	elseif upco<0
	  % The downward continuation matrix
	  A=diag((1+abs(upco)).^[(m:L)+1]);
	end
	
	% Comparisons with Grunbaum only make sense for lowpass
% 	if xver==1 & lp
% 	  % This should give the same result, more or less, less accurate 
% 	  if sord==1
% 	    [a,Vs,c,d,Cs,e,f,g,h,j,D]=sdwcap(TH,L,m,0,-1);
% 	  else
% 	    [a,Vs,c,Cs,e,f,D]=sdwcap2(TH,L,m,0,-1);
% 	  end
% 	  % This should give the eigenvalues again, which we'd had from
%           % orthocheck 
% 	  warning off
% 	  % Check difference integration and kernel eigenvalues
% 	  difer(Vp(:)-diag((C'*D*C)./(C'*C)),[],[],mesg)
% 	  % Check difference integration and diagonalization eigenvalues
% 	  difer(Vp(:)-Vs(:),[],[],mesg)
% 	  % Check difference between the eigenfunctions barring the sign
% 	  % and only wherever the eigenvalues amount to anything
% 	  difer(abs(Cs(:,Vp>1e-6))-abs(C(:,Vp>1e-6)),[],[],mesg)
% 	  warning on
% 	  Vc=diag((C'*A*C*diag(Vp)*C'*A*C));
% 	  Vp0=Vp;
% 	end
% 	
	% Upward continuation from 1 to 1+a or from 1+a to 1:
	% New eigenfunctions, same name
	C{mm}=A*C{mm};
	% Calculate new eigenvalues, same name
	[jk1,jk2,jk3,Vp{mm},nofa]=orthocheck(C{mm},[],TH/180*pi,m,sord,1);

	% Make sure these are sorted, since that's not automatically the case
	% [Vp,ind]=sort(Vp,'descend');
	% C=C(:,ind);
	% Current thinking is: do NOT resort, as you'll want to compare the
	% best at a=0 with whatever it becomes later!
	
% 	if xver==1
% 	  warning off
% 	  % Check difference integration eigenvalues and those from kernel
% 	  difer(Vp{mm}(:)-diag((C{mm}'*D*C{mm})./(C{mm}'*C{mm})),[],[],mesg)
% 	  warning on
% 	  % Check how many Vp>Vp0 
% 	  disp(sprintf('%i/%i eigenvalues greater',sum(Vp{mm}(:)>Vp0(:)), ...
% 		       length(Vp0)))
% 	  disp(sprintf('Shannon number greater: %i',sum(Vp{mm})>sum(Vp0)))
% 	end
	if resc==1
	  % Rescale these functions to have an integral to unity over the
	  % sphere; note: this doesn't make the set orthonormal of course
	  C{mm}=C{mm}*diag(1./nofa);
	end
      end

  end
    for m=0:maxL
      % Distribute this at the right point in the huge matrix
      if m>0
	% Here you supply the negative orders
	G(EM==-m,alpha(2*m):alpha(2*m+1)-1)=C{m+1};
	V(alpha(2*m):alpha(2*m+1)-1)=Vp{m+1};
	MTAP(alpha(2*m):alpha(2*m+1)-1)=-m;
	% It's all neatly ordered here, downgoing within every order
	IMTAP(alpha(2*m):alpha(2*m+1)-1)=1:length(Vp{m+1});
      end
      % Duplicate for the positive order in case the region is axisymmetric
      G(EM==m,alpha(2*m+1):alpha(2*m+2)-1)=C{m+1};
      V(alpha(2*m+1):alpha(2*m+2)-1)=Vp{m+1};
      MTAP(alpha(2*m+1):alpha(2*m+2)-1)=m;
      % It's all neatly ordered here, downgoing within every order
      IMTAP(alpha(2*m+1):alpha(2*m+2)-1)=1:length(Vp{m+1});
    end

    if srt
        [V,isrt]=sort(V,'descend');
        G=G(:,isrt);
    end
    
    % Calculate the Shannon number and compare it to the theory
    N=sum(V);
%    if upco==0
%      difer(N-ldim*sord*(1-cos(TH/180*pi))/2,[],[],mesg);
%    end
    
    % Compute the sum over all orders of the squared coefficients
    % Thus works when they have not been blocksorted yet. 
    GM2AL=repmat(0,ldim,maxL+1);
    for l=0:maxL
      b=(l-1+1)^2+1;
      e=(l+1)^2;
      GM2AL(:,l+1)=sum(G(b:e,:).^2,1)';
    end

    % Make sure that the sum over all degrees is 1 - but I forgot why
    difer(sum(GM2AL,2)-1,[],[],mesg)

    % This is not blockdiagonal, unless you make it thus
    if blox==1
      G=G(blkm,:);
      EM=EM(blkm);
      EL=EL(blkm);
    end
    if ~strcmp(fname,'neveravailable') 
      % Save the results if it isn't a geographical region
      % If the variable is HUGE you must use the -v7.3 flag, if not, you
      % can safely omit it and get more backwards compatibility
      try
    	% If you are running Matlab
        save(fname,'G','V','EL','EM','N','GM2AL','MTAP','IMTAP','-v7.3')
      catch
      % If you are running octave
        save(fname,'G','V','EL','EM','N','GM2AL','MTAP','IMTAP')
      end
    end
  end
end

% Provide output
varns={G,V,EL,EM,N,GM2AL,MTAP,IMTAP};
varargout=varns(1:nargout);
