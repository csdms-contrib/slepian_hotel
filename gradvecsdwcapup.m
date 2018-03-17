function varargout=gradvecsdwcapup(TH,L,m,nth,vcut,grd,method,rnew,rold)
% [E,V,N,th,C,K]=gradvecsdwcapup(TH,L,m,nth,vcut,grd,method,rnew,rold)
%
% Calculates the upward continued and derived polar cap kernel for the 
% gradient vector spherical harmonics Elm
%
% INPUT:
%
% TH         Angular extent of the spherical cap, in degrees
%            OR: Angles of two spherical caps and we want the ring between 
%            them [TH1 TH2] 
% L          Bandwidth (maximum spherical harmonic degree),
% m          Spherical harmonic order of the required data window, -l<m<l
% nth         Number of colatitudes to evaluate (none if 0) [default: 720]
% vcut        Cut-off eigenvalue [default= eps*10]
%             If vcut==-1 automatically recalculates everything
% grd         1 Colatitudes only; returns matrix E
%             2 Colatitude/Longitude; returns cell E
% method      'gl' for Gauss-Legendre via LEGENDREPRODINT [default]
%             'paul' for Paul via PAUL 
% rnew        radius for gradient (at satellite altitude)
% rold        radius for scalar potential (on surface)
%
% OUTPUT:
%
% E           Optimally concentrated tapers, expanded to space domain 
% V           Eigenvalues, sorted
% N           Sum of ALL the eigenvalues (Shannon number for ALL orders)
% C         
% K
%
% See also GRADVECSDWCUP, KERNELEPUP, GRADVECGLMALPHAUP
%
% Last modified by plattner-at-alumni.ethz.ch, 4/22/2015



defval('TH',40)
defval('L',18)
defval('m',0)
defval('nth',0)%,720)
defval('vcut',-1)%,eps*10)
defval('grd',1)
defval('method','gl')
defval('rint',[1 2])

% Defval for output 
defval('N',[]);

% Work with the absolute value of m
mor=m;
m=abs(m);

if(m>max(L))
  error('Order cannot exceed degree')
end



% Figure out if it's lowpass or bandpass
lp=length(L)==1;
bp=length(L)==2;
maxL=max(L);

% The spherical harmonic dimension
ldim=(L(2-lp)+1)^2-bp*L(1)^2;

% Filename of saved data
dirname=fullfile(getenv('IFILES'),'GRADVECSDWCAPUP');
if length(TH)==1
    if lp
      fnpl=fullfile(dirname,sprintf(...
          'GRADVECSDWCAPUP-%g-%i-%i-%i-%g-%g.mat',TH,L,nth,m,rnew,rold));
    elseif bp
      fnpl=fullfile(dirname,sprintf(...
          'GRADVECSDWBLUP-%g-%i-%i-%i-%i-%g-%g.mat',TH,L(1),L(2),nth,m,...
           rnew,rold));
    else
      error('The degree range should be either one or two numbers')
    end
elseif length(TH)==2
    if lp
      fnpl=fullfile(dirname,sprintf(...
          'GRADVECSDWUP-%g-%g-%i-%i-%i-%g-%g.mat',...
          max(TH),min(TH),L,nth,m,rnew,rold));
    elseif bp
      fnpl=fullfile(dirname,sprintf(...
          'GRADVECSDWBLUP-%g-%g-%i-%i-%i-%i-%g-%g.mat',...
          max(TH),min(TH),L(1),L(2),nth,m,...
          rnew,rold));
    else
      error('The degree range should be either one or two numbers')
    end
    
else
    error('Bad choice for TH')
end    
    

% In 7.0.0.19901 (R14), in very rare cases (*-180-0-0) this
% save was buggy, and could not be subsequently loaded.
% Thus used  ~(L==180 & m==0) as an extra condition, now gone.
% Used matzerofix to remediate this.
if exist(fnpl,'file')==2   %&& (vcut>0) & 1==3
  load(fnpl)
  disp(sprintf('%s loaded by GRADVECSDWCAPUP',fnpl))
else
    
  %% Load Bm from Kerneltancapm 
  if length(TH)==2
    TH=[max(TH) min(TH)];
  end
  for whichTH=1:length(TH)
      %[~,~,Bmst{whichTH}]=kerneltancapm(TH(whichTH),L,m);
      %BBm=kernelbm(TH(whichTH),maxL,abs(m));
      %if bp
      %   Bmst{whichTH}=BBm(min(L):max(L),min(L):max(L)); 
      %end      
      Bmst{whichTH}=kernelbm(TH(whichTH),L,abs(m));
      % Make it same size as 
      if m==0        
        if length(L)==1 & min(L)~=0
            B0=zeros(maxL+1-max(m,bp*min(L)),maxL+1-max(m,bp*min(L)));%zeros(ldim);
            B0(2:end,2:end)=Bmst{whichTH};
            Bmst{whichTH}=B0;         
        end            
      end  
  end
  if length(TH)==2
    Bm=Bmst{1}-Bmst{2};
  else
     Bm=Bmst{1}; 
  end
    
  % Get the full Shannon number
  N=ldim*(1-cos(TH(1)/180*pi))/2;
  
  % Approximate Nyquist degree is $\pi/d\theta$
  if nth < (maxL+1) & nth~=0
    error('Sample finer to avoid aliasing')
  end
  
  % Convert to radians
  TH=TH*pi/180;

  for whichK=1:length(TH) % Everything for both K until we get eigenfunctions
  Kst{whichK}=zeros(maxL+1-max(m,bp*min(L)),maxL+1-max(m,bp*min(L)));
  % Initialize kernel
  %K=zeros(maxL+1-max(m,bp*min(L)),maxL+1-max(m,bp*min(L)));

  % Construct kernel: integral of Ylm over the cap
  lmin=max(m,bp*min(L));
  
%   rtop=max(rint);
%   rbot=min(rint);

  switch method
   case 'gl'
    tic
    for lr=lmin:maxL
      for lc=lr:maxL
        % Orthonormalization of Ylm is to unity over unit sphere
        % When TH=180, K should be the identity matrix
        % The pi*(1+(m==0)) is the longitudinal integral
        % Note that this line ALSO would work with 'paul' but then it
        % would be very inefficient in not reusing the database
        Kst{whichK}(lr+1-lmin,lc+1-lmin)=...
            legendreprodint(lr,m,lc,m,cos(TH(whichK)),method)...
            *sqrt(2*lr+1)*sqrt(2*lc+1)/(4*pi)*pi*(1+(m==0));       
        % Symmetrize
        Kst{whichK}(lc+1-lmin,lr+1-lmin)=...
            Kst{whichK}(lr+1-lmin,lc+1-lmin);
      end
    end
    toc
   case 'paul'
    tic
    % Load all the wigner0j symbols up to level 2*maxL
    [~,C0,S0,L0]=zeroj(0,0,2*maxL);
    % Load all the wigner3j symbols up to level 2*maxL
    % If S3 is an empty you just got the symmetrical variety back
    [~,C3,S3,L3]=threej(0,0,2*maxL);
    
    % Get all the integrals of single functions
    Itab=paul(2*maxL,cos(TH(whichK)));

    % We won't load the wigner3j data base since we've symmetrized the
    % storage 
    for lr=lmin:maxL
      for lc=lr:maxL
        % The required degree range
        ELL=(max(abs(lr-lc),2*m):(lr+lc));
        % Note that the zero-bottom symbol needs to have the top row sum
        % even. This is never empty, I think, but if it is, see LEGENDREPRODINT
        ELL=ELL(~mod(ELL+lr+lc,2));
        lELL=length(ELL); ELL=ELL(:);
        
        % The actual wigner0j symbols needed
        w0=zeroj(lr,lc,ELL,L0,2,C0,S0); w0=w0(:);
        % The actual wigner3j symbols needed
        if isempty(S3)
          % You have the symmetrized version, proceed as in THREEJ
          % Find where they're being kept
          [CC,oddperm,phasefix]=wignersort(ELL,lr,lc,-2*m,m,m);         
          % Do the initial evaluation from the loaded variable
          wm=full(C3(CC));
          % Fix the phase
          wm(oddperm)=wm(oddperm).*phasefix;
          % Now fix the triangle condition violations
          wm=wm.*triangle(repmat(lr,lELL,1),repmat(lc,lELL,1),ELL);
          % Now fix the order violations
          wm=wm.*~[lr<abs(m) | lc<abs(m) | ELL(:)<abs(-2*m)];
        else
          % You have the full sparse database ready to pass on to THREEJ
          wm=threej(ELL,lr,lc,-2*m,m,m,L3,[],C3,S3); wm=wm(:);
        end

        % And this here is as in LEGENDREPRODINT
        Q=(-1)^(m+m)*(2*ELL+1).*wm.*w0*sqrt(2-(m==0))*sqrt(2-(m==0))...
          /sqrt(2-[(m+m)==0]);
        % Figure out the right indices accoding to the way the Itab is set up
        indices=ELL.*(ELL+1)/2+m+m+1;       
        
        % The Paul-Gaunt integral with extra adaptation to the YLM
        Kst{whichK}(lr+1-lmin,lc+1-lmin)=Q(:)'*Itab(indices,:)...
            *sqrt(2*lr+1)*sqrt(2*lc+1)/(4*pi)*pi*(1+(m==0));       
        % Symmetrize
        Kst{whichK}(lc+1-lmin,lr+1-lmin)=...
            Kst{whichK}(lr+1-lmin,lc+1-lmin);
      end
    end
  end
  toc
  end % Until here we have two K
  if length(TH)==2
    K=Kst{1}-Kst{2};
  else
      K=Kst{1};
  end

  
  %% Now calculate factor matrices for P and B to generate E
  %bigl=(m:maxL)';
  bigl=(lmin:maxL)';
  facP= sqrt((bigl+1)./(2*bigl+1));
  facB=-sqrt( bigl   ./(2*bigl+1));
  facPmat=spdiags(facP,0,length(bigl),length(bigl));
  facBmat=spdiags(facB,0,length(bigl),length(bigl));
  
  %% And now combine them correctly
  Pm=K;
  K=facPmat*Pm*facPmat' + facBmat*Bm*facBmat'; 
  K=full(K);
  
  
  %% Now upward continue the matrix
  % calculate BKB' = (B(BK)')'
  K=vecupderivative(K,rnew,rold,[],0,lmin:maxL);% This is BK
  K=K';% This is (BK)'
  disp('Multiplication with Bprime')
  K=vecupderivative(K,rnew,rold,[],0,lmin:maxL);% This is B(BK)'
  K=K';% This is (B(BK)')'=BKB';
  % In order to avoid numerical desymmetrification:
  K=(K+K')/2;
  
  
  %% Calculate eigenvalues and eigenvectors; C'*C=I
  [C,V]=eig(K);
  
%   if lp && vcut<=0
%     % Check the partial Shannon number
% %    difer(sum(diag(V))-indeks(nsubm(N,m,3,L),'end'))
%   elseif bp
%     disp('Still need to fix NSUBM for bandpass functions')
%   end

  % Fill in the missing degrees in case it was bandpass
  C=[zeros(lmin-m,size(C,2)) ; C];
  
  % Order eigenvalues and eigenfunctions downgoing
  [V,isrt]=sort(sum(V,1),'descend');
  % V=fliplr(V); C=C(:,fliplr(isrt));
  C=C(:,isrt);
%   unc=fliplr(unc); sdl=fliplr(sdl); com=fliplr(com); sdv=fliplr(sdv); 

  % Only return nonzero "useful" eigenvalues
  C=C(:,V>vcut); 
  V=V(V>vcut);

  Np=[];
  
  
  % Compute spatial functions, colatitudinal part only
  if nth~=0
    % Zonal functions only 
    if m==0
      % Make spatial functions
      % This is SDW (2005) equation (5.10) combined with the sqrt(2-dom) of
      % (5.12) already included!
      [E,th]=pl2th(C,nth,1);
      th=th*180/pi;
      nlon=2*nth-1;
    else
      % This is SDW (2005) equation (5.10) combined with the sqrt(2-dom) of
      % (5.12) already included!
      [E,nlon,lat]=plm2th(C,nth,m,1);
      th=linspace(0,180,size(E,1));
    end    
    % Make E start with a positive lobe and ajust C too
    % Don't take first sample as for m~=0 it is 0
    for index=1:size(E,2)
      C(:,index)=C(:,index)*sign(E(2,index));
      E(:,index)=E(:,index)*sign(E(2,index));
    end
  else
    E=0;
    th=0;
    nlon=0;
  end
  % In 7.0.0.19901 (R14), in very rare cases (3-180-0-0) this
  % save is buggy, and cannot be subsequently loaded
  % Output in degrees
  save(fnpl,...
       'E','V','L','TH','C','th','nlon','K')
end


varns={E,V,N,th,C,K};
varargout=varns(1:nargout);
  