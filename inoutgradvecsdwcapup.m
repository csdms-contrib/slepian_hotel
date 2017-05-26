function [V,C]=inoutgradvecsdwcapup(TH,Lin,Lout,m,rsat,rpotin,rpotout)
% [V,C]=inoutgradvecsdwcapup(TH,Lin,Lout,m,rsat,rpotin,rpotout)
%
% INPUT:
%
% TH         Angular extent of the spherical cap, in degrees
%            OR: Angles of two spherical caps and we want the ring between 
%            them [TH1 TH2] 
% Lin        Maximum spherical-harmonic degree for inner sources
%            or passband [Lintmin Lintmax]
% Lout       Maximum spherical-harmonic degree for outer sources
% rsat       Satellite reference altitude for gradient
% rpotin     Reference radius for inner sources potential
% rpotout    Reference radius for outer sources potential
%
% OUTPUT:
%
% C     Matrix containing vector Slepian function coefficients for the Elm
%       and Flm vector spherical harmonics IN ADDMOUT FORMAT
% V     Eigenvalues (conditioning values)
%
% See also gradvecsdwcapup, inoutgradvecglmalphaup
%
% Last modified by plattner-at-alumni.ethz.ch, 5/11/2016

if length(Lout)>1
    error('No Bandpass for external field')
end


maxLin=max(Lin);
lp=length(Lin)==1;
bp=length(Lin)==2;

if bp
    Lin=[min(Lin) max(Lin)];
end

defval('nth',0)%,720)
defval('vcut',-1)%,eps*10)
defval('grd',1)
defval('method','gl')

% Work with the absolute value of m
mor=m;
m=abs(m);

if(m>max(maxLin,Lout))
  error('Order cannot exceed degree')
end


% Filename of saved data
dirname=fullfile(getenv('IFILES'),'INOUTGRADVECSDWCAPUP');
if length(TH)==1
    if lp    
        fnpl=fullfile(dirname,sprintf(...
          'INOUTGRADVECSDWCAPUP-%g-%i-%i-%i-%g-%g-%g.mat',...
          TH,Lin,Lout,m,rsat,rpotin,rpotout));
    elseif bp
      fnpl=fullfile(dirname,sprintf(...
          'INOUTGRADVECSDWCAPUP-%g-%i_%i-%i-%i-%g-%g-%g.mat',...
          TH,Lin(1),Lin(2),Lout,m,rsat,rpotin,rpotout));
    else
      error('The degree range should be either one or two numbers')
    end 
elseif length(TH)==2
    if lp
        fnpl=fullfile(dirname,sprintf(...
          'INOUTGRADVECSDWUP-%g-%g-%i-%i-%i-%g-%g-%g.mat',...
          max(TH),min(TH),Lin,Lout,m,rsat,rpotin,rpotout));
    elseif bp
        fnpl=fullfile(dirname,sprintf(...
          'INOUTGRADVECSDWUP-%g-%g-%i_%i-%i-%i-%g-%g-%g.mat',...
          max(TH),min(TH),Lin(1),Lin(2),Lout,m,rsat,rpotin,rpotout));
    else
      error('The degree range should be either one or two numbers')
    end 
else
    error('Bad choice for TH')
end    
    

if exist(fnpl,'file')==2   %&& (vcut>0) & 1==3
  load(fnpl)
  disp(sprintf('%s loaded by inoutgradvecsdwcapup.m',fnpl))
else
    % First load Pm and Bm for L=max(Lin,Lout)
    L=max(maxLin,Lout);    
    if length(TH)==2
        %Pm=kernelpm(max(TH),L,m)-kernelpm(min(TH),L,m);
        [~,~,~,~,~,~,~,~,~,~,Pm1]=sdwcap(max(TH),L,m,0,-1);
        [~,~,~,~,~,~,~,~,~,~,Pm2]=sdwcap(min(TH),L,m,0,-1);        
        Pm=Pm1-Pm2;
        clear Pm1
        clear Pm2
        Bm=kernelbm(max(TH),L,m)-kernelbm(min(TH),L,m);                
    else       
        %Pm=kernelpm(TH,L,m);
        [~,~,~,~,~,~,~,~,~,~,Pm]=sdwcap(max(TH),L,m,0,-1);
        Bm=kernelbm(TH,L,m);
    end    
             
    % Now take different matrix subsets for EE, FF, and EF  
    % If m>Lin or m>Lout, then some of these matrices will simply be empty. 
    
    % For EE
    PmEE=Pm(1:maxLin-m+1,1:maxLin-m+1);
    if m==0
        % In this case we need to add the m=0 zero rows and columns
        BmEE=zeros(size(PmEE));
        BmEE(2:end,2:end)=Bm(1:maxLin,1:maxLin);         
    else
        BmEE=Bm(1:maxLin-m+1,1:maxLin-m+1);
    end

    % For FF
    if m==0
        PmFF=Pm(2:Lout+1,2:Lout+1);
        BmFF=Bm(1:Lout,1:Lout);
    else
        PmFF=Pm(1:Lout-m+1,1:Lout-m+1);
        BmFF=Bm(1:Lout-m+1,1:Lout-m+1);
    end
    
    
    % For EF
    if m==0
        PmEF=Pm(1:maxLin+1,2:Lout+1);
        BmEF=zeros(maxLin+1,Lout);
        BmEF(2:end,:)=Bm(1:maxLin,1:Lout);
    else
        PmEF=Pm(1:maxLin-m+1,1:Lout-m+1);
        BmEF=Bm(1:maxLin-m+1,1:Lout-m+1);
    end        
    

    
    % Now calculate factor matrices for P and B to calculate E and F
    theLE=(m:maxLin)';
    theLF=(max(m,1):Lout)';
    
    facPE= sqrt( (theLE+1)./(2*theLE+1) );
    facPF= sqrt( (theLF  )./(2*theLF+1) );
    facBE=-sqrt( (theLE  )./(2*theLE+1) );
    facBF= sqrt( (theLF+1)./(2*theLF+1) );
    
    facPEmat=spdiags(facPE,0,length(facPE),length(facPE));
    facPFmat=spdiags(facPF,0,length(facPF),length(facPF));
    facBEmat=spdiags(facBE,0,length(facBE),length(facBE));
    facBFmat=spdiags(facBF,0,length(facBF),length(facBF));
    
    KEE=  facPEmat*PmEE*facPEmat + facBEmat*BmEE*facBEmat;
    KFF=  facPFmat*PmFF*facPFmat + facBFmat*BmFF*facBFmat;
    KEF=  facPEmat*PmEF*facPFmat + facBEmat*BmEF*facBFmat;
    
    
    % Now upward continue the matrices:
    % calculate BKB' = (B(BK)')'
    % First KEE  
    %lminin=max(Lin,0);
    %lminout=max(Lout,0);
    if m<=maxLin % Otherwise this matrix is empty
        KEE=vecupderivative(KEE,rsat,rpotin,maxLin,0,theLE);
        KEE=KEE';
        KEE=vecupderivative(KEE,rsat,rpotin,maxLin,0,theLE);
        KEE=KEE';
    end   
       
    if m<=Lout % Otherwise this matrix is empty
        KFF=outupderivative(KFF,rsat,rpotout,Lout,0,theLF);
        KFF=KFF';
        KFF=outupderivative(KFF,rsat,rpotout,Lout,0,theLF);
        KFF=KFF';
    end
    
    if m<=maxLin & m<=Lout % Otherwise this matrix is empty
        KEF=vecupderivative(KEF,rsat,rpotin,maxLin,0,theLE);
        KEF=KEF';
        KEF=outupderivative(KEF,rsat,rpotout,Lout,0,theLF);
        KEF=KEF';
    end
    
    % Plattner 5/11/2016: Now delete the low degrees
    if length(Lin)==2 & m<min(Lin);
        KEE=KEE(min(Lin)-m+1:end,min(Lin)-m+1:end);
        KEF=KEF(min(Lin)-m+1:end,:);
    end    
    
    % Now assemble the entire matrix
    K = [KEE KEF;KEF' KFF];           
    
    % And make sure it is symmetric
    try
        fprintf('Numerical asymmetry %g\n',norm(K-K'));
    catch
        fprintf('Numerical asymmetry %g\n',norm(full(K-K')));
    end
    K=(K+K')/2;
    
    % Calculate eigenvalues and eigenvectors; C'*C=I
    [C,V]=eig(full(K));
    
    [V,isrt]=sort(sum(V,1),'descend');
    
    C=C(:,isrt);
    
    % Plattner 5/11/2016: Now add zeros for the low degrees
    if length(Lin)==2 & m<min(Lin);
        C=[zeros(min(Lin)-m,size(C,2));C];
    end
    
    save(fnpl,'C','V');
    
end  
