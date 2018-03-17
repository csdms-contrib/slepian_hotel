function K=kernelepupcap(L,TH,rnew,rold,rotcoords)
% K=kernelepupcap(L,TH,rnew,rold,rotcoords)
%
% Build the full cap Kernel matrix and rotate it to the right center
% location.
%
% INPUT:
%
% L          maximum spherical-harmonic degree
% TH         cap semi-opening angle in degrees
% rnew       satellite altitude
% rold       planet radius
% rotcoords  [longitude, colatitude] of cap center, in degrees
%
% Last modified by plattner-at-alumni.ethz.ch, 3/17/2018

%% Initialize
defval('rotcoords',[])

% This may need to be done in a sparse way
K=zeros((L+1)^2);

theLs=0:L;
% These are the sizes of each block
Lblocksizes=2*theLs+1;
% These are the indices where the blocks start
Lblockstarts=[1,1+cumsum(Lblocksizes)];
Lblockstarts=Lblockstarts(1:end-1);

%% Build the kernel matrix
m=0;
[~,~,~,~,~,Km]=gradvecsdwcapup(TH,L,m,[],[],[],[],rnew,rold);
posn=Lblockstarts;
K(posn,posn)=Km;

for m=1:L
    % Load all the individual m kernels
    [~,~,~,~,~,Km]=gradvecsdwcapup(TH,L,m,[],[],[],[],rnew,rold);
    % The negative ones
    % We need to take away those blocks for m>L, therefore the (m+1:end) 
    posn=Lblockstarts(m+1:end) + 2*m-1;
    K(posn,posn)=Km;
    % The positive ones
    posn=Lblockstarts(m+1:end) + 2*m;
    K(posn,posn)=Km;
end

%% Rotate if necessary
if any(rotcoords)
    disp('First rotation')
    tic
    Krot=rotatematrix(K,0,-rotcoords(2),-rotcoords(1),L)'; 
    toc
    disp('Second rotation')
    tic
    K=rotatematrix(Krot,0,-rotcoords(2),-rotcoords(1),L)'; 
    toc
end

% Remove numerical asymmetries 
K=(K+K')/2;

end

function Kout=rotatematrix(Kin,alpha,beta,gamma,L)
parfor j=1:length(Kin)
    [dems,dels,mz,lmcosi,mzi,mzo,bigm,bigl,rinm,ronm]=addmon(L);
    lmcosi(:,3:4)=reshape(insert(Kin(:,j),0,mzi),2,length(dems))';
    lmcosirot=plm2rot(lmcosi,alpha,beta,gamma);
    h=reshape(lmcosirot(:,3:4),1,2*length(dems));
    h=h(mzo);
    Kout(:,j)=h(:);     
end
end


