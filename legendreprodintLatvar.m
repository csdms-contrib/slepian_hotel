function varargout=legendreprodintLatvar(L1,m1,L2,m2,x0,method,rplanet,rsatfun)
% in=legendreprodintLatvar(L1,m1,L2,m2,x0,method,rplanet,rsatfun)
%
% This function is designed for the potential field at satellite altitude
% with a satellite radial position that varies with latitude, given as
% a function handle rsatfun
% 
% Evaluates the integral of the product of two Schmidt semi-normalized
% real Legendre polynomials P_lm(x)P_l'm'(x)dx from x to 1, where that
% means using Matlab's LEGENDRE([],[],'sch'), see XLM, YLM, PLM etc.
%
% \int_{0}^{\theta_0}P_{L1,m1}(\cos(\theta))...
%    P_{L2,m2}(\cos(\theta))\sin(\theta)\,d\theta
% or indeed \int_{x0}^{1}P_{L1,m1}(x)P_{L2,m2}(x)\,dx
%
% The normalization is such that the integration amounts to
% (4-2*(m==0))/(2l+1) over the entire interval from -1 to 1,
% This normalizes the spherical harmonics to 4\pi/(2l+1).
% Note Schmidt contains the sqrt(2) multiplying Xlm.
%
% INPUT:
%
% L1,L2        Angular degrees of the polynomials, L1,L2>=0
% m1,m2        Angular orders of the polynomials, 0<=m<=L
% x0           Single point with lower integration limit
% method      'automatic' Using analytical formula if possible (default)
%             'dumb' Forcing usage of dumb semi-analytical formula
%             'gl'   Exact result by Gauss-Legendre integration, when m1~=m2
%             'paul' By Wigner expansion and the method of Paul (1978)
% rplanet     planetary radial position
% rsatfun     function handle for satellite radial position depending
%             on sin(latitude) or cos(colatitude)
% savename    unique name you want to use for this function handle to
%             save and restore the calculated kernels  
%            
%
% OUTPUT:
%
% in           The integrated product. 
%
% EXAMPLE:
%
% legendreprodint('demo1') Wigner recursion vs. Gauss-Legendre, L1=L2, m=0
% legendreprodint('demo2') Wigner recursion vs. Gauss-Legendre, L1~=L2, m=0
% legendreprodint('demo3') Dumb summation vs. Gauss-Legendre, L1=L2, m=0
% legendreprodint('demo4') Paul recursion vs. Gauss-Legendre
% legendreprodint('demo5') Verify some analytical formulas
%
% Last modified by plattner-at-princeton.edu, 05/24/2011
% Last modified by fjsimons-at-alum.mit.edu, 06/12/2015

if ~isstr(L1)
  defval('L1',1)
  defval('m1',0)
  defval('L2',2)
  defval('m2',0)
  defval('x0',0)
  defval('method','automatic')

  if length(x0)~=1 && ~strcmp(method,'paul')
    error('Not for multiple limits')
  end
  % Standard spherical harmonics restrictions using LEGENDRE or LIBBRECHT 
  if m1>L1 | m2>L2 | m1<0 | m2<0
    error('Positive order must be smaller or equal than degree')
  end


      % If 'gl' requested by the user or required by the problem
      % Using Gauss-Legendre integration from
      % http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
      % Formulate the integrand as an inline function, anonymous is better
      integrand=inline(sprintf(...
          ['rindeks((%s(x)/%g).^(-%i-1).*legendre(%i,x,''sch''),%i).*',...
           'rindeks((%s(x)/%g).^(-%i-1).*legendre(%i,x,''sch''),%i)'],...
           func2str(rsatfun),rplanet,L1,L1,m1+1,...
           func2str(rsatfun),rplanet,L2,L2,m2+1));
      % Calculate the Gauss-Legendre coefficients    
      % Watch out: multiply if m is not 0
      [w,xgl,nsel]=gausslegendrecof(max(L1+L2,...
                                        max(200*((m1*m2)~=0),...
                                            1000*mod(m1+m2,2))));
      % For l=1 and m=0 this is not even close to enough nodes
      % Calculate integral
      in=gausslegendre([x0 1],integrand,[w(:) xgl(:)]);
      % disp(sprintf('Gauss-Legendre with %i points',nsel))
     
end  
  varns={in};
  varargout=varns(1:nargout);
  
