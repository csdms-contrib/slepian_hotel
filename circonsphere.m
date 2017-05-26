function [lon,cola]=circonsphere(clon,ccola,diam,radius,res)
% [lon,cola]=circonsphere(clon,ccola,diam,radius)
%
% Returns the corrdinates of the outline of a circle on the sphere
%
% INPUT: 
%
% clon    longitude of the center (0<=clon<2*pi)
% ccola   colatitude of the center (0<=ccola<=pi)
% diam    diameter of the circle (in same units as radius)
% radius  radius of the sphere (in same units as diam)
% res     how many points on the cicrle? Default:100
%
% OUTPUT:
%
% lon     longitudinal values of the circle outline points
% cola    colatitudinal values of the circle outline points
%
% Last modified by plattner-at-alumni.ethz.ch, 02/12/2014

defval('res',100)

% Start with the center of the circle in xyz coordinates:
[xc,yc,zc]=sph2cart(clon,pi/2-ccola,radius);

% Projection of center of the circle onto the sphere 
% in cartesian coordinates
cs=[xc;yc;zc];
% Normalize
cs=cs/norm(cs);

% Generate two orthogonal vectors randomly:
q1=rand(3,1);
q1=q1-dot(q1,cs)*cs;
q1=q1/norm(q1);
q2=cross(cs,q1);
q2=q2/norm(q2);


% The parameterization for the circle
t=linspace(0,2*pi,res);

% The opening angle of the circle
alpha=diam/2/radius;

% Circle radius in cartesian coords
a=sin(alpha)*radius;
% The center of the circle in 3D (inside the ball)
center=cs*cos(alpha)*radius;

poscirc=center*ones(1,length(t))+q1*sin(t)*a + q2*cos(t)*a;

% Now transform back to spherical coordinates
[lon,ele,rad]=cart2sph(poscirc(1,:),poscirc(2,:),poscirc(3,:));
cola=pi/2-ele;

% Check if rad==radius everywhere if you want
