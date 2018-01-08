function varargout=GMTexportfield2(coef,fname,planetrad,alt,cmpwrite,area,cutcirc,nanrng,res,onorout)
% [datB,lon,lat]=GMTexportfield2(coef,fname,planetrad,alt,cmp,area,cutcirc,nanrng,res,onorout)
%
% Calculates a gradient field from a spherical harmonic potential field and 
% exports it as GMT readable grid.
%
% Use this instead of GMTexportfield for GMT newer than 5.4.2.
% The order of the latitude points is changed.  
%
% INPUT:
%
% coef          The spherical harmonic coefficients of the potential
% fname         Name for the output grid file
% planetrad     Radius of the planet
% alt           altitude for which you want to show the field
% cmpwrite      component: 
%               1 Radial (outward)
%               2 Southward (colatitudinal)
%               3 Eastward (longitudinal)
%               OR you can choose several: [1 2 3] (default: all)
% area          lon lat range for which to plot the field 
%               [lonmin latmax lonmax latmin]
% cutcirc       [TH clon ccola] parameters for the region you want to cut
%               out: TH is polar cap opening angle (onesided, like in 
%               Plattner & Simons (2014), ACHA), clon is longitude, ccola
%               is colatitude. ALL IN DEGREES
% nanrng        turn all values that are between -nanrng and nanrng into
%               nans (might help to make stuff invisible)
% res           Resolution in degrees per pixel
% onorout       Is the coef in ADDMON (0 or nothing), or ADDMOUT (1)?

%
% OUTPUT:
% 
% datB          Calculated data, set to zero outside of cutcirc, and maybe
%               stuff turned into nans if nanrng is on
% lon,lat       Longitudes/latitudes of the model evaluation
%
% Last modified by plattner-at-alumni.ethz.ch, 5/4/2015


%if ~ischar(coef)

defval('planetrad',3376) % Mars polar radius
defval('alt',0)
defval('cmpwrite',[1 2 3])
defval('area',[0 90 360 -90])
defval('cutcirc',[])
defval('nanrng',0)
defval('res',0.2)
defval('onorout',0)


Lmax=sqrt(length(coef))-1;

[demsz,delsz,mz,lmc,mzin,mzo]=addmon(Lmax);
coefBalt=vecupderivative(coef,planetrad+alt,planetrad,Lmax,0);
lmcosiBalt=coef2lmcosi(coefBalt,onorout);
%lmcosiBalt=[delsz demsz reshape(insert(coefBalt,0,mzin),2,length(demsz))']; 
[datB,lon,lat]=elm2xyz(lmcosiBalt,res,area);    


if ~isempty(cutcirc)
    
     cTH=cutcirc(1);
     clon=cutcirc(2);
     ccola=cutcirc(3);
        
    % To find the indices of the points we want: simply rotate all the
    % points to the North Pole and then keep the ones with colatitudes
    % smaller than the opening angle.
    [thetap,phip]=rottp((90-lat)*pi/180,lon*pi/180,clon*pi/180,ccola*pi/180,0);
    % Now find out which ones are close enough to the north pole
    index=find(thetap<=cTH*pi/180);
    % And set everything to zero that is outside of that circle    
    for cmp=1:3
        datBcut{cmp}=zeros(size(datB{cmp}));
        datBcut{cmp}(index)=datB{cmp}(index);
    end
    datB=datBcut;
    
end

if nanrng   
   for cmpind=1:3
      indeggs=find(abs(datB{cmpind})<=nanrng);
      datB{cmpind}(indeggs)=nan;
   end
end



for cmpind=1:length(cmpwrite)
    fprintf('Data range = %g to %g nT\n',...
        min(min(datB{cmpwrite(cmpind)})),max(max(datB{cmpwrite(cmpind)})))
    %grdwrite2p(lon,flipud(lat),flipud(datB{cmpwrite(cmpind)}),...
    %           sprintf('%s_cmp%d.nc',fname,cmpwrite(cmpind)));
    grdwrite2p(lon,lat,datB{cmpwrite(cmpind)},...
        sprintf('%s_cmp%d.nc',fname,cmpwrite(cmpind)));  
end


varns={datB,lon,lat};
varargout=varns(1:nargout);

