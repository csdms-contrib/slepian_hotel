function [binmean,binsig,longrid,colagrid,npoints,rmean,rsig]=dataBinning(data,lon,cola,r,lores,colares)
% [binmean,binsig,longrid,latgrid,npoints,rmean,rsig]=dataBinning(data,lon,lat,r,res)
%
% Bins the data and calculates the mean and std deviation for each bin. The
% bins are defined by their centerpoints, returned as longrid,latgrid
%
% INPUT:
%
% data      the data to be binned
% lon,cola  longitudes, colatitudes of the data points, both in degrees,
%           0<=lon<=360, 0<=cola<=180
% r         radial position of data points
% lores     Longitudinal Step size in binning grid 
%           (e.g. a cell every lores=0.5 degrees)
% colares   Latitudinal Step size in binning grid 
%
% OUTPUT:
%
% binmean   binmean{1}: radial mean, 
%           binmean{2}: colatitudinal mean,
%           binmean{3}: longitudinal mean
% binsig    same as binmean but for standard deviation
% longrid,
% colagrid  centerpoints of binning grid
% npoints   number of points in each binning cell 
% rmean     mean radial location in each binning cell
% rsig      std deviation in each binning cell
%
% Last modified by plattner-at-alumni.ethz.ch, 5/26/2015


% We need to make the grid. In the future I might want to allow
% variable longitudinal resolution depending on latitude. For now just
% simple
longrid=round(min(lon)):lores:round(max(lon));
colagrid=round(min(cola)):colares:round(max(cola));

% Initialize all the stuff
for cmp=1:3
   binmean{cmp}=zeros(length(longrid),length(colagrid)); 
   binsig{cmp}=zeros(length(longrid),length(colagrid)); 
end
npoints=zeros(length(longrid),length(colagrid));
rmean=zeros(length(longrid),length(colagrid));
rsig=zeros(length(longrid),length(colagrid));

% Now we fill up the matrix of vectors bindices by iterating through all of
% the data points and placing them in the right bin (appending their
% indices to the right vectors in bindices)

% Turn every lon/lat position into a coordinate
lonshift=min(longrid/lores)-1;
colashift=min(colagrid/colares)-1;
indlon=round(lon/lores)-lonshift;
indcola=round(cola/colares)-colashift;

for icola=1:length(colagrid)
    indicola=find(indcola==icola);
    % Reduce all the data sets to this index set
    for cmp=1:3
        datcola{cmp}=data{cmp}(indicola);
    end
    rcola=r(indicola);
    indloncola=indlon(indicola);
    % Now, tease out for each cell along this latitude line
    for ilo=1:length(longrid)
        indis=find(indloncola==ilo);                        
        for cmp=1:3  
            binmean{cmp}(ilo,icola)=mean(datcola{cmp}(indis));
            binsig{cmp}(ilo,icola)=stdp(datcola{cmp}(indis));
        end
        npoints(ilo,icola)=length(indis);
        rmean(ilo,icola)=mean(rcola(indis));
        rsig(ilo,icola)=stdp(rcola(indis));
    end
end

% Need to flip them all
npoints=npoints';
rmean=rmean';
rsig=rsig';
for cmp=1:3
    binmean{cmp}=binmean{cmp}';
    binsig{cmp}=binsig{cmp}';
end


end

function y=stdp(v)
    % This function because Octave returns [] for v=[] but Matlab returns NaN
    if isempty(v)
        y=NaN;
    else
        y=std(v);
    end
end


