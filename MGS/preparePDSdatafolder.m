function preparePDSdatafolder(parts,years,fromdir,todir)
% preparePDSdatafolder(parts,years,fromdir,todir)
%
% Copies the .STS data files into a single directory
%
% INPUT:
%
% parts     names of the folders (one struct entry per year)
% years     vector containing the years to be loaded
% fromdir   directory where the downloaded data files are
% todir     directory to which you want to store the files
%
% Last modified by plattner-at-alumni.ethz.ch, 6/27/2016
% 
% Example:
%
% years=1997:1999;
% parts{1}={'244_273SEP','274_304OCT','305_334NOV','335_365DEC'};
% parts{2}={'001_031JAN','032_059FEB','060_090MAR','091_120APR','121_151MAY',...
%    '152_181JUN','182_212JUL','213_243AUG','244_273SEP','274_304OCT','305_334NOV','335_365DEC'};
% parts{3}={'001_031JAN','032_059FEB','060_090MAR'};
% fromdir='/yourfolder/Mars/MGS_Data/MGS-M-MAG-3-PREMAP_FULLWORD-RES-MAG-V1.0/DATA/MARS/';
% todir='/yourfolder/Mars/MGS_Data/PremappingData/';
% preparePDSdatafolder(parts,years,fromdir,todir)



for i=1:length(years)
    for j=1:length(parts{i})
        fd=sprintf('%s%d/%s/MAG_PCENTRIC/*.STS',fromdir,years(i),parts{i}{j});
        system(sprintf('cp %s %s',fd,todir));
    end
end
