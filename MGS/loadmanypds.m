function [posX,posY,posZ,magX,magY,magZ,sam,sap,sao,date,...
    RMSBvalsX,RMSBvalsY,RMSBvalsZ,BstatX,BstatY,BstatZ,...
    BdynX,BdynY,BdynZ,nerrors]=loadmanypds(whichones,directory,downsample)
% function [posX,posY,posZ,magX,magY,magZ,sam,sap,sao,date,...
%    RMSBvalsX,RMSBvalsY,RMSBvalsZ,BstatX,BstatY,BstatZ,...
%    BdynX,BdynY,BdynZ,nerrors]=loadmanypds(whichones,directory,downsample)
%
% Loads a series of .STS files and combines their data.
% 
% INPUT:
%
% whichones     Vector of .STS filename numbers
% directory     path to where the files can be found
% downsample    take every datum (=1) or every second (=2)... default is 1
%
% OUTPUT:
%
% posX etc.     The combined data. See loadpds
% nerrors       nerrors(1): number of not Mars files (they weren't read)
%               nerrors(2): number of not planetocentric files 
%                           (they weren't read)
% 
% Example: 
% whichones = [97257:97365 98001:98365 99001:99088];
% directory = '/yourfoldername/Mars/MGS_Data/PremappingDataRaw';
% [posX,posY,posZ,magX,magY,magZ,sam,sap,sao,date,...
%    RMSBvalsX,RMSBvalsY,RMSBvalsZ,BstatX,BstatY,BstatZ,...
%    BdynX,BdynY,BdynZ,nerrors]=loadmanypds(whichones,directory);
%
% Last modified by plattner-at-alumni.ethz.ch, 06/27/2016

defval('downsample',1);

nerrors=[0 0];

posX=[];
posY=[];
posZ=[];
magX=[];
magY=[];
magZ=[];
sam=[];
sap=[];
sao=[];
date=[];
RMSBvalsX=[];
RMSBvalsY=[];
RMSBvalsZ=[];
BstatX=[];
BstatY=[];
BstatZ=[];
BdynX=[];
BdynY=[];
BdynZ=[];


%h=waitbar(0,'Loading many days');
for i=1:length(whichones)    
    filename=fullfile(directory,sprintf('%05d.STS',whichones(i)));
    try
    [posXi,posYi,posZi,magXi,magYi,magZi,sami,sapi,saoi,datei,...
        RMSBvalsXi,RMSBvalsYi,RMSBvalsZi,BstatXi,BstatYi,BstatZi,...
        BdynXi,BdynYi,BdynZi,alerttype]...
        =loadpds(filename);    
    posX=[posX;posXi(1:downsample:end)];
    posY=[posY;posYi(1:downsample:end)];
    posZ=[posZ;posZi(1:downsample:end)];
    magX=[magX;magXi(1:downsample:end)];
    magY=[magY;magYi(1:downsample:end)];
    magZ=[magZ;magZi(1:downsample:end)];
    sam=[sam;sami(1:downsample:end)];
    sap=[sap;sapi(1:downsample:end)];
    sao=[sao;saoi(1:downsample:end)];
    date=[date;datei(1:downsample:end,:)];
    RMSBvalsX=[RMSBvalsX;RMSBvalsXi(1:downsample:end)];
    RMSBvalsY=[RMSBvalsY;RMSBvalsYi(1:downsample:end)];
    RMSBvalsZ=[RMSBvalsZ;RMSBvalsZi(1:downsample:end)];
    BstatX=[BstatX;BstatXi(1:downsample:end)];
    BstatY=[BstatY;BstatYi(1:downsample:end)];
    BstatZ=[BstatZ;BstatZi(1:downsample:end)];
    BdynX=[BdynX;BdynXi(1:downsample:end)];
    BdynY=[BdynY;BdynYi(1:downsample:end)];
    BdynZ=[BdynZ;BdynZi(1:downsample:end)];
    
    if alerttype==1
        nerrors(1)=nerrors(1)+1;
    elseif alerttype==2
        nerrors(2)=nerrors(2)+1;
    end
        
    catch
        fprintf('day %d does not exist\n',whichones(i));
    end
    %waitbar(i/length(whichones),h);
end
%delete(h)