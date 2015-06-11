%% Configuration File
% What you need to use these scripts:
% Filepath to shape file
% Filepath template to tiff files
% Filter for region within shape file
% Range of years for which you have tiff files

cd(fileparts(mfilename('fullpath')))

%% Entering the directories of input
A=uigetdir('.','Enter the shape file directory');
shapepath = strcat(A,'/nalcmsmx05gw.shp');

A=uigetdir('.','Enter the temperature directory');
temppath=strcat(A,'/tmax_avg_%d_%02d');

A=uigetdir('.','Enter the precipitation directory');
prcppath=strcat(A,'/prcp_ttl_%d_%02d');

A=uigetdir('.','Enter the solar radiation directory');
sradpath = strcat(A,'/srad_avg_%d_%02d');

filter = 'Matorral tropical o subtropical';

years = 1980:2013;  
months = { 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December' };
                    
lowpctl = 5;
highpctl = 95;

time = @(y, m) (y - years(1))*12 + m; %To index into array

totalTime = 12*length(years);