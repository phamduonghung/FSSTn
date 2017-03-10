% Copyright (c) Haizhao Yang, Stanford University, 2014

global SetPath
global CSPT
global MSPT
		
type = computer;

if strcmp(type,'MAC2'),
  CSPT = ':';
  SetPath = [pwd, CSPT];
  MSPT = ';';
elseif isunix,
  % Mac OS X returns isunix=1
  CSPT = '/';
  SetPath = [pwd, CSPT];
  MSPT = ':';
elseif strcmp(type(1:2),'PC');
  CSPT = '\';	  
  SetPath = [pwd, CSPT];  
  MSPT = ';';
end

disp('Begin to set MATLAB path...')

back = CSPT;
tempPath = path;
front = [MSPT SetPath];
tempPath = [tempPath front];

tempPath = [tempPath front 'results' back];
tempPath = [tempPath front 'comparison' back];

% Hau-Tieng's code
tempPath = [tempPath front 'sst' back];

tempPath = [tempPath front 'sst' back 'tool' back];

tempPath = [tempPath front 'sst' back 'tool' back 'blending' back];

% Multitaper
tempPath = [tempPath front 'tftb' back];

tempPath = [tempPath front 'tftb' back 'data' back];

tempPath = [tempPath front 'tftb' back 'demos' back];

tempPath = [tempPath front 'tftb' back 'mfiles' back];

tempPath = [tempPath front 'tftb' back 'tests' back];

tempPath = [tempPath front 'tftb' back 'multiTaper' back];

% ConceFT code
tempPath = [tempPath front 'conceFT' back];

tempPath = [tempPath front 'conceFT' back 'Conceft' back];

tempPath = [tempPath front 'conceFT' back 'Conceft' back 'Morse' back];

% EMD
tempPath = [tempPath front 'FastEMD' back];

% SynLab
tempPath = [tempPath front 'SynLab' back];

tempPath = [tempPath front 'SynLab' back 'results' back];

tempPath = [tempPath front 'SynLab' back 'Source' back];

tempPath = [tempPath front 'SynLab' back 'Source' back 'SS_CT_2D' back];
tempPath = [tempPath front 'SynLab' back 'Source' back 'SS_CT_2D' back 'demo' back];
tempPath = [tempPath front 'SynLab' back 'Source' back 'SS_CT_2D' back 'src' back];

tempPath = [tempPath front 'SynLab' back 'Source' back 'SS_WP_1D' back];
tempPath = [tempPath front 'SynLab' back 'Source' back 'SS_WP_1D' back 'demo' back];
tempPath = [tempPath front 'SynLab' back 'Source' back 'SS_WP_1D' back 'src' back];


path(tempPath);

disp('Begin to compile MEX files...');
rootDir = pwd;
cd(['SynLab' CSPT 'Source' CSPT 'SS_CT_2D' CSPT 'src' CSPT]);
mex SS_polar.c;
mex SS_polar_old.c;
mex SS_polar_v2.c;
mex SS_polar_v1.c;
cd(rootDir);

disp('Path set!');

clear tempPath front back
clear SetPath MATLABVERSION CSPT
clear type MSPT

