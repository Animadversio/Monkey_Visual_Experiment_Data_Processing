% Fast Sorting of today's experiments
function batchProcessPL2(ephysFNs, varargin)
if nargin == 0 % get today's animal
    animal = input("Animal:",'s');
    date = input("Date:",'s');
    if isempty(date)
        date = datetime();
    end
    searchstr = string([animal,'-',datestr(date,"ddmmyyyy"),'*.pl2']);
    ephysFNs = string(ls("N:\Data-Ephys-Raw\"+searchstr));
end
if nargin <= 1, sdf = 'sdf'; end 
preMeta = arrayfun(@(fn)struct('ephysFN',fn{1}(1:end-4),'sdf',sdf), ephysFNs);
Project_General_copyMissingFiles(preMeta)
for iFn = 1:numel(ephysFNs)
ephysFN = ephysFNS(iFn);
sortPL2(ephysFN, varargin);
end
end

function sortPL2(ephysFN, varargin)
p = inputParser;
defaultEquipment =  'OMNIPLEX';
defaultExpControl = 'ML';
% defaultBaseline =   90;
% defaultRasterWindow =   [0 200];
% defaultExpControlName = '170619_ringo_screening.bhv2';

addRequired(p,  'ephysFN')
% addOptional(p,  'expControlFN',defaultExpControlName)
addParameter(p, 'expControl', defaultExpControl)
addParameter(p, 'equipment', defaultEquipment)
% addParameter(p, 'baselineWindowLength', defaultBaseline )
% addParameter(p, 'rasterWindow', defaultRasterWindow )
addParameter(p, 'sdf', 'sdf') % 'raster'

parse(p,ephysFN, varargin{:})

meta = p.Results;
meta = myPaths(meta);

if strcmp(meta.sdf, 'sdf')
bigMatPath = fullfile(meta.pathMat,[meta.ephysFN '.mat']);
elseif strcmp(meta.sdf, 'raster')
bigMatPath = fullfile(meta.pathMat,[meta.ephysFN '_spike.mat']);
else
error;
end
bigMatrixExists = exist(bigMatPath,'file');

[spikeChans,lfpChans,timeline,spikeID] = plxread_fullExperiment_vcrp(meta);
% Save the source of the channel per isolated unit/hash
if ~spikeID(1) 
    meta.spikeID = spikeID + 1;
else
    meta.spikeID = spikeID ;
end
% read words/bits from Plexon file
Trials = plxread_loadWordsInPlexonFile_v2(meta) ;
if strcmp(meta.sdf, 'sdf')
savefast(bigMatPath,'meta','Trials','spikeChans','lfpChans','timeline');
elseif strcmp(meta.sdf, 'raster')
savefast(bigMatPath,'meta','Trials','spikeChans','lfpChans','timeline');
else
error();
end
end