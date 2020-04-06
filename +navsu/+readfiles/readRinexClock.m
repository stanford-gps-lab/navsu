function [epochs, clk, clkSig] = readRinexClock(filename, nMaxPrn,const)
% readRinexClock
% DESCRIPTION:
%   Parses a RINEX precise clock file, typically .clk files. 
%   Formerly called Read_GLO_CLK.
% INPUT:
%   filename    - path and name of the file to be parsed. 
% OPTIONAL INPUTS:
%   nMaxPrn     - maximum PRN number.  Use to set the size of the output
%                 matrices.
%   const       - single character indicating desired constellation to
%                 parse: 'G','R','E','C','S'
%   
% OUTPUT:
%   epochs      - vector of GPS epochs of the output clock biases
%   clk         - [nPrn x nEpochs] matrix of clock biases in seconds
%   clkSig      - [nPrn x nEpochs] matrix of clock bias sigmas in seconds
%
% See also: navsu.ftp.download, navsu.svOrbitClock

if nargin < 2
    nMaxPrn = 0;
end

if nargin < 3
    const = 'R';
end

temp = importdata(filename);
textData = temp.textdata;

% Pull PRN List
inds = [strfind(textData(:,1),'PRN LIST')];
inds = find(~cellfun(@isempty,inds));
prns = [];
for i = 1:length(inds)
    linetext = textData{inds(i)};
    gInds = strfind(linetext,const);
    prns = [prns; arrayfun(@(x) str2double(linetext(x+(1:2))), gInds)'];
end
prns = prns(~isnan(prns));

% Ensure there are PRNs available
if isempty(prns)
    fprintf(2, 'Error no GLONASS PRNs identified in file: %d\n', filename);
    epochs = [];
    clk = [];
    clkSig = [];
    return
end

% find time interval
tempEpochs = sort(unique(navsu.time.cal2epochs(temp.data(~isnan(temp.data(:,1)),1:6))));
dEpochs = mode(diff(tempEpochs));

% Initialize output variables
epochs = (0:dEpochs:(86400-dEpochs))';
clk = NaN(max(prns),length(epochs));
if nMaxPrn > 0
    clk = NaN(nMaxPrn,length(epochs));
end
clkSig = clk;

indAS = find(strcmp(textData(:,1),'AS'));
indOffset = size(temp.textdata,1)-size(temp.data,1);
for i = 1:length(prns)
    indG = intersect(find(strcmp(textData(:,2),[const num2str(prns(i), '%02d')])),indAS)-indOffset;
    
    year = temp.data(indG,1);
    month = temp.data(indG,2);
    day = temp.data(indG,3);
    hours = temp.data(indG,4);
    minutes = temp.data(indG,5);
    seconds = temp.data(indG,6);
    
    jdx = round(hours*(3600/dEpochs) + minutes*(60/dEpochs) + seconds/dEpochs) + 1;
    indG(isnan(jdx)) = [];
    jdx(isnan(jdx)) = [];
    clk(prns(i),jdx) = temp.data(indG,8);
    k = find(mod(minutes,5) == 0 & mod(seconds, 60) ==0);
    
%     clk_sig(prns(i),jdx(k)) = temp.data(indG(k),9);
end


end
