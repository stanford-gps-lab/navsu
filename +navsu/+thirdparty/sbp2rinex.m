function status = sbp2rinex(filename, varargin)
%% Calls the Swift Binary Protocol to RINEX conversion exe
% Mandatory input:
%   filename      - SBP or SBP JSON filename
% Optional inputs as varargin:
%  'interval'     - desired output interval of observations [s]
%  'outdir'       - output directory for parsed file. nominally same as
%                   input directory
%  'span'         - time span [h] but seems to require a start epoch
%  'epochStart'   - GPS epoch of start time of output
%  'epochEnd'     - GPS epoch of end time of output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% The help info from SBP2RINEX exe %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts Swift Navigation receiver SBP binary and SBP JSON raw data logs to RINEX files.
%
%  Usage:
%
%      sbp2rinex [option ...] file
%  SBAS message file complies with RTKLIB SBAS/LEX message format.
%
%  Options [default]:
%
%      file         receiver log file
%      -ts y/m/d h:m:s  start time [all]
%      -te y/m/d h:m:s  end time [all]
%      -tr y/m/d h:m:s  approximated time for RTCM/CMR/CMR+ messages
%      -ti tint     observation data interval (s) [all]
%      -span span   time span (h) [all]
%      -r format    log format type:
%                   sbp       - Swift Navigation SBP
%                   json      - Swift Navigation SBP-JSON
%      -ro opt,opt  receiver option(s) (use comma between multiple opt):
%                   CONVBASE  - convert base station observations
%                   EPHALL    - include all ephemeris
%                   INVCP     - invert carrier phase polarity
%                   OBSALL    - include observations with RAIM flag set
%      -f freq      number of frequencies [2]
%      -hc comment  rinex header: comment line
%      -hm marker   rinex header: marker name
%      -hn markno   rinex header: marker number
%      -ht marktype rinex header: marker type
%      -ho observ   rinex header: oberver name and agency separated by /
%      -hr rec      rinex header: receiver number, type and version separated by /
%      -ha ant      rinex header: antenna number and type separated by /
%      -hp pos      rinex header: approx position x/y/z separated by /
%      -hd delta    rinex header: antenna delta h/e/n separated by /
%      -v ver       rinex version [3.00]
%      -oi          include iono correction in rinex nav header [off]
%      -ot          include time correction in rinex nav header [off]
%      -ol          include leap seconds in rinex nav header [off]
%      -halfc       half-cycle ambiguity correction [off]
%      -mask   [sig[,...]] signal mask(s) (sig={G|R|E|J|S|C|I}L{1C|1P|1W|...})
%      -nomask [sig[,...]] signal no mask (same as above)
%      -x sat       exclude satellite
%      -y sys       exclude systems (G:GPS,R:GLO,E:GAL,J:QZS,S:SBS,C:BDS,I:IRN)
%      -d dir       output directory [same as input file]
%      -c staid     use RINEX file name convention with staid [off]
%      -o ofile     output RINEX OBS file
%      -n nfile     output RINEX NAV file
%      -g gfile     output RINEX GNAV file
%      -h hfile     output RINEX HNAV file
%      -q qfile     output RINEX QNAV file
%      -l lfile     output RINEX LNAV file
%      -s sfile     output SBAS message file
%      -trace level output trace level [off]
%
%  If not any output file is specified, default output files <file>.obs,
%  <file>.nav, <file>.gnav, <file>.hnav, <file>.qnav, <file>.lnav and
%  <file>.sbs are used. Empty output files are deleted after processing.
%
%  If log format type is not specified, format type is recognized by the input
%  file extension as following:
%      *.sbp        Swift Navigation SBP
%      *.json       Swift Navigation SBP-JSON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addParameter('interval', []);
p.addParameter('outdir',[]);
p.addParameter('span',[]);
p.addParameter('epochStart',[]);
p.addParameter('epochEnd',[]);

parse(p, varargin{:});
res        = p.Results;
interval   = res.interval;
span       = res.span;
outdir     = res.outdir;
epochStart = res.epochStart;
epochEnd   = res.epochEnd;

% The exe is expected to be in the +thirdparty folder
s = what('utility');
pathExe = [s.path '\+thirdparty\sbp2rinex.exe'];

% Go through each of the options and add text to the system inputs
optionText = [];

if ~isempty(interval)
    optionText = [optionText ' -ti ' num2str(interval)];
end

if ~isempty(outdir)
    optionText = [optionText ' -d "' outdir '"'];
end

if ~isempty(epochStart)
   % convert from GPS epoch to calendar date and time;
    %-ts y/m/d h:m:s  start time [all]
   [yr,mn,dy,hr,min,sec] = navsu.time.epochs2cal(epochStart);
   optionTexti = [' -ts ' num2str(yr) '/' num2str(mn,'%02i') '/' num2str(dy,'%02i') ' ' ...
       num2str(hr,'%02i') ':' num2str(min,'%02i') ':' num2str(round(sec),'%02i')];
   optionText = [optionText optionTexti];
end

if ~isempty(epochEnd)
   % convert from GPS epoch to calendar date and time;
    %-te y/m/d h:m:s  start time [all]
   [yr,mn,dy,hr,min,sec] = navsu.time.epochs2cal(epochEnd);
   optionTexti = [' -te ' num2str(yr) '/' num2str(mn,'%02i') '/' num2str(dy,'%02i') ' ' ...
       num2str(hr,'%02i') ':' num2str(min,'%02i') ':' num2str(round(sec),'%02i')];
   optionText = [optionText optionTexti];
end

if ~isempty(span)
    optionText = [optionText ' -span ' num2str(span) ''];
end


filenameText = [' "' filename '"' ];

status = system(['"' pathExe '"' optionText filenameText]);


end







