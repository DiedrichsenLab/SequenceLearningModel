function [T]=slm_genTrial(projName,trial,varargin)
% function [T]=slm_genTrial(projName,trial,varargin)
%
% Function that creates trial structure T for trial of the project specified by projName
% Usage: [T]=slm_genTrial('sr2_rt',1:10,1,1);

if ~exist('varargin','var')
    s=1; % subj num
    b=1; % block num
else
    s=varargin{1}; % subj num
    b=varargin{2}; % block num
end

switch projName
    
    case 'sr1' % simple sequence, 4 presses, free RT
        [T]=slm_genTrial_sr1(trial,s,b);
    
    case 'sr2_rt' % srtt, single press, forced RT
        [T]=slm_genTrial_sr2_rt(trial,s,b);
        
    case 'sr2_rs' % random sequence, 5 presses, forced RT
        [T]=slm_genTrial_sr2_rs(trial,s,b);
        
    case 'SeqEye' % longer sequence, 14 presses, free RT
        [T]=slm_genTrial_SeqEye(trial);
    
    otherwise
        error('Undefined project or genTrial function for this project')

end

end

function [T]=slm_genTrial_sr1(trial,s,b)

% set path to target folder
path2tgt='/Volumes/motorcontrol/robotcode/projects/SequenceRepetition/sr1/target';

% read in example target file from this project
tgt=dir(fullfile(path2tgt,sprintf('sr1_s%02d_b%02d.tgt',s,b))); % select tgt file
D=dload(fullfile(path2tgt,tgt(1).name)); % load the tgt file
if numel(trial)>numel(D.seqNum); trial=1:numel(D.seqNum); end
T=getrow(D,trial); % select how many trials of this block

% ensure "mandatory" fields are present (common fields across experiments)
T.TN=(1:numel(trial))'; % keep track of how many trials
T.stimulus=string(T.cueP); % visual cue (instruction to participant) 
for t=1:numel(trial); T.stimulus(t,1:4)=[str2double(T.stimulus{t}(1)),str2double(T.stimulus{t}(2)),str2double(T.stimulus{t}(3)),str2double(T.stimulus{t}(4))]; end
T.stimulus=str2double(T.stimulus);
T.forcedPressTime=(ones(numel(trial),1)*Inf); % forced-RT time in forced RT paradigm (with respect to beginning of trial: Inf means free RT)
T.stimTime=zeros(numel(trial),1); % when the stimulus appears (with respect to beginning of trial: 0 means right away)
T.numPress=ones(numel(trial),1)*numel(T.stimulus(1,:)); % how many finger presses are required in this task

end

function [T]=slm_genTrial_sr2_rt(trial,s,b)

% set path to target folder
path2tgt='/Volumes/motorcontrol/robotcode/projects/SequenceRepetition/sr2/target';

% read in example target file from this project
tgt=dir(fullfile(path2tgt,sprintf('sr2_rt_s%02d_b%02d.tgt',s,b))); % select tgt file
D=dload(fullfile(path2tgt,tgt(1).name)); % load the tgt file
if numel(trial)>numel(D.seqNum); trial=1:numel(D.seqNum); end
T=getrow(D,trial); % select how many trials of this block

% ensure "mandatory" fields are present (common fields across experiments)
T.TN=(1:numel(trial))'; % keep track of how many trials
T.stimulus=string(T.cueP); % visual cue (instruction to participant) 
for t=1:numel(trial); T.stimulus(t,1:5)=[str2double(T.stimulus{t}(1)),str2double(T.stimulus{t}(2)),str2double(T.stimulus{t}(3)),str2double(T.stimulus{t}(4)),str2double(T.stimulus{t}(5))]; end
T.stimulus=str2double(T.stimulus);
T.forcedPressTime=(ones(numel(trial),1)*3200); % forced-RT time in forced RT paradigm (with respect to beginning of trial: Inf means free RT)
T.stimTime=ones(numel(trial),numel(T.stimulus(1,:))).*(T(1).forcedPressTime-T.prepTime); % when the stimulus appears (with respect to beginning of trial: 0 means right away)
T.numPress=ones(numel(trial),1); % how many finger presses are required in this task

end

function [T]=slm_genTrial_sr2_rs(trial,s,b)

% set path to target folder
path2tgt='/Volumes/motorcontrol/robotcode/projects/SequenceRepetition/sr2/target';

% read in example target file from this project
tgt=dir(fullfile(path2tgt,sprintf('sr2_rs_s%02d_b%02d.tgt',s,b))); % select tgt file
D=dload(fullfile(path2tgt,tgt(1).name)); % load the tgt file
if numel(trial)>numel(D.seqNum); trial=1:numel(D.seqNum); end
T=getrow(D,trial); % select how many trials of this block

% ensure "mandatory" fields are present (common fields across experiments)
T.TN=(1:numel(trial))'; % keep track of how many trials
T.stimulus=string(T.cueP); % visual cue (instruction to participant) 
for t=1:numel(trial); T.stimulus(t,1:5)=[str2double(T.stimulus{t}(1)),str2double(T.stimulus{t}(2)),str2double(T.stimulus{t}(3)),str2double(T.stimulus{t}(4)),str2double(T.stimulus{t}(5))]; end
T.stimulus=str2double(T.stimulus);
T.forcedPressTime=(ones(numel(trial),1)*3200); % forced-RT time in forced RT paradigm (with respect to beginning of trial: Inf means free RT)
T.stimTime=ones(numel(trial),numel(T.stimulus(1,:))).*(T(1).forcedPressTime-T.prepTime); % when the stimulus appears (with respect to beginning of trial: 0 means right away)
T.numPress=ones(numel(trial),1)*numel(T.stimulus(1,:)); % how many finger presses are required in this task

end

function [T]=slm_genTrial_SeqEye(trial)

% set path to target folder
path2tgt='/Volumes/motorcontrol/robotcode/projects/SeqEye_Horizon/target';

% read in example target file from this project
tgt=dir(fullfile(path2tgt,'*1.tgt')); % select tgt file
D=dload(fullfile(path2tgt,tgt(1).name)); % load the tgt file
if numel(trial)>numel(D.seqNum); trial=1:numel(D.seqNum); end
T=getrow(D,trial); % select how many trials of this block

% ensure "mandatory" fields are present (common fields across experiments)
T.TN=(1:numel(trial))'; % keep track of how many trials
T.stimulus=string(T.cueP); % visual cue (instruction to participant) 
for t=1:numel(trial); T.stimulus(t,1:5)=[str2double(T.stimulus{t}(1)),str2double(T.stimulus{t}(2)),str2double(T.stimulus{t}(3)),str2double(T.stimulus{t}(4)),str2double(T.stimulus{t}(5))]; end
T.stimulus=str2double(T.stimulus);
T.forcedPressTime=(ones(numel(trial),1)*Inf); % forced-RT time in forced RT paradigm (with respect to beginning of trial: Inf means free RT)
T.stimTime=zeros(numel(trial),1); % when the stimulus appears (with respect to beginning of trial: 0 means right away)
T.numPress=ones(numel(trial),1)*numel(T.stimulus(1,:)); % how many finger presses are required in this task

end

