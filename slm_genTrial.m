function [T]=slm_genTrial(projName,trial,varargin)
% function [T]=slm_genTrial(projName,trial,varargin)
%
% Function that creates trial structure T for trial of the project specified by projName
% Usage: [T]=slm_genTrial('sr2_rt',1:10,1,1);
 
if isempty(varargin)
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
        [T]=slm_genTrial_SeqEye(trial,s,b);
    
    otherwise
        error('Undefined project or genTrial function for this project')
 
end
 
end
 
function [T]=slm_genTrial_sr1(trial,s,b)
 
% set path to target folder
path2tgt='/Volumes/motorcontrol/robotcode/projects/SequenceRepetition/sr1/target';
 
% read in example target file from this project
tgt=dir(fullfile(path2tgt,sprintf('sr1_s%02d_b%02d.tgt',s,b))); % select tgt file
if isempty(tgt); error('Target file not found for this combination of path, subject, and block number! Make sure that the inputs are correct.'); end
D=dload(fullfile(path2tgt,tgt(1).name)); % load the tgt file
if numel(trial)>numel(D.seqNum); trial=1:numel(D.seqNum); end
T=getrow(D,trial); % select how many trials of this block
 
% ensure "mandatory" fields are present (common fields across experiments)
T.TN=(1:numel(trial))'; % keep track of how many trials
T.stimulus=double(num2str(T.cueP))-double('0'); % visual cue (instruction to participant)
T.forcedPressTime=nan(numel(trial),1); % forced-RT time in forced RT paradigm (with respect to beginning of trial: Inf means free RT)
T.stimTime=zeros(numel(trial),1)*numel(T.stimulus(1,:)); % when the stimulus appears (with respect to beginning of trial: 0 means right away)
T.numPress=ones(numel(trial),1)*numel(T.stimulus(1,:)); % how many finger presses are required in this task
T.Horizon=nan(size(T.stimulus)); % % Horizon feature added. Default to full horizon
T.pressTime = nan(size(T.Horizon)); % to be filled with press times

end
 
function [T]=slm_genTrial_sr2_rt(trial,s,b)
 
% set path to target folder
path2tgt='/Volumes/motorcontrol/robotcode/projects/SequenceRepetition/sr2/target';
 
% read in example target file from this project
tgt=dir(fullfile(path2tgt,sprintf('sr2_rt_s%02d_b%02d.tgt',s,b))); % select tgt file
if isempty(tgt); error('Target file not found for this combination of path, subject, and block number! Make sure that the inputs are correct.'); end
D=dload(fullfile(path2tgt,tgt(1).name)); % load the tgt file
if numel(trial)>numel(D.seqNum); trial=1:numel(D.seqNum); end
T=getrow(D,trial); % select how many trials of this block
lastToneTime=2400;
 
% ensure "mandatory" fields are present (common fields across experiments)
T.TN=(1:numel(trial))'; % keep track of how many trials
T.stimulus=double(num2str(T.cueP))-double('0'); % visual cue (instruction to participant) 
T.forcedPressTime=(ones(numel(trial),1)*lastToneTime); % forced-RT time in forced RT paradigm (with respect to beginning of trial: Inf means free RT)
T.stimTime=repmat(T.forcedPressTime-T.prepTime,1,size(T.stimulus,2)); % when the stimulus appears (with respect to beginning of trial: 0 means right away)
T.numPress=ones(numel(trial),1); % how many finger presses are required in this task
T.Horizon=nan(size(T.stimulus)); % % Horizon feature added. Default to full horizon
T.pressTime = nan(size(T.Horizon)); % to be filled with press times

end
 
function [T]=slm_genTrial_sr2_rs(trial,s,b)
 
% set path to target folder
path2tgt='/Volumes/motorcontrol/robotcode/projects/SequenceRepetition/sr2/target';
 
% read in example target file from this project
tgt=dir(fullfile(path2tgt,sprintf('sr2_rs_s%02d_b%02d.tgt',s,b))); % select tgt file
if isempty(tgt); error('Target file not found for this combination of path, subject, and block number! Make sure that the inputs are correct.'); end
D=dload(fullfile(path2tgt,tgt(1).name)); % load the tgt file
if numel(trial)>numel(D.seqNum); trial=1:numel(D.seqNum); end
T=getrow(D,trial); % select how many trials of this block
T.lastToneTime=2400;
 
% ensure "mandatory" fields are present (common fields across experiments)
T.TN=(1:numel(trial))'; % keep track of how many trials
T.stimulus=double(num2str(T.cueP))-double('0'); % visual cue (instruction to participant)
T.forcedPressTime=(ones(numel(trial),1)*lastToneTime); % forced-RT time in forced RT paradigm (with respect to beginning of trial: Inf means free RT)
T.stimTime=repmat(T.forcedPressTime-T.prepTime,1,size(T.stimulus,2)); % when the stimulus appears (with respect to beginning of trial: 0 means right away)
T.numPress=ones(numel(trial),1)*numel(T.stimulus(1,:)); % how many finger presses are required in this task
T.Horizon=nan(size(T.stimulus)); % % Horizon feature added. Default to full horizon
T.pressTime = nan(size(T.Horizon)); % to be filled with press times

end
 
function [T]=slm_genTrial_SeqEye(trial,s,b)
 
% set path to target folder
path2tgt='/Volumes/motorcontrol/robotcode/projects/SeqEye_Horizon/target';

% set subj name and block type
subjID={'AJ1','AS1','CB1','CS1','DW1',... % subject identifier: the order matters!
    }; % Remaining subjects to be added
blockType='CLAT'; % block type options: 'CLAT'|'CT'|'IM'|'RAND' Right now this has to be specified here manually! Consider looping through in slm_testModel

% read in example target file from this project
tgt=dir(fullfile(path2tgt,sprintf('%s_tgtFiles/%s_%s_B%d.tgt',subjID{s},subjID{s},blockType,b))); % select tgt file
if isempty(tgt); error('Target file not found for this combination of path, subject, and block number! Make sure that the inputs are correct.'); end
D=dload(fullfile(tgt(1).folder,tgt(1).name)); % load the tgt file
D.seqNum=D.seqNumb; D.cueP=num2str(D.cueP);
if numel(trial)>numel(D.seqNum); trial=1:numel(D.seqNum); end
T=getrow(D,trial); % select how many trials of this block
 
% ensure "mandatory" fields are present (common fields across experiments)
T.TN=(1:numel(trial))'; % keep track of how many trials
T.stimulus=double(num2str(T.cueP))-double('0'); % visual cue (instruction to participant)
T.forcedPressTime=nan(numel(trial),1); % forced-RT time in forced RT paradigm (with respect to beginning of trial: Inf means free RT)
T.stimTime=zeros(numel(trial),numel(T.stimulus(1,:))); % when the stimulus appears (with respect to beginning of trial: 0 means right away)
T.numPress=ones(numel(trial),1)*numel(T.stimulus(1,:)); % how many finger presses are required in this task
for t=trial
    hrzn = T.Horizon(t); % Horizon feature added. stimTime will be the actual time that the stimulus came on.
    T.Horizon(t,1:numel(T.stimulus(1,:))) = ones(1,numel(T.stimulus(1,:)))*hrzn; T.Horizon(t,1:hrzn) = NaN;
end
T.stimTime(~isnan(T.Horizon)) = NaN; % set the stimTime for presses that are not shown at time 0 to NaN (depending on horizon size). These will be filled with press times
T.pressTime = nan(size(T.Horizon)); % to be filled with press times

end

