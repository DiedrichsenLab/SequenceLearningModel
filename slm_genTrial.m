function [T]=slm_genTrial(projName,nTrials)
% function [T]=slm_genTrial(projName,nTrials)
%
% Function that creates trial structure T for nTrials of the project specified by projName
% Usage: [T]=slm_genTrial('sr2_rt',3);

switch projName
    
    case 'sr1' % simple sequence, 4 presses, free RT
        [T]=slm_genTrial_sr1(nTrials);
    
    case 'sr2_rt' % srtt, single press, forced RT
        [T]=slm_genTrial_sr2_rt(nTrials);
        
    case 'sr2_rs' % random sequence, 5 presses, forced RT
        [T]=slm_genTrial_sr2_rs(nTrials);
        
    case 'SeqEye' % longer sequence, 14 presses, free RT
        [T]=slm_genTrial_SeqEye(nTrials);
    
    otherwise
        error('Undefined project or genTrial function for this project')

end

end

function [T]=slm_genTrial_sr1(nTrials)

% set path to target folder
path2tgt='/Volumes/motorcontrol/robotcode/project/SequenceRepetition/sr1/target';

% read in example target file from this project
tgt=dir(fullfile(path2tgt,'sr1*1.tgt')); % first block/subject, just pick one
D=dload(fullfile(path2tgt,tgt(1).name)); % load the tgt file
T=getrow(D,1:nTrials); % select how many trials of this block

% ensure "mandatory" fields are present (common fields across experiments)
T.TN=(1:nTrials)'; % keep track of how many trials
T.stimulus=string(T.cueP(1:nTrials)); % visual cue (instruction to participant) 
T.stimTime=T.prepTime; % when the stimulus appears (with respect to beginning of trial: 0 means right away)
T.numPress=ones(nTrials,1)*size(T.stimulus{1},2); % how many finger presses are required in this task
T.forcedPressTime=[ones(nTrials,1)*Inf, ones(nTrials,1)*Inf]; % required response window in forced RT paradigm (with respect to beginning of trial: Inf means free RT)

end

function [T]=slm_genTrial_sr2_rt(nTrials)

% set path to target folder
path2tgt='/Volumes/motorcontrol/robotcode/project/SequenceRepetition/sr2/target';

% read in example target file from this project
tgt=dir(fullfile(path2tgt,'sr2_rt*1.tgt')); % first block/subject, just pick one
D=dload(fullfile(path2tgt,tgt(1).name)); % load the tgt file
T=getrow(D,1:nTrials); % select how many trials of this block

% ensure "mandatory" fields are present (common fields across experiments)
T.TN=(1:nTrials)'; % keep track of how many trials
T.stimulus=string(T.cueP(1:nTrials)); % visual cue (instruction to participant) 
T.stimTime=2400-T.prepTime; % when the stimulus appears (with respect to beginning of trial: 0 means right away)
T.numPress=ones(nTrials,1); % how many finger presses are required in this task
T.forcedPressTime=[(ones(nTrials,1)*2400)-100, ones(nTrials,1)*2400+100]; % % required response window in forced RT paradigm (with respect to beginning of trial: Inf means free RT)

end

function [T]=slm_genTrial_sr2_rs(nTrials)

% set path to target folder
path2tgt='/Volumes/motorcontrol/robotcode/project/SequenceRepetition/sr2/target';

% read in example target file from this project
tgt=dir(fullfile(path2tgt,'sr2_rs*1.tgt')); % first block/subject, just pick one
D=dload(fullfile(path2tgt,tgt(1).name)); % load the tgt file
T=getrow(D,1:nTrials); % select how many trials of this block

% ensure "mandatory" fields are present (common fields across experiments)
T.TN=(1:nTrials)'; % keep track of how many trials
T.stimulus=string(T.cueP(1:nTrials)); % visual cue (instruction to participant) 
T.stimTime=2400-T.prepTime; % when the stimulus appears (with respect to beginning of trial: 0 means right away)
T.numPress=ones(nTrials,1)*size(T.stimulus{1},2); % how many finger presses are required in this task
T.forcedPressTime=[(ones(nTrials,1)*2400)-100, ones(nTrials,1)*2400+100]; % required response window in forced RT paradigm (with respect to beginning of trial: Inf means free RT)

end

function [T]=slm_genTrial_SeqEye(nTrials)

% set path to target folder
path2tgt='/Volumes/motorcontrol/robotcode/project/SeqEye/target';

% read in example target file from this project
tgt=dir(fullfile(path2tgt,'*1.tgt')); % first block/subject, just pick one
D=dload(fullfile(path2tgt,tgt(1).name)); % load the tgt file
T=getrow(D,1:nTrials); % select how many trials of this block

% ensure "mandatory" fields are present (common fields across experiments)
T.TN=(1:nTrials)'; % keep track of how many trials
T.stimulus=string(T.cueP(1:nTrials)); % visual cue (instruction to participant) 
T.stimTime=zeros(nTrials,1); % when the stimulus appears (with respect to beginning of trial: 0 means right away)
T.numPress=ones(nTrials,1)*size(T.stimulus{1},2); % how many finger presses are required in this task
T.forcedPressTime=[ones(nTrials,1)*Inf, ones(nTrials,1)*Inf]; % required response window in forced RT paradigm (with respect to beginning of trial: Inf means free RT)

end

