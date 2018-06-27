slm_NoiselessFitModel('FitMTRT' , Dall , 'planFunc' , 'exp')
slm_NoiselessFitModel('FitMTRT' , Dall , 'planFunc' , 'ramp')
slm_NoiselessFitModel('FitMTRT' , Dall , 'planFunc' , 'logistic');
slm_NoiselessFitModel('FitIPIRT' , Dall , 'planFunc' , 'arbitrary', 'Horizon' , [2],'initPlan' , 'uniform' , 'NameExt' , 'H2');
slm_NoiselessFitModel('FitIPIRT' , Dall , 'planFunc' , 'arbitrary', 'Horizon' , [2],'initPlan' , 'logistic' ,'NameExt' , 'H2logistic');
slm_NoiselessFitModel('FitIPIRT' , Dall , 'planFunc' , 'arbitrary', 'Horizon' , [2],'initPlan' , 'box_exp' ,'NameExt' , 'H2boxexp');
slm_NoiselessFitModel('FitIPIRT' , Dall , 'planFunc' , 'arbitrary', 'Horizon' , [2],'initPlan' , 'exp' ,'NameExt' , 'H2exp');



for box = 2:9
    NameExt = ['_Box' , num2str(box)];
%     slm_NoiselessFitModel('Fit' , Dall , 'planFunc' , 'box_logistic' , 'MsetField' , {'Box' , box}, 'NameExt' , NameExt);
    slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'box_logistic' , 'MsetField' , {'Box' , box}, 'NameExt' , NameExt);
end

for box = 2:9
    NameExt = ['_Box' , num2str(box)];
%     slm_NoiselessFitModel('Fit' , Dall , 'planFunc' , 'box_exp' , 'MsetField' , {'Box' , box}, 'NameExt' , NameExt);
    slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'box_exp' , 'MsetField' , {'Box' , box}, 'NameExt' , NameExt);
end

for box = 2:9
    NameExt = ['_Box' , num2str(box)];
%     slm_NoiselessFitModel('Fit' , Dall , 'planFunc' , 'box_ramp' , 'MsetField' , {'Box' , box}, 'NameExt' , NameExt);
    slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'box_ramp' , 'MsetField' , {'Box' , box}, 'NameExt' , NameExt);
end


%% 



slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'exp')
slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'logistic');
slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'ramp')


slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'box_logistic')
slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'box_ramp')
slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'box_exp');% , 'NumPresses' , 6)


slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'arbitrary' , 'Horizon' , 2, 'NameExt' , 'H2');
slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'arbitrary' , 'Horizon' , 2, 'NameExt' , 'H2logistic');
slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'arbitrary' , 'Horizon' , 2, 'NameExt' , 'H2boxexp');
slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'arbitrary' , 'Horizon' , 2, 'NameExt' , 'H2exp');




