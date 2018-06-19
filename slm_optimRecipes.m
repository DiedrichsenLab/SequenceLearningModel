slm_NoiselessFitModel('Fit' , Dall , 'planFunc' , 'exp')
slm_NoiselessFitModel('Fit' , Dall , 'planFunc' , 'ramp')
slm_NoiselessFitModel('Fit' , Dall , 'planFunc' , 'logistic');

slm_NoiselessFitModel('Fit' , Dall , 'planFunc' , 'box_logistic')
slm_NoiselessFitModel('Fit' , Dall , 'planFunc' , 'box_exp')
slm_NoiselessFitModel('Fit' , Dall , 'planFunc' , 'box_ramp')



slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'exp')
slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'logistic');
slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'ramp')

slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'box_logistic')
slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'box_ramp')
slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'box_exp')


slm_NoiselessFitModel('Fit' , Dall , 'planFunc' , 'exp','NameExt' , 'dummy' , 'NumPresses' , 5)


X  = squeeze(SIM.X(:,:,3));
B  = squeeze(SIM.B(:,:,3));
figure('color' , 'white')
colorz = {'b' ,'r' ,'g' ,'m' ,'c'};
for p = 1:5
    subplot(5,1,p)
    plot(X(p,:) , 'color' , colorz{p} , 'LineWidth' , 3)
    hold on
    line()
end