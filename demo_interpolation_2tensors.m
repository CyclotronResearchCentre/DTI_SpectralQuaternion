%% DEMO 1: INTERPOLATION BETWEEN TWO OR FOUR TENSORS

% In this demo, we show how to interpolate between two tensors, using SQ or
% Log-Euclidean averaging. Some features of the tensors are shown, and
% different rescaling are tested.

% ________________________________________________________
% Copyright (C) 2013 University of Liege, Belgium
    
% Written by A. Collard & C. Phillips, 2013.
% Dept of Electrical Engineering and Computer Science &
% Cyclotron Research Centre, University of Liege, Belgium
% ________________________________________________________

clear all
close all
clc

%% Interpolations of two tensors

% Construction of the two 'extremities' tensors
alpha= (1/360)*2*pi;

R= [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
S1= difftensor(R*diag([10 1 1])*R');

alpha= (63/360)*2*pi;

R= [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
S2= difftensor(R*diag([40 4 1])*R');

% t : interpolation parameter

t= 0:0.125:1;

%% First case : with no rescaling of the orientation interpolation for SQ.

% Comparison of the Log-Euclidean and Spectral-Quaternion interpolation

first_comment = ['The first example will perform an interpolation between 2 tensors,'...
'using the Log-Euclidean method and the Spectral-Quaternion method'];
disp(first_comment);
display('with no rescaling of the orientation interpolation.');

display('For each curve, the evolution of the principal features of');
display('the tensors are also computed and illustrated.');

display('See the file to understand the used commands.');



% Slog : will contain the Log-Euclidean interpolation
S_LogE= difftensor;

% Sspe : will contain the Spectral quaternion interpolation
S_SQ= difftensor;

for ii = 1 : length(t)
    
    S_LogE(ii) = wmean2dt(S1,S2,t(ii),'LogE');
    S_SQ(ii) = wmean2dt(S1,S2,t(ii),'SQ','no');
end

% Illustration of interpolations

graph_display(S_LogE,0.3);
title('Log-Euclidean interpolation');
axis off

graph_display(S_SQ,0.3);
title('Spectral-Quaternion interpolation - no rescaling');
axis off

% Evolution of the determinant

det_loge= getDet(S_LogE);
det_spe= getDet(S_SQ);

% Determine the limits of the y-axis
minDet = min(min(det_loge),min(det_spe));
maxDet = max(max(det_loge),max(det_spe));

% Evolution of the Fractional Anisotropy

FA_loge= getFA(S_LogE);
FA_spe= getFA(S_SQ);

% Determine the limits of the y-axis
minFA = min(min(FA_loge),min(FA_spe));
maxFA = max(max(FA_loge),max(FA_spe));


% Evolution of the Hilbert Anisotropy

HA_loge= getHA(S_LogE);
HA_spe= getHA(S_SQ);

% Determine the limits of the y-axis

minHA = min(min(HA_loge),min(HA_spe));
maxHA = max(max(HA_loge),max(HA_spe));

% Evolution of the angular difference between first eigenvectors

phi_loge= zeros(1,length(t));

pv1= S1.EigVectors(:,1);

for ii = 1 : length(t)
    pvi= S_LogE(ii).EigVectors(:,1);
    phi_loge(ii)= acos(abs(pv1'*pvi));
end

phi_spe= zeros(1,length(t));

pv1= S1.EigVectors(:,1);
for ii = 1 : length(t)
    pvi= S_SQ(ii).EigVectors(:,1);
    phi_spe(ii)= mod(acos(abs(pv1'*pvi)), 2*pi);
end

% Determine the limits of the y-axis

minPhi = min(min(phi_loge),min(phi_spe));
maxPhi = max(max(phi_loge),max(phi_spe));

% Summarize those results in a figure
figure
% Ensure a clear presentation of the figure with subplots
position = [ 500 500 1200 300 ];
h = figure(3);
set(h, 'Position', position);

subplot(141)
plot(t,det_loge,'-ob','linewidth',2);
hold on
plot(t,det_spe,'-or','linewidth',2);
hold off
title('Det')

subplot(142)
plot(t,FA_loge,'-ob','linewidth',2);
hold on
plot(t,FA_spe,'-or','linewidth',2);
hold off
axis([0 1 minFA maxFA]);
title('FA');

subplot(143)
plot(t,HA_loge,'-ob','linewidth',2);
hold on
plot(t,HA_spe,'-or','linewidth',2);
hold off
axis([0 1 minHA maxHA]);
title('HA');

subplot(144)
plot(t,phi_loge,'-ob','linewidth',2);
hold on
plot(t,phi_spe,'-or','linewidth',2);
hold off
axis([0 1 minPhi maxPhi]);
title('\phi');
legend('Log-E','SQ');

display('Execution paused. Press any key to continue.');
pause;


%% Second case : with 'kappa' rescaling 
display('In the second example, a rescaling of the orientation interpolation');
display('is used. The evolutions of the other features (anisotropy, determinant,...)');
display('are not impacted.')

% Those weights are proposed in 
% A. Collard, S.Bonnabel, C. Phillips and R. Sepulchre, 'Anisotropy
% preserving DTI processing', 2013;
% in Equation (40).


% Sspe : will contain the Spectral quaternion interpolation
S_SQ= difftensor;

for ii = 1 : length(t)
    S_SQ(ii) = wmean2dt(S1,S2,t(ii),'SQ','kappa');
end

% Illustration of interpolations

graph_display(S_SQ,0.3);
title('Spectral-Quaternion interpolation - \kappa rescaling');
axis off

% Evolution of det, HA, FA are identical than with no rescaling.

% Evolution of the angular difference between first eigenvectors

phi_spe= zeros(1,length(t));

pv1= S1.EigVectors(:,1);
for ii = 1 : length(t)
    pvi= S_SQ(ii).EigVectors(:,1);
    phi_spe(ii)= mod(acos(abs(pv1'*pvi)), 2*pi);
end

% Determine the limits of the y-axis

minPhi = min(min(phi_loge),min(phi_spe));
maxPhi = max(max(phi_loge),max(phi_spe));

figure
plot(t,phi_loge,'-ob','linewidth',2);
axis([0 1 minPhi maxPhi]);
hold on
plot(t,phi_spe,'-or','linewidth',2);
legend('LogE','SQ with \kappa');
title('Orientation evolution - Angular difference');

display('Execution paused. Press any key to continue.')
pause
%% Third case : with 'anisotropic' rescaling 
display('The third example uses another rescaling of the orientation.');
display('In this case, the weights are multiplied by the ratio between the');
display('anisotropy of the considered tensor and the anisotropy of the average');
display('tensor.');
% In this case, the weights are multiplied by the ratio between the
% anisotropy of the considered tensor and the anisotropy of the average
% tensor.


% Sspe : will contain the Spectral quaternion interpolation
S_SQ= difftensor;

for ii = 1 : length(t)

    S_SQ(ii) = wmean2dt(S1,S2,t(ii),'SQ','HA');
end

% Illustration of interpolations

graph_display(S_SQ,0.3);
title('Spectral-Quaternion interpolation - HA rescaling');
axis off

% Evolution of det, HA, FA are identical than with no rescaling.

% Evolution of the angular difference between first eigenvectors


phi_spe= zeros(1,length(t));

pv1= S1.EigVectors(:,1);
for ii = 1 : length(t)
    pvi= S_SQ(ii).EigVectors(:,1);
    phi_spe(ii)= mod(acos(abs(pv1'*pvi)), 2*pi);
end

% Determine the limits of the y-axis

minPhi = min(min(phi_loge),min(phi_spe));
maxPhi = max(max(phi_loge),max(phi_spe));

figure
plot(t,phi_loge,'-ob','linewidth',2);
axis([0 1 minPhi maxPhi]);
hold on
plot(t,phi_spe,'-or','linewidth',2);
legend('Log-E','SQ with HA');
title('Orientation evolution - Angular difference');    



