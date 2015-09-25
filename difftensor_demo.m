%% DEMO 1: INTERPOLATION BETWEEN TWO OR FOUR TENSORS

% This function will show how to compute an interpolating curve between two
% tensors, and how to interpolate a 2x2 tensor field.


% ________________________________________________________
% Copyright (C) 2013 University of Liege, Belgium
    
% Written by A. Collard & C. Phillips, 2013.
% Dept of Electrical Engineering and Computer Science &
% Cyclotron Research Centre, University of Liege, Belgium
% ________________________________________________________

clear all
close all
% clc

%% Interpolations of two tensors

% Construction of the two 'extremities' tensors
alpha= (1/360)*2*pi;

R= [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
S1= difftensor(R*diag([10 1 1])*R');

alpha= (63/360)*2*pi;

R= [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
S2= difftensor(R*diag([40 4 1])*R');

S = difftensor;
S(1) = S1;
S(2) = S2;

% t : interpolation parameter

t = 0:0.125:1;

%% First case : with no rescaling of the orientation interpolation for SQ.

% Comparison of the Log-Euclidean and Spectral-Quaternion interpolation

display('The first example will perform an interpolation between 2 tensors,');
display('using the Log-Euclidean method and the Spectral-Quaternion method');
display('with no rescaling of the orientation interpolation.');
display('For each curve, the evolution of the principal features of');
display('the tensors are also computed and illustrated.');
display('See the file to understand the used commands.');

% Slog : will contain the Log-Euclidean interpolation
S_LogE = difftensor;

% Sspe : will contain the Spectral quaternion interpolation
S_SQ = difftensor;

for ii = 1 : length(t)
    w = [(1-t(ii)) t(ii)];
    S_LogE(ii) = wmean(S,w,'LogE');
    S_SQ(ii) = wmean(S,w,'SQ','no');
end

% Illustration of interpolations
graph_display(S_LogE,0.3);
title('Log-Euclidean interpolation');
axis off

graph_display(S_SQ,0.3);
title('Spectral-Quaternion interpolation - no rescaling');
axis off

% Evolution of the determinant
det_loge = getDet(S_LogE);
det_spe = getDet(S_SQ);

% Determine the limits of the y-axis
minDet = min(min(det_loge),min(det_spe));
maxDet = max(max(det_loge),max(det_spe));

% Evolution of the Fractional Anisotropy
FA_loge = getFA(S_LogE);
FA_spe = getFA(S_SQ);

% Determine the limits of the y-axis
minFA = min(min(FA_loge),min(FA_spe));
maxFA = max(max(FA_loge),max(FA_spe));

% Evolution of the Hilbert Anisotropy
% HA_loge = getHA(S_LogE);
% HA_spe = getHA(S_SQ);
HA_loge = [S_LogE.HA];
HA_spe = [S_SQ.HA];

% Determine the limits of the y-axis
minHA = min(min(HA_loge),min(HA_spe));
maxHA = max(max(HA_loge),max(HA_spe));

% Evolution of the angular difference between first eigenvectors
phi_loge = zeros(1,length(t));
pv1 = S1.EigVectors(:,1);

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
S_SQ = difftensor;
for ii = 1 : length(t)
    w = [(1-t(ii)) t(ii)];
    S_SQ(ii) = wmean(S,w,'SQ','kappa');
end

% Illustration of interpolations
graph_display(S_SQ,0.3);
title('Spectral-Quaternion interpolation - \kappa rescaling');
axis off

% Evolution of det, HA, FA are identical than with no rescaling.
% Evolution of the angular difference between first eigenvectors

phi_spe = zeros(1,length(t));
pv1 = S1.EigVectors(:,1);
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
S_SQ = difftensor;
for ii = 1 : length(t)
    w = [(1-t(ii)) t(ii)];
    S_SQ(ii) = wmean(S,w,'SQ','HA');
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

display('Execution paused. Press any key to continue.');
pause

%% Interpolations of four tensors
display('In the first example, the four tensors at the corner of the field');
display('are averaged to fill the field. The weights are related to the location');
display('of the tensors. No rescaling is used for the SQ method.');

% Construction of the four corner tensors 
Sc = difftensor;
Sc(1) = eye(3);

alpha = pi/4;
R = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1];
Sc(2)= R*diag([2 0.2 0.1])*R';

alpha = pi/3;
R = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1];
Sc(3) = R*diag([2 1.9 0.1])*R';

alpha = 0;
R = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1];
Sc(4) = R*diag([2 1.9 0.1])*R';

% Initialization 
Ssq = difftensor; % will contain the SQ average
Sloge = difftensor; % will contain the Log-E average

% Fill the array of tensors
for ii = 1 : 10
    for jj = 1:10
        % Computation of the weights (depend upon the location in the
        % array)
        terme_i = ones(4,1)-[0 0 1 1]'+ ((ii-1)/9)*[-1 -1 1 1]';
        terme_j = ones(4,1)-[0 1 0 1]'+ ((jj-1)/9)*[-1 1 -1 1]';
        w = terme_i.*terme_j;
        
        Ssq(ii,jj) = wmean(Sc,w); % SQ method, no rescaling (default values)
        Sloge(ii,jj)= wmean(Sc,w,'LogE'); % Log-E method
    end
end

% Fields of tensors
graph_display(Ssq,0.4,'Hacol');
view(0,90);
axis off,
title('Spectral-Quaternion - no rescaling');

graph_display(Sloge,0.4,'Hacol');
view(0,90);
axis off
title('Log-Euclidean');

% Fields of Principal Direction Diffusion (PDD)
vect_display(Ssq);
view(0,90);
axis off
title('Spectral-Quaternion - no rescaling');

vect_display(Sloge);
view(0,90);
axis off
title('Log-Euclidean');

% Contour maps of the anisotropy
FA_SQ_no = getFA(Ssq);
FA_LE = getFA(Sloge);

figure
[C,h] = contour(FA_LE,'linewidth',2);
clabel(C,h,'FontSize',11);
title('FA LogE');

figure
[C,h] = contour(FA_SQ_no,'linewidth',2);
clabel(C,h,'FontSize',11);
title('FA SQ - no rescaling');

display('Execution paused. Press any key to continue.');
pause
%% SQ - autre ponderation
display('In the second example, the kappa weights are used.');

Ssq_kappa = difftensor;
for ii = 1 : 10
    for jj = 1:10
        terme_i = ones(4,1)-[0 0 1 1]'+ ((ii-1)/9)*[-1 -1 1 1]';
        terme_j = ones(4,1)-[0 1 0 1]'+ ((jj-1)/9)*[-1 1 -1 1]';
        w = terme_i.*terme_j;
        Ssq_kappa(ii,jj) = wmean(Sc,w,'SQ','kappa');
    end
end

graph_display(Ssq_kappa,0.4,'Hacol');
view(0,90);
axis off
title('SQ- \kappa rescaling');

vect_display(Ssq_kappa);
view(0,90);
axis off
title('SQ with \kappa rescaling');

display('Execution paused. Press any key to continue.');
pause

%% Rotation and scaling of each tensor of an array

display('Using the rotate function, a rotation can be applied to each tensor');
display('of an array. The scale function enables to rescale the tensors (all');
display('the tensors are multiplied by a given scalar s');

R = [0 -1 0;1 0 0; 0 0 1];

Snew = rotate(Ssq_kappa, R);
Snew = scale(Snew,2);

graph_display(Snew,0.4,'Hacol');
view(0,90);
axis off
title('Rotation and scaling of each tensor of the preceding field');

vect_display(Snew);
view(0,90);
axis off
title('Rotation and scaling of each tensor of the preceding field')








