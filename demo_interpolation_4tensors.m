%% DEMO 2: INTERPOLATION BETWEEN FOUR TENSORS

% In this demo, we show how to compute the mean of many tensors (4 in this
% case). The evolutions of the features of tensors are also shown.

% ________________________________________________________
% Copyright (C) 2013 University of Liege, Belgium
    
% Written by A. Collard & C. Phillips, 2013.
% Dept of Electrical Engineering and Computer Science &
% Cyclotron Research Centre, University of Liege, Belgium
% ________________________________________________________

clear all
close all
clc

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
        terme_i= ones(4,1)-[0 0 1 1]'+ ((ii-1)/9)*[-1 -1 1 1]';
        terme_j= ones(4,1)-[0 1 0 1]'+ ((jj-1)/9)*[-1 1 -1 1]';
        w= terme_i.*terme_j;
        
        Ssq(ii,jj) = wmean(Sc,w,'SQ','no'); % SQ method, no rescaling 
        Sloge(ii,jj)= wmean(Sc,w,'LogE'); % Log-E method
    end
end

% Fields of tensors
graph_display(Ssq,0.4,'Hacol');
view(0,90);
axis off
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
title('FA LE');

figure
[C,h] = contour(FA_SQ_no,'linewidth',2);
clabel(C,h,'FontSize',11);
title('FA SQ- no rescaling');

display('Execution paused. Press any key to continue.');
pause
%% SQ - autre ponderation

Ssq_kappa = difftensor;
for ii = 1 : 10
    for jj = 1:10
        terme_i= ones(4,1)-[0 0 1 1]'+ ((ii-1)/9)*[-1 -1 1 1]';
        terme_j= ones(4,1)-[0 1 0 1]'+ ((jj-1)/9)*[-1 1 -1 1]';
        w= terme_i.*terme_j;
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








