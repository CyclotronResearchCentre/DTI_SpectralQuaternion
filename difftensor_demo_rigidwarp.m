%% DEMO 3: Rigid warp example

% In this example, we show how to use the tools defined in the 'difftensor'
% object to perform a (spatial) transformation of a DT image.

% A simple rotation of the image is presented, but this method could be
% used for any linear transformation.

% Possible variations are :
% - a different reorientation strategy (the presented one is the 'Finite
% Strain')
% - other rescaling for the SQ method
% - linear transformation including translation, shear, and so on.

% ________________________________________________________
% Copyright (C) 2014 University of Liege, Belgium
    
% Written by A. Collard & C. Phillips, 2014.
% Dept of Electrical Engineering and Computer Science &
% Cyclotron Research Centre, University of Liege, Belgium
% ________________________________________________________

clear all
close all
clc


%% Loading of the example image (slice of a real DTI)

a = load('slice_true_im1.mat');
Im_t = a.Im_t; % Im_t is a cell array : each array contains a 3x3 symmetric 
                % positive definite matrices (a diffusion tensor)
d = difftensor(Im_t);   % the difftensor function can handle such a cell array and
                        % transform it in a 'difftensor' array.

% Illustration of the 'true' image - only a small part of the image is
% shown.

graph_display(d(12:33,7:30),5);  
view(0,90);
axis off
title('True image');

%% Rigid rotation of the whole image (no rescaling)

% As an example, a rigid rotation of the whole image will be performed. Any
% linear transformation can be performed using this method. The linear
% transformation should be given as a 3x3 matrix M, and a translation
% vector, t.

alpha = 2*pi/15;
% All the image will be rotated of 2pi/15 rad.

M = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
% No translation of the image
t = [0 0 0]';

% Find R, the rotation matrix used for the reorientation
% Finite Strain Strategy (see 'Alexander et al, Spatial transformations of
% diffusion tensor magnetic resonance images, IEEE Trans. on Med. Imaging,
% 20, (2001)'.

R= (M*M')^(-1/2)*M; % this formula applies to any linear transfo matrix.


% "zero" difftensor, equivalent to a zero in a scalar image

F= difftensor(0.005*eye(3));

% Initialization of the transformed images

ds = size(d);
d_trans_SQ= difftensor;
d_trans_LE = difftensor;
l = ds(1);
lm = ceil(l/2);
m= ds(2);
mm = ceil(m/2);
if length(ds)<3
    h=1;
else
	h = ds(3);
end
hm= ceil(h/2);


for ii=1:l
    for jj=1:m
        for kk=1:h
            
            % find coordinates 
            x= [ii-lm,jj-mm,kk-hm]';
            M1x= M*x;
            Mtx= M1x+t;
            
             % the center of the image is used as the 'zero' coordinate.
            Mtx = Mtx+ [lm mm hm]';
            
            % Interpolation is used to compute the value of the tensor at
            % the transformed coordinates
            
        i1= floor(Mtx(1)); i2= ceil(Mtx(1));
        
        j1= floor(Mtx(2)); j2= ceil(Mtx(2));
        
        k1= floor(Mtx(3)); k2= ceil(Mtx(3));
        
        if i1<=0 ||i2> l|| j1<=0 ||j2> m || k1<=0 || k2>h % outside the image
            d_trans_SQ(ii,jj,kk)=F;
            d_trans_LE(ii,jj,kk)=F;
        elseif Mtx(1)>=1 && Mtx(1)<=l && Mtx(2)>=1 && Mtx(2)<=m && Mtx(3)>=1 && Mtx(3)<=h %inside
            
            % the interpolation is based on the 8 nearest neighbours
        ind_i= [i1 i1 i2 i2 i1 i1 i2 i2]';
        ind_j= [j1 j2 j1 j2 j1 j2 j1 j2]';
        ind_k= [k1 k1 k1 k1 k2 k2 k2 k2]';
    
        ind= (ind_k-ones(8,1))*m*l+ (ind_j-ones(8,1))*l+ ind_i;
    
        
        L= d(ind);
        
        
        xx= Mtx(1)-i1;
        yy= Mtx(2)-j1;
        zz= Mtx(3)-k1;
        terme_i= ones(8,1)-[0 0 1 1 0 0 1 1]'+ xx*[-1 -1 1 1 -1 -1 1 1]';
        terme_j= ones(8,1)-[0 1 0 1 0 1 0 1]'+ yy*[-1 1 -1 1 -1 1 -1 1]';
        terme_k= ones(8,1)-[0 0 0 0 1 1 1 1]'+ zz*[-1 -1 -1 -1 1 1 1 1]';

        w= terme_i.*terme_j.*terme_k;
        
        WAsq = mean(L,w,'SQ');
        d_trans_SQ(ii,jj,kk)= WAsq;
        
        WA_le = mean(L,w,'LogE');
        d_trans_LE(ii,jj,kk)= WA_le;
                   
        
        end
        
    
        end
    end
end

% The reorientation applies to the whole image
S_Sq= rotate(d_trans_SQ,R);
Sloge = rotate(d_trans_LE,R);

% Illustration of the results
graph_display(S_Sq(12:33,7:30),5); 
view(0,90);
axis off
title('SQ- \kappa rescaling');

graph_display(Sloge(12:33,7:30),5); 
view(0,90);
axis off
title('Log-E');




