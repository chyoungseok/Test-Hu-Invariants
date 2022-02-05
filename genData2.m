%% generate pattern
close all; clear; clc;

% x1 = [0 0 0 0 0 0 0 0 0 0 0;
%       0 1 1 1 0 0 0 0 0 0 0;
%       0 1 0 1 0 0 0 0 0 0 0;
%       0 1 0 1 1 1 1 1 1 1 0;
%       0 0 0 0 0 0 0 0 0 0 0];
% 
% x2 = [0 0 0 0 0 0 0 0 0 0 0;
%       0 0 0 0 1 1 1 0 0 0 0 ;
%       0 0 0 0 1 0 1 0 0 0 0;
%       0 1 1 1 1 0 1 1 1 1 0;
%       0 0 0 0 0 0 0 0 0 0 0];
% 
% x3 = [0 0 0 0 0 0 0 0 0 0 0;
%       0 0 0 0 0 0 0 1 1 1 0;
%       0 0 0 0 0 0 0 1 0 1 0;
%       0 1 1 1 1 1 1 1 0 1 0;
%       0 0 0 0 0 0 0 0 0 0 0];
% 
% subplot(1, 3, 1); imshow(x1, 'InitialMagnification', 10000); hold on; title("x1", 'FontSize', 25)
% subplot(1, 3, 2); imshow(x2, 'InitialMagnification', 10000); hold on; title("x2", 'FontSize', 25)
% subplot(1, 3, 3); imshow(x3, 'InitialMagnification', 10000); hold on; title("x3", 'FontSize', 25)
% 
% X = [x1; x2; x3];

template = zeros(5,50);
X = {};
pat_len = 4; % length of pattern where the values are 1

iter = 0; % counting the iteration of while structure
while iter*pat_len+pat_len-1 < 50
    x = template;
    x(2, 2+iter*pat_len:2+iter*pat_len+pat_len-1) = 1;
    x(3, 2+iter*pat_len) = 1; x(3, 2+iter*pat_len+pat_len-1)=1;
    x(4, 2:2+iter*pat_len) = 1; x(4, 2+iter*pat_len+pat_len-1:end-1) = 1;
    X{iter+1} = x;

    iter = iter + 1;
end

for i = 1 : iter
    subplot(subplot(iter/2, 2, i))
    if rem(i,2) == 1
        temp = X{(i+1)/2};
        imshow(temp, 'InitialMagnification', 100000); hold on; 
        title(['X', num2str((i+1)/2)], 'FontSize',20)
    else
        temp = X{i/2+6};
        imshow(temp, 'InitialMagnification', 100000); hold on; 
        title(['X', num2str(i/2+6)], 'FontSize',20)
    end
end

% X = template;
% X(2, 2:2+pat_len-1) = 1;
% X(3, 2) = 1; X(3, 2+pat_len-1) = 1;
% X(4,2) = 1; X(4, 2+pat_len-1: end-1) = 1;
% imshow(X, 'InitialMagnification', 10000);

moments = [];
moments_norm = [];
Y = {X{1}, X{6}, X{12}};
N = iter;
%% 
for i = 1 : iter   
    image = X{i};
    [height, width] = size(image);
    
    %% Calculate the required parameters
    % define a co-ordinate system for image 
    xgrid = repmat((-floor(height/2):1:ceil(height/2)-1)',1,width);
    ygrid = repmat(-floor(width/2):1:ceil(width/2)-1,height,1);
    
    [x_bar, y_bar] = centerOfMass(image,xgrid,ygrid);
    
    % normalize coordinate system by subtracting mean
    xnorm = x_bar - xgrid;
    ynorm = y_bar - ygrid;
    
    %% Calculate the central moments
    % central moments
    mu_11 = central_moments( image ,xnorm,ynorm,1,1);
    mu_20 = central_moments( image ,xnorm,ynorm,2,0);
    mu_02 = central_moments( image ,xnorm,ynorm,0,2);
    mu_21 = central_moments( image ,xnorm,ynorm,2,1);
    mu_12 = central_moments( image ,xnorm,ynorm,1,2);
    mu_03 = central_moments( image ,xnorm,ynorm,0,3);
    mu_30 = central_moments( image ,xnorm,ynorm,3,0);
    
    %% Calculate Hu's Invariant moments
    %central_moment = [mu_11, mu_20, mu_02, mu_21, mu_12, mu_03, mu_30];
    %calculate first 8 hu moments of order 3
    I_one   = mu_20 + mu_02;
    I_two   = (mu_20 - mu_02)^2 + 4*(mu_11)^2;
    I_three = (mu_30 - 3*mu_12)^2 + (mu_03 - 3*mu_21)^2;
    I_four  = (mu_30 + mu_12)^2 + (mu_03 + mu_21)^2;
    I_five  = (mu_30 - 3*mu_12)*(mu_30 + mu_12)*((mu_30 + mu_12)^2 - 3*(mu_21 + mu_03)^2) + (3*mu_21 - mu_03)*(mu_21 + mu_03)*(3*(mu_30 + mu_12)^2 - (mu_03 + mu_21)^2);
    I_six   = (mu_20 - mu_02)*((mu_30 + mu_12)^2 - (mu_21 + mu_03)^2) + 4*mu_11*(mu_30 + mu_12)*(mu_21 + mu_03);
    I_seven = (3*mu_21 - mu_03)*(mu_30 + mu_12)*((mu_30 + mu_12)^2 - 3*(mu_21 + mu_03)^2) + (mu_30 - 3*mu_12)*(mu_21 + mu_03)*(3*(mu_30 + mu_12)^2 - (mu_03 + mu_21)^2);
    I_eight = mu_11*(mu_30 + mu_12)^2 - (mu_03 + mu_21)^2 - (mu_20 - mu_02)*(mu_30 + mu_12)*(mu_21 + mu_03);
    
    %% Apply log, and view the results.
    hu_moments_vector = [I_one, I_two, I_three,I_four,I_five,I_six,I_seven,I_eight];
    hu_moments_vector_norm= -sign(hu_moments_vector).*(log10(abs(hu_moments_vector)));
    moments = [moments; hu_moments_vector];
    moments_norm = [moments_norm; hu_moments_vector_norm];
end

figure;
for m_i = 1 : 8
    subplot(4, 2,m_i)
    plot(moments(:,m_i), 'x-'); hold on;
    title(['I', num2str(m_i)], 'FontSize', 20);
end
% legend(["I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8"]);
% legend(["I3", "I4", "I5", "I6", "I7", "I8"]);