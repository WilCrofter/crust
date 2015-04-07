% This is Poisson solver
% Last edit?: 4/06/2015
clear all;
Probe_L = 58; %mm
N_elements = 128;
dis_trans_receiv = 36;%mm
typ_dis_trans_receiv = 36;%mm

Speed_Backgrnd = 1010;%m/s

% Error (mm)
Error_max = 2; % mm

if Probe_L/2 ~= floor(Probe_L/2)
    error('Probe_L should be even')
end
if N_elements/2 ~= floor(N_elements/2)
    error('Probe_L should be even')
end

% Grid definition
gridx0 = [-Error_max-Probe_L/2:1:Probe_L/2+Error_max]';
gridy0 = [0;1];
gridz0 = [-Error_max:1:typ_dis_trans_receiv+Error_max]';

I_bg = 1/Speed_Backgrnd*ones(length(gridz0)-1,length(gridx0)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit-receive geometry
pos_probe_elem = linspace(-(Probe_L/N_elements)*(N_elements/2-0.5),(Probe_L/N_elements)*(N_elements/2-0.5),N_elements)';
% load('Trans_pose2.mat')
% trans_pose = Trans_pose;
trans_pose = [pos_probe_elem,zeros(N_elements,1)];
recev_pose = [pos_probe_elem,dis_trans_receiv*ones(N_elements,1)];
% 
% trans_pose2 = [pos_probe_elem,zeros(N_elements,1)];
% recev_pose2 = [pos_probe_elem,dis_trans_receiv*ones(N_elements,1)];

%%%%%%%%%%%%%%%%%%%%SM & tof%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_t2 = getSysMat_multi_src_pos (trans_pose2, recev_pose2, gridx0, gridz0); 

SM_t = getSysMat_multi_src_pos (trans_pose, recev_pose, gridx0, gridz0); 
load('TOF_exp_seg');
tof = TOF_seg_1d;
% load('TOF_Lei.mat')
% tof = tof_Lei;
% t=tof;
% t = reshape(t,128,128);
% t = t';
% tof = t(:);
load('TOF_exp_bg.mat')
tof_bg = TOF_time_1d;
% t = tof_bg;t = reshape(t,128,128);
% t = t';
% tof_bg = t(:);
% load('TOF_bg_auto.mat')
% tof_bg = t;
% t = tof_bg;t = reshape(t,128,128);
% t = t';
% tof_bg = t(:);



N_iter = 100;
b_pre = I_bg(:);
tof = tof_bg - tof;
% writerObj = VideoWriter(['Movie_EM','.avi']);
% open(writerObj);
for m = 1 : N_iter
    SUM = SM_t'*(tof./(SM_t*b_pre));
    b = b_pre .* SUM ./ sum(SM_t,1)';
    imagesc(1./reshape(I_bg(:)-b,size(I_bg))); colorbar;title(['iteration ',num2str(m)]);
%     imagesc(1./reshape(b,size(I_bg))); colorbar;title(['iteration ',num2str(m)]);
    drawnow;
    b_pre = b;
%     Frame = getframe(gcf);
%     writeVideo(writerObj,Frame); 
end
% close(writerObj);
