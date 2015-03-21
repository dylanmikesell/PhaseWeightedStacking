
% demo by Haney 4/15/2014

clear all

% time this
tic

addpath('functions/');

% load 10 Redoubt LPs from April 4, 2009
a = load('10evts_rd02z.txt');
[a1 a2] = size(a);

% number of events to use
Nv = 10;
% time sample rate
dt = 1/50;
% time vector
tvec = [0:(a2-1)]*dt;
% power of stacking as in Schimmel et al. (2010; GJI)
pwr = 2;

% S-transform of linear stack
[stranl,fvec] = S_transform_FD_fullspec(mean(a(1:Nv,:)),dt);

% make phase weight
sumr = zeros(a2,a2);
for ii=1:Nv
    % S-transform of ii-th seismogram
    [stran,fvec] = S_transform_FD_fullspec(a(ii,:),dt);
    % Equation 6 in Schimmel et al. (2010; GJI)
    sumr = sumr + (stran./abs(stran)).*exp(i*2*pi*(fvec'*tvec));
    ii
end
% Equation 6 in Schimmel et al. (2010; GJI)
sumr = abs(sumr/Nv).^pwr;

% apply phase weight to linear stack to get phase-weighted-stack
% Equation 7 in Schimmel et al. (2010; GJI)
stranpws = sumr.*stranl;

% transform back to the time domain
pws = S_transform_inverse_fullspec(stranpws,fvec);

% make a plot
figure
fsize = 16;
subplot(2,1,1)
plot(tvec,mean(a(1:Nv,:))); axis([ 0 81.9 -3*(10^-7) 3*(10^-7) ])
set(gca,'Fontsize',fsize,'FontWeight','bold');
ylabel(' Amp. (m/s) ','FontSize',fsize,'FontWeight','bold');
%xlabel(' Time (s) ','FontSize',fsize,'FontWeight','bold');
title(' Linear Stack ','FontSize',fsize,'FontWeight','bold');

subplot(2,1,2)
plot(tvec,pws); axis([ 0 81.9 -3*(10^-7) 3*(10^-7) ])
set(gca,'Fontsize',fsize,'FontWeight','bold');
ylabel(' Amp. (m/s) ','FontSize',fsize,'FontWeight','bold');
xlabel(' Time (s) ','FontSize',fsize,'FontWeight','bold');
title(' Phase-Weighted Stack ','FontSize',fsize,'FontWeight','bold');

orient landscape
print(gcf,'-dpsc2','pws_example_redoubt.ps');

toc

