function st = S_transform_inverse_fullspec(S,fvec)
%
% This script calculates the inverse S-transfrom of a signal by integrating
% over the time window and inverse Fourier transforming.
%
% USAGE: st = S_transform_inverse_fullspec(S,fvec)
%
% INPUT:
%   S    = time-frequency matrix of the complex S-transform coefficients
%           (row=f,col=t)
%   fvec = Frequency vector (symmetric contain both positive and negative
%           frequencies)
% OUPUT:
%   st   = time series trace after inverse S-transform
%
% EXAMPLE:
%
% Time and frequency sampling parameters
% fmax = 100;                   % [Hz]  Maximum frequency in S-transform
% dt   = 1/(2*fmax);            % [s]   Sample interval in time domain
% T    = 1;                     % [s]   Length of time domain trace
% df   = 1/T;                   % [Hz]  Sample interval in frequency domain
% npts = round(1/dt/df);        % [int] From sampling theorem
% tvec = 0:dt:(npts-1)*dt;      % [s]   Time vector
% 
% % Properties of Ricker wavelet
% fc    = 25;   % [Hz] center frequency
% tau   = 0.35; % [s] time delay
% r_amp = 1;    % [a.u.] ricker wavelet amplitude
% 
% % Ricker wavelet in time and frequency domains
% st    = rickerTD(r_amp,fc,tau,tvec); % Time domain
% 
% % S-transform of time domain trace and the inverse S-transform
% [sS,fvec] = S_transform_FD_fullspec(st,dt);    % matrix(freq,time)
% st_inv    = S_transform_inverse_fullspec(sS,fvec); % inverse trace
% 
% % Plot and compare different approaches
% 
% figure;
% plot(tvec,st,'-k'); hold on;
% plot(tvec,st_inv,'--r'); legend('Original','Inverse');
% ylabel('Amplitude (a.u.)');  xlabel('Time (s)'); title('Ricker wavelet');
% 
% % plot S-transform of Ricker wavelet trace
% figure;
% imagesc(tvec,fvec,abs(sS)); xlabel('Time (s)'); set(gca,'YDir','Normal');
% ylabel('Frequency (Hz)'); title('S-transform'); c=colorbar;
%
%
% DISCLAIMER:
% The accompanying program is intended for the use by members of the
% applied geophysics group of TU Delft only. THE PROGRAM IS PROVIDED ON AN
% "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, EITHER
% EXPRESS OR IMPLIED INCLUDING, WITHOUT LIMITATION, ANY WARRANTIES OR
% CONDITIONS OF TITLE, NON-INFRINGEMENT, MERCHANTABILITY OR FITNESS FOR A
% PARTICULAR PURPOSE.
%
% AUTHOR:
% Dylan Mikesell, mikesell@mit.edu, January 2014
%
% Code mostly from an example by Robert Glenn Stockwell's function
% inverse_st.m ($Revision: 1.0 $  $Date: 2004/10/10  $)
%
% Reference is "Localization of the Complex Spectrum: The S-Transform" from
% IEEE Transactions on Signal Processing, vol. 44., number 4, April 1996,
% pages 998-1001.

% Dimensions of S-transform matrix
[nf,npts] = size(S);

% Check that matrix is in correct dimensions
if nf ~= numel(fvec)
    S = S';
    [nf,npts] = size(S);
end

% integrate over time-axis to get FFT spectrum
spec_full = sum(S,2); 

% the time series is the inverse fft of this
ts = ifft(fftshift(spec_full));

% and take the real part, the imaginary part should be zero.
st = real(ts);

% if odd number of points add 1 sample to match original trace
if rem(npts,2) ~= 0
    st(end+1)=0;
end

return
