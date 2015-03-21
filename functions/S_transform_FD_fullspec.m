function [S,FVEC] = S_transform_FD_fullspec(s,dt,k)
%
%
% This routine follows the example given in Appendix A of C.R. Pinnegar's
% paper "Time-frequency and time-time filtering with the S-transform and
% TT-transform", Digital Signal Processing 15 (2005) 604-620,
% doi:10.1016/j.dsp.2005.02.002. I have only changed some notation for ease
% of understanding the code and added a few checks.
% 
% USAGE: [S,FVEC] = S_transform_FD_fullspec(st,dt,k)
% 
% INPUT:
%   s    = time series
%   dt   = sample interval (s) (default=1)
%   k    = integer value for number of periods to make width of Gaussian
%   (default=2)
% OUPUT:
%   S    = time-frequency matrix of the complex S-transform coefficients
%           (row=f,col=t)
%   FVEC = frequency vector of associated frequencies (symmetric containing
%           both positive and negative frequencies
% 
% EXAMPLE:
%
% fmax  = 50;                   % [Hz] Nyquist frequency for S-transform
% dt    = 1/2/fmax;             % (s) standard dt from sampling theory
% Tmax  = 5;                    % (s) trace length
% npts  = floor(Tmax/dt)+1;     % number of time samples
% tvec  = (0:npts-1).*dt;       % time vector
% 
% % make trace
% delay = 2;                            % ricker wavelet delay
% fc    = 10;                           % center frequency
% ampc  = 1;                            % ricker amplitude 
% h     = rickerTD(ampc,fc,delay,tvec); % my trace
% 
% [S,fvec] = S_transform_FD_fullspec(h,dt); % frequency domain S-transform
%
% imagesc(tvec,fvec,abs(S));
% set(gca,'YDir','Normal');  c=colorbar; ylabel('Frequency (Hz)');
% title('S-transform'); xlabel('Time (s)');
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

%--------------------------------------------------------------------------
% set defaults
if nargin < 2
    dt = 1; % (s)
    k  = 2; % default is for 2 period window
elseif nargin < 3
    k  = 2; % default is for 2 period window
end
%--------------------------------------------------------------------------
npts=numel(s); % number of points in trace
% make an even number of points
oddflag = 0;
if mod(npts,2)~=0
    s       = s(1:end-1);
    npts    = numel(s);
    oddflag = 1;
end

mid = npts/2;
idx = [0:mid-1, -mid:-1];

FVEC = idx./dt/npts;     % frequency vector
S    = zeros(npts,npts); % allocate S-transform matrix

st = fft(s);     % trace FFT
st = [st st st]; % concatenate trace FFT for sliding to avoid self-aliasing

% loop through frequencies using parallel or serial loops
%
% Parallel option commented out by Haney 4/15/2014
%
%if matlabpool('size')==0
%     disp('Running S-transform in serial');
    for f=1:npts
        % moving Gaussian window function
        W = exp( -2*pi^2*(idx.^2)./ ( (k/2)*(f-mid-1).^2 ) );
        % slide through FFT and apply window
        S(f,:) = ifft( st(f+mid:f+mid+npts-1).*W);
    end
%else
%     disp('Running S-transform in paralell');
%     parfor f=1:npts
%         % moving Gaussian window function
%         W = exp( -2*pi^2*(idx.^2)./ ( (k/2)*(f-mid-1).^2 ) );
%         % slide through FFT and apply window
%         S(f,:) = ifft( st(f+mid:f+mid+npts-1).*W);
%     end
% end

% this is the f=0 s-transform which is just the mean.
S(mid+1,:) = ifft( st(npts+1:2*npts).*[1 zeros(1,npts-1)] );

if oddflag % add zeros to last time sample is odd npts
    S(:,npts+1) = complex(zeros(npts,1));
end

FVEC=fftshift(FVEC); % arrange frequency vector to be symmetric

%--------------------------------------------------------------------------
return