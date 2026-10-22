function [H_f,F_f,phase_corr]=gtf_pp(G_t,Gradient_Amplitude,sig,DwellTime,cal_fac,carrier,Tukey_slope_len,rel_wnd_len ...
                                        ,phase_corr_man,integ_fd)
e = exp(1) ;

size_sig  = size(sig) ;
nmax      = size_sig(2);

% Preparing the sampling rate conversion
% Sample rate of the excitation waveform (chirp) is always 100 kHz (corresponds to a dwell time of 10 µs).
if 1
  % down scaling (Exitation waveform  is decimated, faster)
  fs       = 1/DwellTime;
  upsf     = 1;
  dnsf     = round(DwellTime/10e-6) ;
  adc_len  = size_sig(1) ;
  switch(dnsf)
      case ( 1)
        fft_len = 200000; % 0.5 Hz
      case ( 2)
        fft_len = 100000; % 0.5 Hz
      case ( 4)
        fft_len =  50000; % 0.5 Hz
      case ( 5)
        fft_len =  40000; % 0.5 Hz
      case ( 8)
        fft_len =  25000; % 0.5 Hz
      case (10)
        fft_len =  20000; % 0.5 Hz
      otherwise
        error('Dwelltime not supported');
  end
else
  % upscaling (ADC data is interpolated, slower)
  fs      = 100000         ;
  upsf    = DwellTime*fs   ;
  dnsf    =      1         ;
  adc_len = size_sig(1)*upsf    ;
  fft_len = 200000         ;   % 0.5 Hz  (>= adc_len)
end


fc_idx    = fft_len/2+1;
freq      = ( (-fs/2):(fs/fft_len):(fs/2-fs/fft_len) );
F_f       = freq(fc_idx:1:end)';


% Create FFT-Window
wnd_len    = round(adc_len*rel_wnd_len) ;
TUKEY      = tukeywin(wnd_len, Tukey_slope_len/(wnd_len/fs/2));  % 0-> rect  ... 1-> Hann
fft_window = [TUKEY ; zeros((adc_len-wnd_len),1) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Excitation waveform G(t) resp. G(f) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G_t = G_t / max(abs(G_t))      ; % normalization to 1
G_t = G_t * Gradient_Amplitude ;
G_t = G_t(1:dnsf:end)          ; %
G_t =[G_t ;  zeros(fft_len-length(G_t),1)] ; % Zero-Padding

% FFT
G_f    = fftshift(fft(G_t/adc_len)) ;
G_f    = G_f( fc_idx:1:end) *2      ; %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Ch. 1 Carrier signal from ADC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AM cancellation, normalization to 1
carrier = carrier./abs(carrier);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ch. n DSB data from ADC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply calibration factor
sig = sig * cal_fac;

% Mixing to baseband
sig_c  = sig .* conj(carrier);

% Averaging of repetitions
if nmax==1
  sig_ca =     (sig_c  ) ;
else
  sig_ca = mean(sig_c.').';
end

% Upsampling, if requested
if upsf>1    % upsampling by upsf
  sig_ca  = resample(sig_ca,upsf,1);
end

% Signal fenstern
kwin    = 1./mean(fft_window)   ;
sig_caw = sig_ca .* fft_window  ;

% remove_DC cmponent
sig_caw = sig_caw -mean(real(sig_caw))*kwin*fft_window -1i*mean(imag(sig_caw))*kwin*fft_window;

% phase correction
if phase_corr_man~=0
  % manual phase correction
  phase_corr_deg = phase_corr_man;
else
  % data driven phase correction
  phi = zeros(1,198);
  r2i = zeros(1,198);
  for n=1:1:198 % 1 deg steps
    phi(n) = (-99+1.0*(n-1));
    tmpxx  = abs(sig_caw).*e.^(1i*(angle(sig_caw)-phi(n)*pi/180));
    r2i(n) = ( (sum(imag(tmpxx.').^2)) / (sum(real(tmpxx.').^2)) );
  end
  r2i_mx_idx = find(r2i==min(r2i)) ;
  r2i_mx_idx = r2i_mx_idx(1) ;
  if 1
    v   = r2i([r2i_mx_idx-2 r2i_mx_idx-1 r2i_mx_idx r2i_mx_idx+1 r2i_mx_idx+2]) ; % neighborhood of the minimum
    x   = -2:1:2;
    xq  = -2:0.001:2 ; % Interpolation in 0.01 deg steps
    vq  = interp1(x,v,xq,'spline');
    vq_mx_idx = find(vq==min(vq));
    vq_mx_idx = vq_mx_idx(1) ;
    phase_corr_deg = phi(r2i_mx_idx) + xq(vq_mx_idx)   ;
  end
end
phase_corr = phase_corr_deg *pi/180;
sig_caw    = abs(sig_caw).*e.^(1i*(angle(sig_caw)-phase_corr));

% FFT
sig_cawz = [sig_caw ; zeros(fft_len-length(sig_caw),1)] ; % Zero-Padding
dsb_f    = fftshift(fft(sig_cawz)/adc_len);

% Combining  sidebands
lsb_f = flip(conj(dsb_f(2:1:fc_idx))) ;
usb_f =           dsb_f(fc_idx:1:end) ;
P_f   = lsb_f + usb_f                 ;


f0   = 1000/(2*pi) ; %159.15 Hz
fvec = -1i*(0:(fs/fft_len):(fs/2-fs/fft_len))'/f0 ;
% Integration, if requested
if integ_fd >=(+1)
  for n=1:1:integ_fd
    P_f = P_f ./ fvec ;
  end
end
% Differentiation, if requested
if integ_fd <=(-1)
  for n=1:1:-integ_fd
    P_f = P_f .* fvec ;
  end
end

%disp(['angle(P_f(100)) = ' num2str(angle(P_f(50))*180/pi)  ])


% P(f) / G(f) = H(f) Transfer function
H_f = P_f ./  G_f  ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


