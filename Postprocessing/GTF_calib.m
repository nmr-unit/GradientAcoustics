clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The calib_data_x.mat files were created from the corresponding Siemens twix files
% using the function mapVBVD, which can be found in the repository
% https://github.com/pehses/mapVBVD.

fname       =  { "calib_data_1.mat" ...  % Calibration prior to measurement
               ; "calib_data_2.mat" ...  % Calibration after measurement
               }

Data_Ch     =  4 ;  % Microphone at boc input channel 4

calib_level =  1 ;  % The calibration level is 1 Pa (resp. 94 dB SPL) or 10 Pa (resp. 114 dB SPL)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load Mat files
cdata = [];
for n=1:1:numel(fname)
  load(fname{n})
  cdata = cat(3,cdata,data) ;
end

size_cdata = size(cdata);
n_start    = 1             ;
n_end      = size_cdata(3) ;
adc_len    = size_cdata(1) ;

% Assign the numbers of the physical input channels to the numbering in the twix file
Ch_order  = get_ch_order(data);
Rx_Ch_Ref = Ch_order(1);
Rx_Ch_Sig = Ch_order(Data_Ch);

sigref     = squeeze(cdata(:,Rx_Ch_Ref,n_start:n_end)) ;
sigref     = sigref./abs(sigref); % AM cancellation, normalization to 1
sig        = squeeze(cdata(:,Rx_Ch_Sig,n_start:n_end)) ;


% Flat top window
% Heinzel, G., Rüdiger, A. and Schilling, R. (2002).
% Spectrum and spectral density estimation by the Discrete Fourier transform (DFT),
% including a comprehensive list of window functions and some new flat-top windows.
% Max Planck Institute for Gravitational Physics (Albert Einstein Institute)
z = 2*pi*(0:adc_len-1)/(adc_len-1);
HFT144D    = 1 - 1.96760033*cos(z)   + 1.57983607*cos(2*z) ...
               - 0.81123644*cos(3*z) + 0.22583558*cos(4*z) ...
               - 0.02773848*cos(5*z) + 0.00090360*cos(6*z);
fft_window = HFT144D.';


% Main loop
for n=1:1:n_end
  sig_c    = sig(:,n) .* conj(sigref(:,n)) ;
  sig_c    = sig_c -mean(real(sig_c)) -1i*mean(imag(sig_c));  % remove DC components
  sig_cw   = sig_c.*fft_window ;
  sig_fmag = abs(fftshift(fft(sig_cw)))/adc_len;
  tmp_max  = max(sig_fmag);
  sig_fmax(n) = 2 * tmp_max(1) ;
end

% divide by sqrt(2) to convert peak values ​​into effective (RMS) values
sig_eff = sig_fmax/sqrt(2) ;
cal_vec = calib_level./sig_eff ;
cal_fac = mean(cal_vec) ;

disp('')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Calibration factor for input channel ' num2str(Data_Ch) ': ' num2str(cal_fac) ' Pa/(signal-value [a.u.])'])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('')

% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



