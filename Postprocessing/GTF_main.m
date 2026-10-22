clear all;

if exist('OCTAVE_VERSION', 'builtin') > 0
  pkg load signal;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The file ADC_raw_data.mat was created from the corresponding Siemens twix file
% using the function mapVBVD, which can be found in the repository
% https://github.com/pehses/mapVBVD.
MAT_fname = 'ADC_raw_data.mat' ;

% Settings from the sequence
Gradient_Amplitude  = 10     ; % mT/(m/s)
DwellTime           = 80e-6  ; % s

% physical data channel to process (e. g. microphone at box input # 4 )
Data_Ch          =      4 ; % 2, 3, 4 , ...

% Calibration factor in the unit Pa/(signal-level [a.u.])
% corresponds to the output of the script GTF_CALIB.m
cal_fac          =    445.2 ;

% integrate or differentiate the frequency domain signal
% set this flag to +2 if a conversion of acceleration to displacement is
% required (e.g., in case of an acceleration sensor)
integ_fd         =  0    ;  %  0 -> off
                            % +n -> n times integration     in frequency domain
                            % -n -> n times differentiation in frequency domain
% FFT window
Tukey_slope_len  = 12.5e-3 ; % s
rel_wnd_len      =  1      ; % 0.2 .. 1

% A gradient pulse packet consists of 6 repetitions following the sequence
% of gradient directions (+X, +Y, +Z, -X, -Y, -Z).
n_start          =    1 ; % number of the first gradient pulse packet to process
n_pptp           =   50 ; % number of gradient pulse packets to process

phase_corr_man   =    0     ;  %  ==0 -> data driven phase correction
                               %  ~=0 -> manual phase correction (deg)

% The gradient excitation waveform G(t) is stored in 16-bit integer format,
% the sampling rate is always 100 kHz.
fname_G_t = 'Chirp-600ms-4kHz.dat' ;
fs_G_t    = 100e3 ; % Hz

% Plot settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set_fig_nbr =  1  ;

FN_A        = 'sans'; % Font name
FS_A        =    18 ; % Font size

total_plot_width = 360 / 0.352778 ;     % 1 point = 0.3527777778 mm
plot_box_height  =  55 / 0.352778 ;     % 1 point = 0.3527777778 mm

fmin        =   50     ;
fmax        = 4050   ;

plot_lin    =    1    ; % 0 -> linear axes
                        % 1 -> x-axis log frequency, y-axis dB

if plot_lin==1
  ymin=0; ymax=+2.5;
  text_loc   = [fmin,ymax+ 0.55*FS_A*(ymax-ymin)/plot_box_height ];
  ylabel_str = '{\it{|H_{ac,\lambda}(f)|}} [Pa/(mT/m)]';
else
  ymin=-85; ymax=+15;
  text_loc   = [fmin,ymax+ 0.55*FS_A*(ymax-ymin)/plot_box_height ];
  ylabel_str = '{\it{|H_{ac,\lambda}(f)|}} [dB]';
end
xlabel_str = '{\it f } [Hz]' ;

% 'forbidden bands' Skyrafit
ffbc = [ 580 1120 ] ; % Center  [Hz]
ffbr = [ 100  220 ] ; % Span    [Hz]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data file etc. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

load(MAT_fname) ; % array 'data' is complex

size_data = size(data) ;

% Assign the numbers of the physical input channels to the numbering in the twix file
Ch_order  = get_ch_order(data);
Rx_Ch_Ref = Ch_order(1);
Rx_Ch_Sig = Ch_order(Data_Ch);


% Load the gradient excitation waveform
G_t = load(fname_G_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n_start=n_start:1:1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for gd=1:1:3

  n_end = (n_start-1)*6+gd + (n_pptp-1)*6;

  if size_data(3) < n_end
    while n_end>size_data(3)
      n_pptp = n_pptp-1 ;
      n_end = (n_start-1)*6+gd + (n_pptp-1)*6;
    end
  end

  % positive polarity of the gradient excitation
  carrier_1 =  squeeze(data(:,Rx_Ch_Ref,((n_start-1)*6 +gd  ):6:n_end  )) ;
  sig_1     =  squeeze(data(:,Rx_Ch_Sig,((n_start-1)*6 +gd  ):6:n_end  )) ;
  % negative polarity of the gradient excitation
  carrier_2 =  squeeze(data(:,Rx_Ch_Ref,((n_start-1)*6 +gd+3):6:n_end+3)) ;
  sig_2     =  squeeze(data(:,Rx_Ch_Sig,((n_start-1)*6 +gd+3):6:n_end+3)) ;

  carrier   = [carrier_1  carrier_2] ;
  sig       = [sig_1     -sig_2    ] ;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate the Gradient-to-Noise transfer function H_f %%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [H_f(:,gd),F_f,phase_corr]= gtf_pp(G_t,Gradient_Amplitude,sig,DwellTime, ...
                              cal_fac,carrier,Tukey_slope_len,rel_wnd_len, ...
                              phase_corr_man,integ_fd);

end
time_sum = toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the magnitudes of the transfer functions
gtf_plot_mag_XYZ(set_fig_nbr,H_f,F_f,total_plot_width,plot_box_height,plot_lin, ...
fmin,fmax,ymin,ymax,xlabel_str,ylabel_str,text_loc,FN_A,FS_A,ffbc,ffbr)

% Plot the gradient excitation waveform
figure(set_fig_nbr+1)
plot(( 1:1:length(G_t))/fs_G_t,G_t)
grid on
title(fname_G_t)
xlabel('t [s]')
xlim([0 length(G_t)/fs_G_t+0.0125])
ylabel('[a.u.]')
ylim([-2^15-1000 +2^15+1000])

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



