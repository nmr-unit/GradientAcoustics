function Ch_order=get_ch_order(data)
% The channel with the largest amplitude is the carrier channel.
% Empirically, two cases were identified for channel sorting in the TWIX file.
% The carrier signal (physical input channel 1) corresponds to either position 2 or 4.
% Depending on this, the “Ch_order” vector is assigned,
% which contains the positions of all physical input channels.

size_data = size(data) ;
if size_data(2)~=4
  error('Only 4 channels are supported');
end

tmp_ch=mean(abs(squeeze(data(:,:,1))));
Rx_Ch_Ref = find(tmp_ch==max(tmp_ch)) ;
if Rx_Ch_Ref==4
  Ch_order = [4 3 2 1];
elseif Rx_Ch_Ref==2
  Ch_order = [2 1 4 3];
else
  error('No carrier channel found');
end


