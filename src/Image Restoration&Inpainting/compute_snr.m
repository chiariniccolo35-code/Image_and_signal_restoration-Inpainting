function snr = compute_snr(sig,ref)
% Compute Signal-to-Noise Ratio (SNR)
%
% Input:
%       sig         Corrupted signal
%       ref         Reference signal
%  
% Output:
%       snr         SNR value

%Ps  = sum( ( ref(:) - mean(ref(:)) ).^2 ); % signal power (no DC)
Ps  = sum( ( ref(:) ).^2 ); % signal power (no DC)
Pn  = sum( (ref(:) - sig(:)).^2 ); % noise power (assuming no DC)

snr     = 10 * log10( Ps / Pn );
