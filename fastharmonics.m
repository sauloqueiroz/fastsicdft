%{
/*
 * Copyright 2024 Saulo Jorge Beltrao de Queiroz
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
%}
%pkg load signal %for Octave only

% SIC DFT algorithm
function xhat = compactsic(x)
  sqrtN = sqrt(length(x)); 
  xhat = zeros(1, sqrtN); % Initialize output DFT vector
  for n = 0:sqrtN-1
    xhat(n+1) = 0;
    for r = 0:sqrtN-1
      xhat(n+1) = xhat(n+1) + x(n + 1 + r*sqrtN);
    end
  end
end


% Use for the computation of a specific frequency
function X_k = compute_dft_frequency(x, N, k)
  X_k = 0;                
  for n = 0:N-1
    X_k = X_k + x(n+1) * exp(-2j * pi * k * n / N);
  end
end

[x, Fs] = audioread('piano_A4_Fs38720Hz.wav'); %fs 7744

%=================================================================>
Nsic = 7744;    % number of points to analyze
Nfft=2^13;
sqrtN=sqrt(Nsic);


%=================================================================> FFT Analysis
%----------------------------> signal windwing to combat spectral leakage and FFT
alpha = 0.6; %for tukey windowing
windowfft = tukeywin(Nfft, alpha); 
x_windowedfft = x(1:Nfft) .* windowfft';
tic
yfft=fft(x_windowedfft);
toc
%=================================================================> Compacted SIC DFT Analysis
%----------------------------> signal windwing to combat spectral leakage and DFT compacting
xhat=compactsic(x(1:Nsic)); 
windowsic = tukeywin(sqrtN, alpha);  % tukey window
xhatwindowed = xhat(1:sqrtN) .* windowsic'; %windowing the signal to be transformed

Y=zeros(1,Nsic);
Y(Y == 0) = 10e-5;%noise floor. zeroes cause problem in semilogy plotting

% This block computes only the fundamental and and the H harmonics. Thus, the
% number of complex multiplications (H+1)*sqrtN is far lower than the 2^13-point FFT execution.
% Note, however, that the code is not optimimized to hardware as the native FFT implementation is.
% Future work can do this for fairer performance comparisons. 
Y(1)=compute_dft_frequency(xhatwindowed, sqrtN, 0); %DC frequency
nHarmonics=6;  %1 fundamental frequency + number of harmonics (5, in this case).
%for k = 1:sqrtN-1 %we dont need all sqrtN freqs. Uncomment this otherwise.
for k = 1:nHarmonics %includes fundamental frequency i.e. 1st iteration computes fund. freq.
    Y(k*sqrtN + 1) = compute_dft_frequency(xhatwindowed, sqrtN, k); 
    %Y(sicIndexes(k+1) + 1) = compute_dft_frequency(xhatwindowed, sqrtN, k); % avoid index multiplication as above
end


%uncomment this to run FFT on the signal compacted by our proposed algorithm (SIC)
%{
next_power_of_two = 2^ceil(log2(sqrtN));
pad=zeros(next_power_of_two - sqrtN, 1);
xhatwindowed_padded = [xhatwindowed, pad'];
Xhat=fft(xhatwindowed_padded);

Y(1) = Xhat(1);
for k = 2:sqrtN-1
  Y((k-1)*sqrtN + 1) = Xhat(k);
end
%}

%=================================================================> Plotting
deltaFfft = Fs / Nfft;
deltaFsic = Fs / Nsic;

frequenciessic = (0:Nsic/2-1) * deltaFsic;  % SIC Frequency vector
frequenciesfft = (0:Nfft/2-1) * deltaFfft;  % FFT Frequency vector
magnitudefft = abs(yfft(1:Nfft/2)) / Nfft;  % Magnitude of FFT
magnitudesic = abs(Y(1:Nsic/2)) / Nsic;     % Magnitude of SIC

figure;

semilogy(frequenciessic(2:Nsic/2), magnitudesic(2:Nsic/2), 'LineWidth', 3);  
hold on
semilogy(frequenciesfft(1:Nfft/2), magnitudefft(1:Nfft/2));  
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%title('Spectral Analysis of A440 Key of Piano');

l=legend('88-point SIC DFT', '8192-point FFT')
set(l, 'FontSize', 14);  % Change the font size to 14 (adjust as needed)
set(gca, 'FontSize', 14);  % Adjust the font size of the axis labels and tick labels

grid on;
axis([0 3000 10^-5 1])
xticks(0:440:3000)
print('A4pianokey.eps', '-depsc');
