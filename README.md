# Design-of-FIR-Filters-using-hamming-window

# DESIGN OF LOW PASS FIR DIGITAL FILTER 
# VENKATESAN S
# 212223060296
# AIM: 
          
  To generate design of high pass FIR digital filter using SCILAB 

# APPARATUS REQUIRED: 

  PC Installed with SCILAB 

# PROGRAM 
## Low Pass FIR Filter
```
clc;
clear;
close;

// Low Pass FIR Filter using Hamming Window
N = 31;              // Filter length
fc = 0.3;            // Normalized cutoff frequency (fc = Fcutoff / (Fs/2))
n = 0:N-1;
alpha = (N-1)/2;

// Hamming window
w = 0.54 - 0.46*cos(2*%pi*n/(N-1));

// Ideal Low Pass Filter Impulse Response
hd = zeros(1, N);
for i = 1:N
    if (i-1) == alpha then
        hd(i) = 2*fc;
    else
        hd(i) = sin(2*%pi*fc*(i-1-alpha)) / (%pi*(i-1-alpha));
    end
end

// Multiply by window
h = hd .* w;

// Frequency Response
[H, f] = frmag(h, 512);

// Plot
figure;
subplot(2,1,1);
plot(f, 20*log10(abs(H)));
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
title('LOW PASS FIR FILTER (Hamming Window)');

subplot(2,1,2);
plot(f, atan(imag(H), real(H)));
xlabel('Normalized Frequency');
ylabel('Phase (radians)');
title('Phase Response');
```
## High Pass FIR Filter
```
clc;
clear;
close;

// High Pass FIR Filter using Hamming Window
N = 31;              // Filter length
fc = 0.3;            // Normalized cutoff frequency
n = 0:N-1;
alpha = (N-1)/2;

// Hamming window
w = 0.54 - 0.46*cos(2*%pi*n/(N-1));

// Ideal High Pass Filter Impulse Response
hd = zeros(1, N);
for i = 1:N
    if (i-1) == alpha then
        hd(i) = 1 - 2*fc;
    else
        hd(i) = -sin(2*%pi*fc*(i-1-alpha)) / (%pi*(i-1-alpha));
    end
end

// Multiply by window
h = hd .* w;

// Frequency Response
[H, f] = frmag(h, 512);

// Plot
figure;
subplot(2,1,1);
plot(f, 20*log10(abs(H)));
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
title('HIGH PASS FIR FILTER (Hamming Window)');

subplot(2,1,2);
plot(f, atan(imag(H), real(H)));
xlabel('Normalized Frequency');
ylabel('Phase (radians)');
title('Phase Response');
```

# OUTPUT
## Low Pass FIR Filter
<img width="1920" height="1080" alt="Screenshot 2025-11-18 112528" src="https://github.com/user-attachments/assets/b3598038-32fa-4ba1-88d4-d5341153fa64" />

## High Pass FIR Filter

<img width="1920" height="1080" alt="Screenshot 2025-11-18 112603" src="https://github.com/user-attachments/assets/9ee89e06-eb6a-4588-bb71-f676b3275c4a" />

# RESULT
Design of high pass FIR digital filter using SCILAB is generated
