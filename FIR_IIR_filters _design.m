%% Q1
% for Impulse Invariance Method
N=8;
omega_c_II=0.859;
[b,a]=butter(N,omega_c_II,'s'); 

% for Bilinear Method
N=9;
omega_c_BL=0.926;
[c,d]=butter(N,omega_c_BL,'s'); 

[H_II_BW,w]=freqs(b,a,1024);
[H_BL_BW,w]=freqs(c,d,1024);

figure
plot(w,abs(H_II_BW))
hold on
plot(w,abs(H_BL_BW))
legend('Impulse Invariance','Bilinear Transform')

%We find the Impulse Invariance and Bilinear filters now 
[b_z,a_z]= impinvar(b,a,1);
[c_z,d_z]= bilinear(c,d,1);

%We plot the IIR Filters
[H_II,w]=freqz(b_z,a_z,1024);
[H_BL,w]=freqz(c_z,d_z,1024);
    
figure
plot(w,abs(H_II))
hold on
plot(w,abs(H_BL))
legend('Impulse Invariance','Bilinear Transform')
%% Q2
M=22;
w_c=0.3;
w_square=rectwin(M);
w_hamming=hamming(M);
w_hann=hann(M);
w_blackman=blackman(M);
w_bartlet=bartlett(M);

b_square=fir1(M-1,w_c,w_square);
b_hamming=fir1(M-1,w_c,w_hamming);
b_hann=fir1(M-1,w_c,w_hann);
b_blackman=fir1(M-1,w_c,w_blackman);
b_bartlet=fir1(M-1,w_c,w_bartlet);

%we create an ideal LPF by taking fir of order 1024
figure 
subplot(2,1,1);
Ideal_LPF=[rectwin(307)' zeros(1,length(w)-307)];
plot(w,Ideal_LPF)
subplot(2,1,2); 
plot(w,angle(Ideal_LPF))


figure
[T_1,w]=freqz(b_square,1,1024);
plot(w,abs(T_1))
title('square')
pks = findpeaks(abs(T_1));
X = ['Peak side-lobe amplitudev of square:  ',num2str(pks(2)),'  MSE of square:  ',num2str(immse(T_1,Ideal_LPF'))];
disp(X)
figure
freqz(b_square,1,1024);
title('square')

figure
[T_2,w]=freqz(b_hamming,1,1024);
plot(w,abs(T_2))
title('hamming')
pks = findpeaks(abs(T_2));
X = ['Peak side-lobe amplitudev of hamming:  ',num2str(pks(2)),'  MSE of hamming:  ',num2str(immse(T_2,Ideal_LPF'))];
disp(X)
figure
freqz(b_hamming,1,1024);
title('hamming')

figure
[T_3,w]=freqz(b_hann,1,1024);
plot(w,abs(T_3))
title('hann')
pks = findpeaks(abs(T_3));
X = ['Peak side-lobe amplitudev of hann:  ',num2str(pks(2)) ,'  MSE of hann:  ',num2str(immse(T_3,Ideal_LPF'))];
disp(X)
figure
freqz(b_hann,1,1024);
title('hann')

figure
[T_4,w]=freqz(b_blackman,1,1024);
plot(w,abs(T_4))
title('blackman')
pks = findpeaks(abs(T_4));
X = ['Peak side-lobe amplitudev of blackman:  ',num2str(pks(2)),'  MSE of blackman:  ',num2str(immse(T_4,Ideal_LPF'))];
disp(X)
figure
freqz(b_blackman,1,1024);
title('blackman')

figure
[T_5,w]=freqz(b_bartlet,1,1024);
plot(w,abs(T_5))
title('bartlet')
pks = findpeaks(abs(T_5));
X = ['Peak side-lobe amplitudev of bartlet:  ',num2str(pks(2)),'  MSE of bartlet:  ',num2str(immse(T_5,Ideal_LPF'))];
disp(X)
figure
freqz(b_bartlet,1,1024);
title('bartlet')

%% Q3
n=1:29;
M=29;
w_c=0.38*pi;
w=kaiser(29,3.549542) ;
figure
stem(w)
xlabel('n', 'FontSize', 11);
ylabel('w[n]', 'FontSize', 11);
title('Kiesser window')

h_n=sin(w_c*(n-M/2))./(pi*(n-M/2));
figure
stem(h_n)

figure
[H,w] = freqz(h_n,1);
plot(w,abs(H))
xlabel('w', 'FontSize', 11);
ylabel('H(w)', 'FontSize', 11);
figure
freqz(h_n,1)


Ideal_LPF=[rectwin(196)' zeros(1,length(w)-196)];
plot(w,Ideal_LPF)
hold on
plot(w,abs(H))
error_aprrox=abs(H)-Ideal_LPF';
hold off
figure
plot(w,20*log10((error_aprrox)))
xlabel('w', 'FontSize', 11);
ylabel('Error Aprox.(w)', 'FontSize', 11);