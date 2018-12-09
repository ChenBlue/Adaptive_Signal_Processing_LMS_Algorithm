N = 2000; % length of sequence
M = 11; %taps
sim_times = 200; % simulation times
tau = 0; % delays
mu = 0.075; % step size

W = 2.9; %Channel parameter
n = [1 2 3];
h_tmp = 1/2*(1+cos(2*pi*(n-2)/W)); 
h = [0 h_tmp]; % raised cosine channel

U = zeros(1,M); 
mse = zeros(1,N); % mean squared error
y = zeros(1,N); % output of equalizer
e = zeros(1,N); % error
%u_rev = fliplr(u);
for j = 1:sim_times
    d = double(randn(1,N)>0); % input signal
    u_tmp = conv(d,h);
    u = [zeros(1,M-1) u_tmp(1:N)];
    w = zeros(1,M);

    for i=tau+1:N
        U = fliplr(u(i:i+M-1));
        y(i) = U * transpose(w);
        e(i) = d(i-tau)-y(i);
        w = w + mu*e(i)*U;
    end
    mse = mse + e.^2;
end

mse = mse/sim_times;
figure,plot(mse);
set(gca, 'YScale', 'log')
xlabel('Number of adaptation cyckes, n');
ylabel('Mean squarred error');
title('Learning curve');
figure,plot(w);
