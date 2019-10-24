clear all
close all


t = linspace(0,4,1e4);


yy = sin(2*pi*t) + 2*sin(2*pi*2*t) + 3*sin(2*pi*3*t) + 0.5*sin(2*pi*20*t);

figure(1), clf;

hold on
plot(t,yy)




ofyy = fft(yy);

%%
fyy = ofyy(2:end/2);

[~,maxmode] = max(abs(fyy));
freq = (maxmode) / numel(yy) / mean(diff(t));
sprintf('%.8f',freq)
phi = unwrap(angle(fyy(maxmode)))

plot(t(1:numel(ofyy)),ifft(ofyy));
%%
% only keep good modes
fyy_ = ofyy;
[a,b] = sort(fyy_);
good_modes = false(size(fyy_));
good_modes(b(end-5:end)) = true;

fyy_(~good_modes) = 0;

plot(t,ifft(fyy_));

%% only keep first 10 modes

sfyy = ofyy;

sfyy(15:end-15) = 0;
plot(t,ifft(sfyy));
plot(t,abs(ifft(sfyy)) .* sign(real(ifft(sfyy))) );

%% only keep a few modes

clear all

t = linspace(0,4,1e4);

yy = sin(2*pi*t) + 2*sin(2*pi*2*t) + 3*cos(2*pi*3*t) + 0.5 * sin(2*pi*5*t);

ofyy = fft(yy,numel(yy));

modes_to_keep = 10;
fyy = ofyy;
fyy(modes_to_keep+1:end-(modes_to_keep-1)) = 0;
% fyy(1:50)

figure(1),clf
hold on
plot(t,yy)
plot(t,ifft(ofyy))

plot(t, 2./numel(yy) .* (-imag(ofyy(5)) * sin(2*pi*t) - imag(ofyy(9)) * sin(2*pi*2*t) - imag(ofyy(13)) * sin(2*pi*3*t) + ...
    real(ofyy(5)) * cos(2*pi*t) + real(ofyy(9)) * cos(2*pi*2*t) + real(ofyy(13)) * cos(2*pi*3*t))) 



figure(2); clf 
hold on
plot(abs(ofyy(1:50)))
plot(abs(fyy(1:50)))



%% stay on 1 period


clear all

t = linspace(0,1,1e3);
t(end-10:end) = [];

yy = sin(2*pi*t) + 2*sin(2*pi*4*t) + 3*cos(2*pi*8*t) + 0.5 * sin(2*pi*10*t);

ofyy = fft(yy,numel(yy));

modes_to_keep = 10;
fyy = ofyy;
fyy(modes_to_keep+1:end-(modes_to_keep-1)) = 0;
% fyy(1:50)



figure(1),clf
hold on
plot(t,yy)
plot(t,ifft(ofyy))

% dummy_mat = zeros(modes_to_keep,numel(yy));
T = max(t);
pluto = arrayfun(...
    @(m,tt)  real(ofyy(m)) * cos(2*pi*(m-1)*tt/T) - imag(ofyy(m)) * sin(2*pi*(m-1)*tt/T) ,...
    repmat((1:modes_to_keep)',1,numel(t)),...
    repmat(t,modes_to_keep,1) );

pippo = 2./numel(yy) .* sum(pluto);

plot(t, pippo)
% plot(t, 2./numel(yy) .* arrayfun(@(i) sum( real(ofyy(i))*cos(2*pi*i*t)  - imag(ofyy(2)) * sin(2*pi*i*t) ),[2 5 9])) 



figure(2); clf 
hold on
plot(abs(ofyy(1:50)))
plot(abs(fyy(1:50)))