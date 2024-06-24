clear all
close all
clc

set(0, 'DefaultLineLineWidth', 2);


ncfile = '../data/wec_as_multiport_clean.nc';
% ncdisp(ncfile);

Ainf  = squeeze(ncread(ncfile,'added_mass'));
Binf = squeeze(ncread(ncfile,'radiation_damping'));
M = squeeze(ncread(ncfile,'inertia_matrix'));
B = squeeze(ncread(ncfile,'friction'));
g = squeeze(ncread(ncfile,'g'));
K = squeeze(ncread(ncfile,'hydrostatic_stiffness'));

w = squeeze(ncread(ncfile,'omega'));
f = squeeze(ncread(ncfile,'freq'));

Zi = (M + Ainf).*1j.*w + Binf + B + K./(1j.*w);

%% Plot Buoy Impedance
figure
subplot(2,1,1)
semilogx(f,20*log10(abs(Zi)))

subplot(2,1,2)
semilogx(f,angle(Zi))

figure
subplot(2,1,1)
semilogx(f,real(Zi))

subplot(2,1,2)
semilogx(f,imag(Zi))

Temp = squeeze(ncread(ncfile,'excitation_force'));
Hex = Temp(:,1) + 1j*(Temp(:,2));

figure
subplot(2,1,1)
semilogx(f,abs(Hex))

subplot(2,1,2)
semilogx(f,angle(Hex))

%% PTO Param
N = 12.4666;
Kt = 6.1745;
Rw = 0.5;
Lw = 0;
Jd = 2;
Bd = 1;
Kd = 0;

Zd = Jd.*(1j*w) + Bd + Kd./(1j*w);
Zw = Lw.*(1j*w) + Rw;

Z11 = Zd*N^2;
Z12 = Kt*N;
Z21 = Kt*N;
Z22 = -Zw;

%% Controller to maximize mechanical vs electrical power
ZloadMech = Z22 + Z12.*Z21./(conj(Zi)-Z11);
Zout = -Z22 + Z12.*Z21./(Zi + Z11);
Zload = conj(Zout);

% Transducer Gain
Gt_elec = transducerGain(Zi,Z11,Z12,Z21,Z22,Zload);
Gt_mech = transducerGain(Zi,Z11,Z12,Z21,Z22,ZloadMech);

[maglm,phaselm] = bode(frd(ZloadMech,w),w);
[magle,phasele] = bode(frd(Zload,w),w);

figure
subplot(3,1,1)
plot(f,squeeze(maglm),'b--',f,squeeze(magle),'r')
legend('Mechanical','Electrical')
grid on
ylabel('Magnitude')

subplot(3,1,2)
plot(f,squeeze(phaselm),'b--',f,squeeze(phasele),'r')
legend('Mechanical','Electrical','Location','southeast')
grid on
ylabel('Angle [deg]')

subplot(3,1,3)
plot(f,Gt_mech,'b--',f,Gt_elec,'r')
ylim([-1 1])
legend('Mechanical','Electrical')
grid on
ylabel({'Transducer'; 'Power Gain'})
xlabel('Frequency [Hz]')

sgtitle('Zload')

%%
% Zout = -Z22 + Z12.*Z21./(Zi + Z11);
% Ga = availableGain(Zi,Z11,Z21,Zout);
%
% subplot(3,1,3)
% hold on
% area(f,Ga,'FaceAlpha',0.2)
% legend('Mechanical','Electrical','Available')

%% Tune Controller
mycolor = {[216/255 27/255 96/255];...
    [30/255 136/255 229/255];
    [255/255 193/255 7/255]};


Zout = -Z22 + Z12.*Z21./(Zi + Z11);
freq_val = [0.475;0.575;0.675];

figure
subplot(3,1,1)
plot(f,squeeze(magle),'k')
hold on
grid on
ylabel('Magnitude')

subplot(3,1,2)
plot(f,squeeze(phasele),'k')
hold on
grid on
ylabel('Angle [deg]')

subplot(3,1,3)
plot(f,Gt_elec,'k')
hold on
grid on
ylabel({'Transducer'; 'Power Gain'})
xlabel('Frequency [Hz]')

% Go = operatingGain(Z11,Z12,Z21,Z22,Zload);
% subplot(4,1,4)
% plot(f,Go,'k')
% hold on

for ii = 1:length(freq_val)
    freq = freq_val(ii);
    [Kp, Ki] = optimalPI(Zout,Zw,Kt,f,freq);
    disp(Ki)
    Ctrl = Kp + Ki./(1j*w);
    Zload_pi = Kt./Ctrl - Zw;

    [magpi,phasepi] = bode(frd(Zload_pi,w),w);


    Gt = transducerGain(Zi,Z11,Z12,Z21,Z22,Zload_pi);

    subplot(3,1,1)
    plot(f,squeeze(magpi),'--','Color',mycolor{ii})
    grid on
    ylabel('Magnitude')

    subplot(3,1,2)
    plot(f,squeeze(phasepi),'--','Color',mycolor{ii})
    grid on
    ylabel('Angle [deg]')

    subplot(3,1,3)
    area(f,Gt,'FaceAlpha',0.4,'FaceColor',mycolor{ii})
    ylim([0 1])
    grid on
    ylabel({'Transducer'; 'Power Gain'})
    xlabel('Frequency [Hz]')
    hold on

    % Go = operatingGain(Z11,Z12,Z21,Z22,Zload_pi);
    % subplot(4,1,4)
    % area(f,Go,'FaceAlpha',0.2)
    % ylim([0 1])
    % grid on
    % ylabel({'Operating'; 'Power Gain'})
    % xlabel('Frequency [Hz]')
    % hold on
end
sgtitle('Zload')

subplot(3,1,1)
legendStrings = "f_{pi} = " + string(freq_val) + " Hz";
legend(['Z_{\rm out}^*';legendStrings])

%% Kd
mycolor = {[179/255 021/255 041/255];...
    [215/255 095/255 076/255];
    [246/255 164/255 130/255]};

figure
subplot(3,2,1)
plot(f,real(conj(Zi)),'k')
hold on

subplot(3,2,3)
plot(f,imag(conj(Zi)),'k')
hold on

N = 12.4666;
Kt = 6.1745;
Rw = 0.5;
Lw = 0;
Jd = 2;
Bd = 1;
Kd = 0;
Zd = Jd.*(1j*w) + Bd + Kd./(1j*w);
Zw = Lw.*(1j*w) + Rw;
Z11 = Zd*N^2;
Z12 = Kt*N;
Z21 = Kt*N;
Z22 = -Zw;
Zout = CalZout(Z11,Z12,Z21,Z22,Zi);
Gt = transducerGain(Zi,Z11,Z12,Z21,Z22,conj(Zout));
subplot(3,2,[5 6])
plot(f,Gt,'k')
hold on

Kd_val = [0 -25 -50];

for ii = 1:numel(Kd_val)

    Kd = Kd_val(ii);

    Zd = Jd.*(1j*w) + Bd + Kd./(1j*w);
    Zw = Lw.*(1j*w) + Rw;

    Z11 = Zd*N^2;
    Z12 = Kt*N;
    Z21 = Kt*N;
    Z22 = -Zw;

    Zout = CalZout(Z11,Z12,Z21,Z22,Zi);
    Zin = CalZin(Z11,Z12,Z21,Z22,conj(Zout));

    figure(6)
    subplot(3,2,1)
    plot(f,real(Zin),'color',mycolor{ii})
    xlabel('Frequency [Hz]')
    ylabel('Real Zin')
    grid on

    subplot(3,2,3)
    plot(f,imag(Zin),'color',mycolor{ii})
    ylim([-10e3 20e3])
    xlabel('Frequency [Hz]')
    ylabel('Imaginary Zin')
    grid on

    subplot(3,2,2)
    plot(f,real(Zout),'color',mycolor{ii})
    hold on
    xlabel('Frequency [Hz]')
    ylabel('Real Zout')
    grid on

    subplot(3,2,4)
    plot(f,imag(Zout),'color',mycolor{ii})
    hold on
    xlabel('Frequency [Hz]')
    ylabel('Imaginary Zout')
    grid on

    Gt = transducerGain(Zi,Z11,Z12,Z21,Z22,conj(Zout));
    subplot(3,2,[5 6])
    plot(f,Gt,'color',mycolor{ii})
    hold on
    xlabel('Frequency [Hz]')
    ylabel({'Transducer'; 'Power Gain'})
    grid on


    figure(10)
    Gt = transducerGain(Zi,Z11,Z12,Z21,Z22,conj(Zout));
    plot(f,Gt,'color',mycolor{ii})
    ylim([0 1])
    grid on
    ylabel({'Transducer'; 'Power Gain'})
    xlabel('Frequency [Hz]')
    hold on
end

figure(6)
subplot(3,2,[5 6])
legendStrings = "K_{d} = " + string(Kd_val);
legend(['Z_{\rm i}^*',legendStrings],'location','northeastoutside')


%% Jd
mycolor = {[016/255 101/255 171/255];...
    [058/255 147/255 195/255];
    [142/255 196/255 222/255]};

figure(7)
subplot(3,2,1)
plot(f,real(conj(Zi)),'k')
hold on

subplot(3,2,3)
plot(f,imag(conj(Zi)),'k')
hold on

N = 12.4666;
Kt = 6.1745;
Rw = 0.5;
Lw = 0;
Jd = 2;
Bd = 1;
Kd = 0;
Zd = Jd.*(1j*w) + Bd + Kd./(1j*w);
Zw = Lw.*(1j*w) + Rw;
Z11 = Zd*N^2;
Z12 = Kt*N;
Z21 = Kt*N;
Z22 = -Zw;
Zout = CalZout(Z11,Z12,Z21,Z22,Zi);
Gt = transducerGain(Zi,Z11,Z12,Z21,Z22,conj(Zout));
subplot(3,2,[5 6])
plot(f,Gt,'k')
hold on

Jd_val = [2,10,15];

for ii = 1:numel(Jd_val)

    Jd = Jd_val(ii);

    Zd = Jd.*(1j*w) + Bd + Kd./(1j*w);
    Zw = Lw.*(1j*w) + Rw;

    Z11 = Zd*N^2;
    Z12 = Kt*N;
    Z21 = Kt*N;
    Z22 = -Zw;

    Zout = CalZout(Z11,Z12,Z21,Z22,Zi);
    Zin = CalZin(Z11,Z12,Z21,Z22,conj(Zout));

    figure(7)
    subplot(3,2,1)
    plot(f,real(Zin),'color',mycolor{ii})
    ylim([0 5000])
    hold on
    xlabel('Frequency [Hz]')
    ylabel('Real Zin')
    grid on

    subplot(3,2,3)
    plot(f,imag(Zin),'color',mycolor{ii})
    ylim([-10e3 10e3])
    hold on
    xlabel('Frequency [Hz]')
    ylabel('Imaginary Zin')
    grid on

    subplot(3,2,2)
    plot(f,real(Zout),'color',mycolor{ii})
    hold on
    xlabel('Frequency [Hz]')
    ylabel('Real Zout')
    grid on

    subplot(3,2,4)
    plot(f,imag(Zout),'color',mycolor{ii})
    hold on
    xlabel('Frequency [Hz]')
    ylabel('Imaginary Zout')
    grid on


    Gt = transducerGain(Zi,Z11,Z12,Z21,Z22,conj(Zout));
    subplot(3,2,[5 6])
    plot(f,Gt,'color',mycolor{ii})
    hold on
    xlabel('Frequency [Hz]')
    ylabel({'Transducer'; 'Power Gain'})
    grid on


    figure(10)
    Gt = transducerGain(Zi,Z11,Z12,Z21,Z22,conj(Zout));
    plot(f,Gt,'--','color',mycolor{ii})
    ylim([0 1])
    grid on
    ylabel({'Transducer'; 'Power Gain'})
    xlabel('Frequency [Hz]')
    hold on
end

figure(7)
subplot(3,2,[5 6])
legendStrings = "J_{d} = " + string(Jd_val);
legend(['Z_{\rm i}^*',legendStrings],'location','northeastoutside')




% %% Jd with PI
% 
% mycolor = {[016/255 101/255 171/255];...
%     [058/255 147/255 195/255];
%     [142/255 196/255 222/255]};
% 
% 
% N = 12.4666;
% Kt = 6.1745;
% Rw = 0.5;
% Lw = 0;
% Jd = 2;
% Bd = 1;
% Kd = 0;
% Zd = Jd.*(1j*w) + Bd + Kd./(1j*w);
% Zw = Lw.*(1j*w) + Rw;
% Z11 = Zd*N^2;
% Z12 = Kt*N;
% Z21 = Kt*N;
% Z22 = -Zw;
% Zout = CalZout(Z11,Z12,Z21,Z22,Zi);
% Gt = transducerGain(Zi,Z11,Z12,Z21,Z22,conj(Zout));
% subplot(3,2,[5 6])
% plot(f,Gt,'k')
% hold on
% 
% Jd_val = [2,10,15];
% 
% figure
% for ii = 1:numel(Jd_val)
% 
%     Jd = Jd_val(ii);
% 
%     Zd = Jd.*(1j*w) + Bd + Kd./(1j*w);
%     Zw = Lw.*(1j*w) + Rw;
% 
%     Z11 = Zd*N^2;
%     Z12 = Kt*N;
%     Z21 = Kt*N;
%     Z22 = -Zw;
% 
%     Zout = CalZout(Z11,Z12,Z21,Z22,Zi);
% 
%     [Kp, Ki] = optimalPI(Zout,Zw,Kt,f,0.575)
%     Ctrl = Kp + Ki./(1j*w);
%     Zload_pi = Kt./Ctrl - Zw;
% 
%     [magpi,phasepi] = bode(frd(Zload_pi,w),w);
% 
%     Gt = transducerGain(Zi,Z11,Z12,Z21,Z22,Zload_pi);
% 
%     plot(f,Gt)
%     hold on
% 
% end





%%


function Zout = CalZout(Z11,Z12,Z21,Z22,Zi)
Zout = - Z22 + Z12.*Z21./(Zi + Z11);

end

function Zin = CalZin(Z11,Z12,Z21,Z22,Zl)
Zin = Z11 + Z12.*Z21./(Zl - Z22);
end

%%
function G = operatingGain(Z11,Z12,Z21,Z22,Zload)
Zin = CalZin(Z11,Z12,Z21,Z22,Zload);
amp = Z21./(Zload - Z22);
G = amp.*conj(amp) .* real(Zload)./real(Zin);
end
%%
function G = transducerGain(Zi,Z11,Z12,Z21,Z22,Zload)
den = (Zload - Z22).*(Zi + Z11) + Z12.*Z21;
G = 4*Z21.*conj(Z21).*real(Zi).*real(Zload)./(den.*conj(den));
end

%%
function G = availableGain(Zi,Z11,Z21,Zout)
amp = (Z21./(Zi + Z11));
G = amp.*conj(amp).*real(Zi)./real(Zout);
end

%%
function [Kp, Ki] = optimalPI(Zout,Zw,Kt,f,freq)
[~,x] = min(abs(f-freq));
Zctrl_opt = Kt./(conj(Zout) + Zw);
Kp = real(Zctrl_opt(x));
Ki = real(Zctrl_opt(x)*(1j*2*pi*f(x)));
end

%%
