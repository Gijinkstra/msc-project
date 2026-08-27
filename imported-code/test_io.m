% test with audio toolbox
clear;

fs = 44100;
% framelength = 1024;

dev = audioDeviceReader('Device', "io|2");
set(dev,'NumChannels',2)
setup(dev);

fileWriter = dsp.AudioFileWriter('test.wav','FileFormat','WAV');

disp('Apply Impact Hammer Now (time window = 10 s)')
tic
while toc < 10
    acquiredAudio = dev();
    fileWriter(acquiredAudio);
end
disp('measurement complete.')

release(dev);
release(fileWriter);

y = audioread('test.wav');

figure(1);
clf;
plot(y);

impactHammerSignal = y(:, 2);
celloBodySignal = y(:, 1);

[hammerFft, freqAx] = singleSidedFft(impactHammerSignal, fs);
[celloFft, freqAxC] = singleSidedFft(celloBodySignal, fs);

figure(2)
tiledlayout(2, 1)
t1 = nexttile;
hold(t1, "on")
plot(celloBodySignal, 'LineWidth', 1.1, 'DisplayName', 'Cello body signal')
plot(impactHammerSignal, 'LineWidth', 1.1, 'DisplayName', 'Impact hammer signal')

t2 = nexttile;
hold(t2, "on")
plot(freqAxC, 20*log10(celloFft), 'LineWidth', 1.1, 'DisplayName', 'Cello body signal')
% plot(freqAx, 20*log10(hammerFft), 'LineWidth', 1.1, 'DisplayName', 'Impact hammer signal')

soundsc(celloBodySignal, fs)