% test with audio toolbox
clear;

fs = 44100;
% framelength = 1024;

dev = audioDeviceReader('Device', 'io|2');
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
