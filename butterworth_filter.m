% clear all, close all, clc
% find data snippet in samples based on electrode 10's neural data
% Create data matrix numElectrodex x numSamplesInSnippet
% Get the video heatmap code section working
%% Load the Neural Data

% Open the neural data file
fileLocation = 'C:\Users\chivu\OneDrive\Documents\MATLAB\research_project_graczyk\NSP Data';
fileName = 'RP01_ICMS_2-Contact_Characterization_2023-01-06-15-06-32-152_NSP1006.ns5';
fileStr = strcat(fileLocation,filesep,fileName);
openNSx(fileStr,'read'); %make sure NPMK folder is added to MATLAB's path
% Open the sensory data file
fileLocationSensory = 'C:\Users\chivu\OneDrive\Documents\MATLAB\research_project_graczyk';
fileNameSensory = 'RP01_ICMS_2-Contact_Characterization_2023-01-06-15-06-32-152.mat';
fileStrSensory = strcat(fileLocationSensory,filesep,fileNameSensory);
load(fileStrSensory);
%% Load neural data
neuralDataUnmapped = (double(NS5.Data)./4).*10^(-6);
% May need to break this up into smaller steps to avoid MATLAB/computer
% overload
%% Organize the time
fs = 30000;
dt = 1/fs;
%TendSamples = length(NS5.Data(1,:));
%TendSec = TendSamples/fs;
TendSec = NS5.MetaTags.DataDurationSec;
TendSamples = NS5.MetaTags.DataPoints;
TstartSec = 0;
TstartSamples = 1;
timeAll = TstartSec:dt:TendSec-dt;
%% Map the neural data from channel-indexing to electrode-indexing
electrodeNumber = zeros(128,1);
zeroChar = double('0');
for i = 1:128
   elecChar =  NS5.ElectrodesInfo(i).Label(7:9);
   elecDouble = double(elecChar);
   elecDouble2 = elecDouble(elecDouble ~= 0);
   numDigits = size(elecDouble2,2);
   if(numDigits == 1)
       elecNumber = elecDouble2-zeroChar;
   elseif(numDigits == 2)
       elecNumber = sum([10 1] .* (elecDouble2-zeroChar));
   elseif(numDigits == 3)
       elecNumber = sum([100 10 1] .* (elecDouble2-zeroChar));
   else
       disp('Error')
   end
   electrodeNumber(i) = elecNumber;
end
elecChannel = [electrodeNumber,[1:128]'];
elecChannelSorted = sortrows(elecChannel,1);
neuralData = neuralDataUnmapped([elecChannelSorted(:,2);129;130],:);
%% Plot the Channel-specific information for visualization
%Stimulated electrode = 10
%Stimulated electrode = 123
%Synch pulse - high during stim, started a little early - channel 130
%Monitor port - stimulation waveform, precise timing - channel 129
for electrode = 130 %1:128
    figure
    plot(timeAll,neuralData(electrode,:));
    xlabel('Time (secs)')
    ylabel('Voltage (uV)')
    titleString = strcat("Voltage Plot for Electrode ",num2str(electrode));
    title(titleString);
    hold on;    
end
synchOutput = neuralData(129,:);
monitorPortOutput = neuralData(130,:);

%% Determine duration of stim artifact
thresh = 0.0000133; % determined from plot of electrode 129
for electrode = 129
    stimDuration = timeAll(neuralData(electrode,:) >= thresh);
end
% use for more precise threshold without outliers
% stimDuration(1,length(stimDuration));
% stimDuration(1,1);
% stimDuration = stimDuration(1,length(stimDuration)) - stimDuration(1,1);

%% Use sample and hold to extract artifact from monitor port
counter = 1;
sample = 0;
holdData = neuralData(130,:);
for counter = 1:length(timeAll)
    if timeAll(1,counter) == stimDuration(1,1)
       if sample == 0
          sample = holdData(1,(counter-1)); % determine sample value
          index = counter; % set the index for first thresh cross
       end
    else
       counter = counter + 1;
    end
end

for counting = 1:length(stimDuration)
    holdData(1,index) = sample;
    index = index + 1;
end

% plot with the sample and hold
for electrode = 130 %1:128
    figure
    plot(timeAll,holdData(1,:));
    xlabel('Time (secs)')
    ylabel('Voltage (uV)')
    titleString = strcat("Voltage Plot for Electrode ",num2str(electrode));
    title(titleString);
    hold on;    
end

%sample = 1.4550e-04
%index = 663886

%% Use monitor port to sample and hold each pulse
clear stimIndex
clear index

threshA = 0.0000133; % determined from plot of electrode 129
threshB = 0.0069129;
count = 1;

% detect all pulses
for electrode = 129
    for i = 1:length(timeAll)
        if neuralData(electrode, i) >= threshA && neuralData(electrode, i) <= threshB
            stimIndex(count) = timeAll(i);
            index(count) = i;
            count = count + 1;
        end
    end
end

% extra detection for pulses that are within neural signal
threshC = 0.0000060; % determined from plot of electrode 129
threshD = 0.0000134;
% do not reset count value - append to end of stimIndex and index vectors
stimExtra = [];
stimExtra(1) = stimIndex(1);
stimExtra(2) = stimIndex(end);

for electrode = 129
    for i = 1:length(timeAll)
        if stimExtra(1) <= timeAll(i) && stimExtra(end) >= timeAll(i)
            if neuralData(electrode, i) >= threshC && neuralData(electrode, i) <= threshD
                stimIndex(count) = timeAll(i);
                index(count) = i;
                count = count + 1;
            end
        end
    end
end


% determine the sampling value
count = 1;
sample = 0;
electrodeNum = 130;

holdData = neuralData(electrodeNum,:);
for counter = 1:length(timeAll)
    if timeAll(1,counter) == stimIndex(1,count)
        if sample == 0
          sample = holdData(1,(counter-1)); % determine sample value
        end
    end
end

% hold for 700 to 800 us to cut out the pulse
for counting = 1:length(stimIndex) % hold for 700 - 800 us (~ 0.0007)
    clear count;
    count = index(counting);
    for i = 1:30
        holdData(1,count) = sample;
        count = count + 1;
     end
end


% plot with the sample and hold
for electrode = electrodeNum %1:128
    figure
    plot(timeAll,holdData(1,:));
    xlabel('Time (secs)')
    ylabel('Voltage (uV)')
    titleString = strcat("Voltage Plot for Electrode ",num2str(electrode));
    title(titleString);
    hold on;    
end

%sample = 1.4550e-04
%index = 663886

% plotted with butterworth filter
F = 250;
Fs = 1000;
[b, a] = butter(4, F/Fs, 'high');
inputSignal = holdData(1, :);
outSignal = filter(b, a, inputSignal);

plot(timeAll, holdData(1, :))
hold on
plot(timeAll, outSignal)
legend('blue line', 'raw data', 'filtered data')
hold off


%% Use butterworth filter

% electrode 10
F = 250;
Fs = 1000;
[b, a] = butter(4, F/Fs, 'high');
inputSignal = neuralData(10, :);
outSignal = filter(b, a, inputSignal);
figure
plot(timeAll, outSignal)

plot(timeAll, neuralData(10, :))
hold on
plot(timeAll, outSignal)
hold off

figure
plot(timeAll, neuralData(10, :))


%% butterworth filter with sample and hold
counter = 1;
sample = 0;
holdData = neuralData(15,:);
for counter = 1:length(timeAll)
    if timeAll(1,counter) == stimDuration(1,1)
       if sample == 0
          sample = holdData(1,(counter-1)); % determine sample value
          index = counter; % set the index for first thresh cross
       end
    else
       counter = counter + 1;
    end
end

for counting = 1:length(stimDuration)
    holdData(1,index) = sample;
    index = index + 1;
end

% plot with the sample and hold
for electrode = 15 %1:128
    figure
    plot(timeAll,holdData(1,:));
    xlabel('Time (secs)')
    ylabel('Voltage (uV)')
    titleString = strcat("Voltage Plot for Electrode ",num2str(electrode));
    title(titleString);
    hold on;    
end

F = 250;
Fs = 1000;
[b, a] = butter(4, F/Fs, 'high');
inputSignal = holdData(1, :);
outSignal = filter(b, a, inputSignal);

plot(timeAll, holdData(1, :))
hold on
plot(timeAll, outSignal)
hold off


%%
timeAll(1, 3:24)

%% Clean up
close all;