function f = Process(app)

    f0 = 1.5e6; % ultrasound frequency, Hz
    
    %limits for frequency filtering, Hz
    lowf = 2.5e6; 
    lowlim = ceil(lowf / f0);
    highf = 6.0e6; 
    highlim = floor(highf / f0);
    
    BW = 100e3; %width of the notch for the notch filter, Hz
    NumPoint = app.Points.Value; %number of treatment points with the same acoustic parameters
    NumPulse = app.Pulses.Value; %number of pulses delivered per point (usually 60)

    CavitationPersistence = zeros(1, NumPoint);
    CavitationMap = zeros(NumPulse, NumPoint); %binary indicator (1=yes, 0=no) of whether cavitation had occured within a certain pulse and point
    PCDNoiseMap = zeros(NumPulse, NumPoint); %broadband noise value for each recorded signal
    ICDoseMEAN = zeros(NumPulse, NumPoint);

    Start = app.Start.Value; %number assigned by the scope to the first signal
    
    %main loop through each layer 
    for nSPOT = 1 : NumPoint

        %set temporary arrays for each nSPOT
        IfCavitation = zeros(1, NumPulse); %store binary indicator of cavitation occurence for each of 60 pulses at this focal point
        PCDNoise = zeros(1, NumPulse); %store broadband noise for each of 60 pulses at this focal point

        %iterating through each pulse
        for nPulse = 1 : NumPulse

            if (Start+(nSPOT-1)*NumPulse+nPulse-1)<10
                Name=[app.directory, '/C1Trace0000',num2str(Start+(nSPOT-1)*NumPulse+nPulse-1),'.trc.txt']
            elseif (Start+(nSPOT-1)*NumPulse+nPulse-1)<100
                Name=[app.directory, '/C1Trace000',num2str(Start+(nSPOT-1)*NumPulse+nPulse-1),'.trc.txt']
            elseif (Start+(nSPOT-1)*NumPulse+nPulse-1)<1000
                Name=[app.directory, '/C1Trace00',num2str(Start+(nSPOT-1)*NumPulse+nPulse-1),'.trc.txt']
            else 
                Name=[app.directory, '/C1Trace0',num2str(Start+(nSPOT-1)*NumPulse+nPulse-1),'.trc.txt']
            end

            MyDataIm = importdata(Name);
            t = MyDataIm(:,1); %the time vectore
            Ft = MyDataIm(:,2); %the signal vector

            %setting up filters to filter the imported signal
            Fs = 1e5 / (t(2) * 1e5 - t(1) * 1e5); %frequency domain range

            C = fir1(1000, 2 * [lowf highf] / Fs); %band pass filter
            FilterX = filtfilt(C, 1, Ft); %execute the bandpass filtering

            for nn = lowlim : 0.5 : highlim
                wo = nn * f0 / (Fs / 2);  bw = BW / (Fs / 2);
                [b,a] = iirnotch(wo, bw);
                FilterX = filtfilt(b, a, FilterX); %execute the notch filtering for every notch
            end
            
            %plot time vs. signal
            l = line(app.FilterX, t, FilterX);

            %TODO: FilterX is the filtered PCD signal; in GUI, we want to plot FilterX vs t for each signal.  
            %to get the frequency spectrum of the signal, you'd need to plot
            %abs(fft(FilterX)) - again, for each signal

            %Rose criterion
            RangeNoiseStart = find(t > 10e-6,1); %start of the electric noise portion of the signal, seconds
            RangeNoiseEnd = find(t > 55e-6,1); %end of the electric noise portion of the signal, seconds   
            RangeFocusStart = find(t > 60e-6,1); %start of the PCD signal from focus
            RangeFocusEnd = find(t > 9.6e-4,1); %end of the PCD signal from focus

            FilterXNoise=abs(FilterX(RangeNoiseStart:RangeNoiseEnd)); %defining what is electric noise and what is its maximum amplitude
            RefThreshold=max(FilterXNoise);

            if max(abs(FilterX(RangeFocusStart:RangeFocusEnd)))>=(2.23*RefThreshold) %comparing the "real" signal to the reference noise amplitude - the actual Rose criterion
                IfCavitation(nPulse) = 1; %binary outcome: 1=yes, cavitation occured, 0=no, cavitation did not occur
            else
                IfCavitation(nPulse) = 0;
            end

            %if cavitation occured, we calculate the cumulative noise
            %amplitude for that pulse, PCDNoise
            if IfCavitation(nPulse) == 1
                FilterXROI = FilterX(RangeFocusStart : RangeFocusEnd);
                PCDNoise(nPulse) = sum(abs(FilterXROI)) / length(FilterXROI);
            else
                IfCavitation(nPulse) = 0;
                PCDNoise(nPulse) = 0;
            end

        end

        CavitationMap(:, nSPOT) = IfCavitation; %adding the binary cavitation readouts for all 60 pulses at this focal point to the overall matrix
        PCDNoiseMap(:, nSPOT) = PCDNoise; %same for cavitation noisenoise

        %calculate mean cavitation noise and persistence for this point
        %TODO we want to show the values of ICDoseMEAN and
        %CavitationPersistence updated for each spot in the indicator
        %windows of the GUI
        PCDNoiseRemoveZero = PCDNoise(find(PCDNoise~=0));
        if isempty(PCDNoiseRemoveZero) == 1
            PCDNoiseRemoveZero = 0;
        end
        ICDoseMEAN(nSPOT) = mean(PCDNoiseRemoveZero); 
        CavitationPersistence(nSPOT) = sum(IfCavitation) / length(IfCavitation) * 100; 

    end

    %overall mean values for persistence, noise and probability that we want to
    %display at the end when the script ends running
    Persist = mean(CavitationPersistence) %mean persistence value
    %STDPersist=std(CavitationPersistence)/sqrt(length(CavitationPersistence))%standard deviation of persistence value (no need to display it just yet)
    MeanNoise = mean(PCDNoiseMap(find(PCDNoiseMap~=0))) %mean broadband noise
    %STDNoise=std(PCDNoiseMap(find(PCDNoiseMap~=0)))%/sqrt(length(ICDoseMEAN(find(ICDoseMEAN~=0))))
    Probability = length(CavitationPersistence(find(CavitationPersistence~=0)))*100/NumPoint %cavitation probability in this experiment

    f = [Persist, MeanNoise, Probability];

end

