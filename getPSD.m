function [psd,freq,dfreq,psdsg] = getPSD(data,dt,gettrend,order,width,datum)
    % getPSD()
    % input
    % data : NxM data array
    % dt   : sampling size along N
    % gettrend : true/false of running trend sgolay analysis
    % order : trend analysis order (used when gettrend = true)
    % width : trend analysis window (used when gettrend = true)
    % datum : AxB data array include 1 sampling time column and M datum
    % data columns. Leave it [] if you want to run without.
    % 
    % output
    % psd : power spectrum density
    % freq : corresponding frequencies
    % dfreq : frequency sampling size
    % psdsg : psd trend
    
    [N,M] = size(data);
    L = round(N/2,0);
    freq = linspace(0,0.5/dt,L)';
    dfreq = mean(diff(freq));
    psd = zeros(size(data));
    for i=1:M
        psd(:,i) = (abs(fft(data(:,i)))/L).^2;
    end
    psd = 2*psd(1:L,:);
    
    % Calibrate out datum
    if ~isempty(datum)
        assert(size(datum,2)-1 == M);
        psd_datum = interp1(datum(:,1),datum(:,2:end),freq,'cubic');
        psd = psd - psd_datum;
    end    
    
    psdsg = zeros(size(psd));
    if gettrend
        for i=1:M
            psdsg(:,i) = sgolayfilt(psd(:,i),order,width);
        end
    end
end

