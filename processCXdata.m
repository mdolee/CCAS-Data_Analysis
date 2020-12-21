function [t,dt,a,v,l,lsg] = processCXdata(CXfile,removetrend,order,tdur)
    %PROCESSCXDATA Summary of this function goes here
    %   Detailed explanation goes here

    [t, a] = readCXfile(CXfile); % Read CX file
    dt = median(diff(t)); % Get the sampling time 
    g = 9.81; % Gravitational constatnt m/s^2

    m = ~isnan(a(:,1)); % Find NaN value in acc
    a = a(m,:) .* g; % Scale acc from g to m/s^2
    t = t(m); 
    t = t - t(1);
    a = a - repmat(mean(a,1),size(a,1),1);
    
    % Estimate velocity by integration of acceleration
    v = cumtrapz(t,a);
    v = v - repmat(mean(v,1),size(v,1),1);
    
    % Estimate displacement by integration of velocity
    l = cumtrapz(t,v);

    % Savitzky-Golay filtering to remove large scale trend which is likely
    % intrinsic to the accelerometer rather than anything physical.
    % order 3 and filter filter duration 1 sec. (Signal processing toolbox
    % needed).
    wind= round(tdur/dt,0);
    if mod(wind,2) == 0
        wind = wind + 1;
    end

    lsg = zeros(size(l));    
    % HL update 11/11/20
    if removetrend
        for i=1:size(l,2)
            lsg(:,i) = sgolayfilt(l(:,i),order,wind);
        end
        % Calibrate out the large scale trend 
        l = l - lsg;
    end
end

