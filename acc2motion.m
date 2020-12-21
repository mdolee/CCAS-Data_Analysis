format long g
scriptFolder = pwd;
dataFolder = 'CCAS_vibration_data_sample/';

%% Process data
tdur = 5.0; order = 3; width = 71; arcsec = 0.000126/2.355;

tic;
dataFile = fullfile(scriptFolder,dataFolder,'CX1_1027 - CX1_1027 - A - 20190726_210802.csv');
[t0,dt0,a0,v0,l0,lsg0] = processCXdata(dataFile,true,order,tdur);
[psd0,freq0,dfreq0,psdsg0] = getPSD(l0,dt0,true,order,width,[]);

dataFile = fullfile(scriptFolder,dataFolder,'Outer_Tower - CX1_2613 - AC - 20201101_015058.csv');
[to,dto,ao,vo,lo,lsgo] = processCXdata(dataFile,true,order,tdur);
[psdo,freqo,dfreqo,psdsgo] = getPSD(lo,dto,true,order,width,[freq0,psdsg0]);

dataFile = fullfile(scriptFolder,dataFolder,'Inner_Tower - CX1_2612 - AC - 20201101_015108.csv');
[ti,dti,ai,vi,li,lsgi] = processCXdata(dataFile,true,order,tdur);
[psdi,freqi,dfreqi,psdsgi] = getPSD(li,dti,true,order,width,[freq0,psdsg0]);
toc;

%% Display results
fig = figure(1);
fig.Position = [150,150,1900,800];
fig.Color = [1 1 1];

% Display the datum power spectrum
ax = subplot(1,2,1);
hold on;
ph = plot(ax,freq0,psd0); 
ph(1).LineStyle = 'none'; ph(1).Marker='o'; ph(1).MarkerFaceColor='w';ph(1).LineWidth=1;ph(1).MarkerSize=3;
ph(2).LineStyle = 'none'; ph(2).Marker='s'; ph(2).MarkerFaceColor='w';ph(2).LineWidth=1;ph(2).MarkerSize=3;
ph(3).LineStyle = 'none'; ph(3).Marker='d'; ph(3).MarkerFaceColor='w';ph(3).LineWidth=1;ph(3).MarkerSize=3;
ph = plot(ax,freq0,psdsg0);
ph(1).LineWidth = 2; ph(2).LineWidth = 2; ph(3).LineWidth = 2;
lh = yline(arcsec^2); lh.LineWidth=2;
hold off;
ax.XScale = 'log';
ax.YScale = 'log';
ax.XLim=[0.1,30];
ax.YLim=[1e-19,1e-6];
ax.XTick=[0.1,0.2,0.3,0.4,0.5,1,2,3,4,5,10,20,30];
ax.XTickLabel = cellstr(string(round(ax.XTick,1))');
ax.GridLineStyle = ':';
legend(["x","y","z","x(filt)","y(filt)","z(filt)","1'' fwhm"]);
th = title(ax,"(A) Pier data (as calibration)",'Units','Normalized', ...
     'Position',[0.17 0.95 0],'FontName','Calibri','FontSize',16);
xlabel(ax,"Frequency [Hz]",'FontName','Calibri','FontSize',12);
ylabel(ax,"PSD [m^2/Hz]",'FontName','Calibri','FontSize',12);
grid on; box on;
pos = ax.Position;
pos(1)=0.05; pos(3) = 0.4; ax.Position = pos;
ax2 = axes('Position',pos,'XAxisLocation','top','YAxisLocation','left','Color','none');
ax2.XScale = ax.XScale;
ax2.YScale = ax.YScale;
ax2.XLim=ax.XLim;
ax2.YLim=ax.YLim;
ax2.YTick=ax.YTick;
ax2.XTick=ax.XTick;
ax2.XTickLabel = cellstr(string(round(1./ax2.XTick/2.*1000,0))');
ax2.XAxis.FontSize = 8;
xlabel(ax2,'Sampling Time [ms]','FontName','Calibri','FontSize',12);

% Display the outer tower power spectrum
ax = subplot(2,2,2);
hold on;
ph = plot(freqo,psdo);
ph(1).LineStyle = 'none'; ph(1).Marker='o'; ph(1).MarkerFaceColor='w';ph(1).LineWidth=1;ph(1).MarkerSize=3;
ph(2).LineStyle = 'none'; ph(2).Marker='s'; ph(2).MarkerFaceColor='w';ph(2).LineWidth=1;ph(2).MarkerSize=3;
ph(3).LineStyle = 'none'; ph(3).Marker='d'; ph(3).MarkerFaceColor='w';ph(3).LineWidth=1;ph(3).MarkerSize=3;
ph = plot(freqo,psdsgo);
ph(1).LineWidth = 2; ph(2).LineWidth = 2; ph(3).LineWidth = 2;
lh = yline(arcsec^2); lh.LineWidth=2;
hold off;
ax.XScale = 'log';
ax.YScale = 'log';
ax.XLim=[0.1,30];
ax.YLim=[1e-15,1e-6];
ax.XTick=[0.1,0.2,0.3,0.4,0.5,1,2,3,4,5,10,20,30];
ax.XTickLabel = cellstr(string(round(ax.XTick,1))');
% ax.XTickLabel = cellstr(string(round(1./ax.XTick/2.*1000,0))');
ax.GridLineStyle = ':';
legend(["x","y","z","x(filt)","y(filt)","z(filt)","1'' fwhm"]);
title(ax,"(B) Outer tower (calibrated)",'Units','Normalized', ...
     'Position',[0.15 0.9 0],'FontName','Calibri','FontSize',16);
xlabel(ax,"Frequency [Hz]",'FontName','Calibri','FontSize',12);
ylabel(ax,"PSD [m^2/Hz]",'FontName','Calibri','FontSize',12);
ax.XAxis.TickLabels = ''
ax.XAxis.Label.Visible = 'off'
grid on; box on;
pos = ax.Position;
pos(1) = 0.5; pos(2) = 0.52; pos(4) = 0.42; pos(3)=0.45; ax.Position = pos;
ax2 = axes('Position',pos,'XAxisLocation','top','YAxisLocation','left','Color','none');
ax2.XScale = ax.XScale;
ax2.YScale = ax.YScale;
ax2.XLim=ax.XLim;
ax2.YLim=ax.YLim;
ax2.YTick=ax.YTick;
ax2.XTick=ax.XTick;
ax2.XTickLabel = cellstr(string(round(1./ax2.XTick/2.*1000,0))');
ax2.XAxis.FontSize = 8;
xlabel(ax2,'Sampling Time [ms]','FontName','Calibri','FontSize',12);

% Display the inner tower power spectrum.
ax = subplot(2,2,4);
hold on;
ph = loglog(freqi,psdi);
ph(1).LineStyle = 'none'; ph(1).Marker='o'; ph(1).MarkerFaceColor='w';ph(1).LineWidth=1;ph(1).MarkerSize=3;
ph(2).LineStyle = 'none'; ph(2).Marker='s'; ph(2).MarkerFaceColor='w';ph(2).LineWidth=1;ph(2).MarkerSize=3;
ph(3).LineStyle = 'none'; ph(3).Marker='d'; ph(3).MarkerFaceColor='w';ph(3).LineWidth=1;ph(3).MarkerSize=3;
ph = loglog(freqi,psdsgi);
ph(1).LineWidth = 2; ph(2).LineWidth = 2; ph(3).LineWidth = 2;
lh = yline(arcsec^2); lh.LineWidth=2;
hold off;
ax.XScale = 'log';
ax.YScale = 'log';
ax.XLim=[0.1,30];
ax.YLim=[1e-15,1e-6];
ax.XTick=[0.1,0.2,0.3,0.4,0.5,1,2,3,4,5,10,20,30];
ax.XTickLabel = cellstr(string(round(ax.XTick,1))');
% ax.XTickLabel = cellstr(string(round(1./ax.XTick/2.*1000,0))');
ax.GridLineStyle = ':';
legend(["x","y","z","x(filt)","y(filt)","z(filt)","1'' fwhm"]);
title(ax,"(C) Inner tower (calibrated)",'Units','Normalized', ...
     'Position',[0.15 0.9 0],'FontName','Calibri','FontSize',16);
xlabel(ax,"Frequency [Hz]",'FontName','Calibri','FontSize',12);
ylabel(ax,"PSD [m^2/Hz]",'FontName','Calibri','FontSize',12);
grid on; box on;
pos = ax.Position;
pos(1) = 0.5; pos(2) = 0.09; pos(4) = 0.42; pos(3)=0.45; ax.Position = pos;
