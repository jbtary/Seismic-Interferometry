% Test nccf function with GOFAR data
clear

% Data files
dir1 = 'G04mat';
dir2 = 'G06mat';
[info1,info2] = format_name_nccf2(dir1,dir2);

% Original sampling frequencies of data streams in Hz
fs1 = 50; fs2 = 50;
% Decimation factor for each sampling frequency (so that they have the same one)
dfac1 = 5; dfac2 = 5;
% New sampling frequencies
fs1 = fs1/dfac1; fs2 = fs2/dfac2;
if fs1 ~= fs2; disp('Different frequencies for the 2 datastreams. Stopping'); return; end

nfiles = 7*24; % # of files to average for nccf computation
nfiles2 = 1*24; % # of files to increment between nccf computations
norm1b = 0; % Flag for 1-bit normalization
funcle = 40; % Nccf length in seconds
segle = 4*funcle*fs1; % Data segment length in samples (adjust it so that there is no loss of data at the end), 50% overlap used

% Band-pass filtering
fNy = fs1/2;
fcl = 1; fch = 2; % Frequency bounds in Hz
if fch - fcl < 1
    h = fdesign.peak('N,F0,BW', 10, fcl+((fch-fcl)/2), fch-fcl, fs1); % peak filter
    Hd = design(h, 'butter');
else
    [b,a]=butter(6,[fcl fch]/fNy,'bandpass'); % Butterworth filter of order 6
end

% Comments:
% Without orientation, only ZZ definitively corresponds to a non-zero term
% of the Green's function for Rayleigh waves (Stehly et al., 2007)

% Loop on number of files (should be the same for both data streams, and 
% same timing as well)
kk = 1;
for ii = 1:nfiles2:length(info1)-nfiles
    
    % 1 - computation of nccf
    ll = 1;
    for jj = ii:ii+nfiles-1
        data1 = load([dir1 '/' info1(jj).name],'data','time');
        datatmp1 = data1.data'; time1 = data1.time; clear data1
        data2 = load([dir2 '/' info2(jj).name],'data','time');
        datatmp2 = data2.data'; time2 = data2.time; clear data2
        
        date1 = time1(1); date2 = time2(1);
        clear time1 time2
        
        % Correction of time delay between data1 and data2 (should be
        % small! Of the order of a few samples max)
        dtdat = (date1 - date2)*86400; % Time delay in seconds
        datatmp2 = delayseq(datatmp2,dtdat,fs2*dfac2); % Correction of delay
        clear dtdat date1 date2
        
        datatmp1 = detrend(datatmp1);
        datatmp2 = detrend(datatmp2);
        datatmp1 = decimate(datatmp1,dfac1);
        datatmp2 = decimate(datatmp2,dfac2);
        
        % Check on the length of both data streams
        if length(datatmp1) ~= length(datatmp2)
            disp(['File ' num2str(jj) ' Size data1z: ' num2str(size(datatmp1))...
                ' size data2z: ' num2str(size(datatmp2))])
            clear datatmp*
            continue
        end
        
        % Filtering
        if fch - fcl < 1
            datatmp1 = filtfilt(Hd.sosMatrix,Hd.ScaleValues,datatmp1);
            datatmp2 = filtfilt(Hd.sosMatrix,Hd.ScaleValues,datatmp2);
        else
            datatmp1 = filtfilt(b,a,datatmp1);
            datatmp2 = filtfilt(b,a,datatmp2);
        end
        
        % One-bit normalization
        if norm1b == 1
            %         data1 = sign(data1);
            %         data2 = sign(data2);
            % Alternative: normalize by envelope (L. Stehly)
            envd1 = hilbert(datatmp1); envd1 = abs(envd1);
            envd2 = hilbert(datatmp2); envd2 = abs(envd2);
            datatmp1 = datatmp1./envd1;
            datatmp2 = datatmp2./envd2;
            clear envd1 envd2
        end
        
%         % Messing up the data with zeros
%         numb1 = randi(length(datatmp1)-101,[10,1]);
%         numb2 = randi(100,[10,1]);
%         for nn = 1:10; datatmp1(numb1(nn):numb1(nn)+numb2(nn)-1,1) = zeros(numb2(nn),1); end;
%         clear numb*
%         numb1 = randi(length(datatmp1)-101,[10,1]);
%         numb2 = randi(100,[10,1]);
%         for nn = 1:10; datatmp2(numb1(nn):numb1(nn)+numb2(nn)-1,1) = zeros(numb2(nn),1); end;
%         clear numb*
        
        % Loop on data segments
        noverlap = round(0.5*segle); % 50% overlap used
        indb = 1:(segle-noverlap):length(datatmp1)-(segle-1);
        inde = segle:(segle-noverlap):length(datatmp1);
        for nn = 1:length(indb)
            
            [xctmp,lag] = xcorr(datatmp1(indb(nn):inde(nn),1),...
                datatmp2(indb(nn):inde(nn),1),'coeff');
            xctmp = xctmp/max(abs(xctmp));
            xcdtmp(:,ll) = xctmp(round(-funcle*fs1) + floor(end/2)+1:...
                round(funcle*fs1) + floor(end/2)+1);
            ll = ll+1;
            
            clear xctmp
        end
        
        clear datatmp* noverlap indb inde nn
    end
    
    xcd(:,kk) = mean(xcdtmp,2,'omitnan'); % Avg nccf
    xcd(:,kk) = xcd(:,kk)/max(abs(xcd(:,kk))); % Normalize
%     xcdsv(:,:,kk) = xcdtmp;
    
    kk = kk+1;
    clear jj ll xcdtmp
    disp([num2str(ii) '/' num2str(length(info1)-nfiles)])
end
lag = lag(round(-funcle*fs1) + floor(end/2)+1:round(funcle*fs1) + floor(end/2)+1);
xc_avg = mean(xcd,2,'omitnan'); % Avg nccf

figure
ax(1) = axes('Position',[0.1 0.1 0.8 0.6]);
imagesc(lag*(1/fs1),1:size(xcd,2),xcd'); colormap gray
xlabel('Lag time (s)'); ylabel('nccf #');
ax(2) = axes('Position',[0.1 0.7 0.8 0.2]);
plot(lag*(1/fs1),xc_avg,'k'); ylabel('Avg. coeff.')
set(gca,'XTickLabel',[])
linkaxes(ax,'x');

clear kk ax h Hd

%% Test clock drift measure with cross-spectrum (nccf.m)

% clear
% load nccf_OBS01_03_CompZ_fcl0_25_fch0_5_1bitNorm
xcdsv = xcd;

% Window of the NCCFs used in the calculation of the time lag
% Values for OBS 1 and 3: winb = 160*fs; wine = 230*fs;
figure(111)
plot(lag*(1/fs),xc_avg,'k'); ylabel('Avg. coeff.'); xlabel('Lag time (s)');
disp('Pick beginning and end of window of wave arrival.')
disp('Keep in mind that nccf might not be symmetric.')
disp('Zoom and press enter')
pause;
[tb,~] = ginput;

winb = round(min(abs(tb)))*fs;
wine = round(max(abs(tb)))*fs;
lag0 = floor(length(lag)/2)+1;
ccthr = 0.85;
indtrack = 1:size(xcd,2); % Track of indexes
nit = 10;
clear tb; close figure 111

for it = 1:nit
    
    xc_avg = mean(xcd,2); % Avg nccf
    xc_pos = xcd(lag0+winb:lag0+wine,:);
    xc_neg = xcd(lag0-wine:lag0-winb,:);
    xc_avg_pos = xc_avg(lag0+winb:lag0+wine,:);
    xc_avg_neg = xc_avg(lag0-wine:lag0-winb,:);
    
    xc_avg_pos = xc_avg_pos - mean(xc_avg_pos);
    xc_avg_neg = xc_avg_neg - mean(xc_avg_neg);
    
    for ii = 1:size(xcd,2)
        
        xc_pos(:,ii) = xc_pos(:,ii) - mean(xc_pos(:,ii));
        
        % Time delay based on quadratic interp of cross-corr function
        [R,L]=xcorr(xc_avg_pos,xc_pos(:,ii),'coeff');
        [Rmax,Imax]=max(R);
        P=polyfit([L(Imax-1) L(Imax) L(Imax+1)], [R(Imax-1) R(Imax) R(Imax+1)],2);
        cdyn(ii,1) = (-P(2)/(2*P(1)))/fs;
        cc(ii,1) = Rmax;
        iv = find(R > ccthr);
        if isempty(iv) == 1
            cc(ii,2) = -99;
        else
            cc(ii,2) = (L(iv(end)) - L(iv(1)))*(1/fs);
        end
        clear R L Rmax Imax P iv
        
        xc_neg(:,ii) = xc_neg(:,ii) - mean(xc_neg(:,ii));
        
        [R,L]=xcorr(xc_avg_neg,xc_neg(:,ii),'coeff');
        [Rmax,Imax]=max(R);
        P=polyfit([L(Imax-1) L(Imax) L(Imax+1)], [R(Imax-1) R(Imax) R(Imax+1)],2);
        cdyn(ii,2) = (-P(2)/(2*P(1)))/fs;
        cc(ii,3) = Rmax;
        iv = find(R > ccthr);
        if isempty(iv) == 1
            cc(ii,4) = -99;
        else
            cc(ii,4) = (L(iv(end)) - L(iv(1)))*(1/fs);
        end
        clear R L Rmax Imax P iv
        
        cdyn(ii,3) = (cdyn(ii,1) - cdyn(ii,2))/2;
        cdyn(ii,4) = (cdyn(ii,1) + cdyn(ii,2))/2;
        
    end
    
    % Find daily NCCFs that are correlated > 0.6 with the averaged one
    if it == 1
        ccmax = max(cc(:,1),cc(:,3));
        iv = find(ccmax > 0.6);
        xcd = xcd(:,iv); indtrack = indtrack(iv);
        cdyn = cdyn(iv,:);
    end
    
    % Correct time shift for nccfs that are > cc avg
    for ii = 1:size(xcd,2)
        xcd(:,ii) = delayseq(xcd(:,ii),cdyn(ii,3),fs);
    end
    
    csv(:,it) = cdyn(:,3);
    
    cdynf = cdyn;
    clear ccmax iv cdyn xc_*
end

xc_avg = mean(xcd,2); % Avg nccf
xc_avg_pos = xc_avg(lag>0,:);
xc_avg_neg = xc_avg(lag<0,:);

xc_avg_pos = xc_avg_pos - mean(xc_avg_pos);
xc_avg_neg = xc_avg_neg - mean(xc_avg_neg);

% Switch up down neg avg nccf (!)
[R,L]=xcorr(xc_avg_pos,flipud(xc_avg_neg),'coeff');
[~,Imax]=max(R);
P=polyfit([L(Imax-1) L(Imax) L(Imax+1)], [R(Imax-1) R(Imax) R(Imax+1)],2);
csta = (-P(2)/(2*P(1)))/fs;
clear R L Rmax Imax P xc_*

cabs = csta + cdynf(:,3);

xcd = xcdsv;
clear xcdsv

% Test avec vda (nccfvda.m)

% clear
% load nccf_OBS01_03_CompZ_fcl0_25_fch0_5_1bitNorm

xc_avg = mean(xcd,2); % Avg nccf

% Parameters
ccthr = 0.85; % Threshold on cc for definition of fit
win = 10*fs; % Window in samples
noverlap = round(0.5*win);

% Window of the NCCFs used in the calculation of the fit (need to be
% customized for each OBS pair). Following values are for OBS 1 and 3
winb = 186*fs;
wine = 201*fs;
lag0 = floor(length(lag)/2)+1;

indb = 1:(win-noverlap):length(xc_avg)-(win-1);
inde = win:(win-noverlap):length(xc_avg);

% Center of windows
tind = win-noverlap+1:win-noverlap:length(xc_avg)-(win-noverlap-1);
iw1 = find(tind>(lag0-wine) & tind<(lag0-winb));
iw2 = find(tind>(lag0+winb) & tind<(lag0+wine));
iw = [iw1 iw2]; clear iw1 iw2

for ii = 1:size(xcd,2)
    
    xcd(:,ii) = xcd(:,ii) - mean(xcd(:,ii));
    
    for jj = 1:length(indb)
        
        % Time delay based on quadratic interp of cross-corr function
        [R,L]=xcorr(xc_avg(indb(jj):inde(jj)),xcd(indb(jj):inde(jj),ii),'coeff');
        [Rmax,Imax]=max(R);
        P=polyfit([L(Imax-1) L(Imax) L(Imax+1)], [R(Imax-1) R(Imax) R(Imax+1)],2);
        dttmp(jj,1) = (-P(2)/(2*P(1)))/fs;
        cc(jj,1) = Rmax;
        iv = find(R > ccthr); % Confidence interval
        if isempty(iv) == 1
            cc(jj,2) = NaN;
        else
            cc(jj,2) = (L(iv(end)) - L(iv(1)))*(1/fs);
        end
        clear R L Rmax Imax P iv
        
    end
    
    % Selection of points with good cross corr, possibility of reweighing
    % the linear fit with confidence interval, use confidence interval in
    % selection also?
    iv = find(cc(:,1)>ccthr);
    ind = intersect(iv,iw);
    if isempty(ind) == 1 || length(ind) == 1
        fit(ii,:) = [NaN NaN];
    else
        fit(ii,:) = polyfit((tind(ind)-(floor(length(xc_avg)/2)+1))'*(1/fs),dttmp(ind),1);
    end
    % Add confidence on regression, taking into account number of points?
    dtsv(:,ii) = dttmp;
    cc1sv(:,ii) = cc(:,1);
    cc2sv(:,ii) = cc(:,2);
    drift(ii,1) = polyval(fit(ii,:),0); % Drift for lag=0
    
    clear cc dttmp jj iv ind
end

%%
for jj = 1:length(indtrack)
    nn = indtrack(jj);
    
    iv = find(cc1sv(:,nn)>ccthr);
    ind = intersect(iv,iw);
    
    figure
    litr = polyval(fit(nn,:),lag*(1/fs));
    ax(1) = axes('Position',[0.1 0.1 0.8 0.3]);
    [axf,h1,h2] = plotyy((tind-(floor(length(xc_avg)/2)+1))*(1/fs),...
        cc1sv(:,nn),lag*(1/fs),litr,'bar','plot');
    h1.BaseValue = 0.8; h1.EdgeColor = 'w'; h1.FaceColor = [.7 .7 .7];
    h2.LineStyle = '--'; h2.Color = 'k';
    axf(1).YLim = [0.7 1]; axf(1).YTick = 0.7:0.1:1;
    axf(2).YColor = [0 0 0]; axf(2).YLim = [min(dtsv(ind,nn))-0.25 max(dtsv(ind,nn))+0.25];
    axf(2).YTick = round(min(dtsv(ind,nn)),1)-0.25:0.2:round(max(dtsv(ind,nn))+0.25,1);
    set(get(axf(1),'Ylabel'),'String','CC coeff.')
    set(get(axf(2),'Ylabel'),'String','\delta t (s)')
    xlim([min(lag)*(1/fs) max(lag)*(1/fs)])
    xlabel('Lag time (s)')
    
    ax(2) = axes('Position',[0.1 0.1 0.8 0.3]);
    plot((tind-(floor(length(xc_avg)/2)+1))*(1/fs),...
        dtsv(:,nn),'k*')
    hold on
    errorbar((tind(ind)-(floor(length(xc_avg)/2)+1))*(1/fs),...
        dtsv(ind,nn),cc2sv(ind,nn)/2,'r*')
    ylabel('dt (s)'); xlabel('Time (s)')
    axis([min(lag)*(1/fs) max(lag)*(1/fs) min(dtsv(ind,nn))-0.25 max(dtsv(ind,nn))+0.25])
    axis off
    
    ax(3) = axes('Position',[0.1 0.45 0.8 0.2]);
    plot(lag*(1/fs),xcd(:,nn),'k'); ylabel(['NCCF #' num2str(nn)])
    set(gca,'XTickLabel',[]); ylabel('nccf');
    xlim([min(lag)*(1/fs) max(lag)*(1/fs)])
    
    ax(4) = axes('Position',[0.1 0.7 0.8 0.2]);
    plot(lag*(1/fs),xc_avg,'k'); ylabel('Avg. coeff.')
    set(gca,'XTickLabel',[]); xlim([min(lag)*(1/fs) max(lag)*(1/fs)])
    title(['NCCF #' num2str(nn)])
    linkaxes(ax,'x');
    
%     filename = ['NCCF_OBS01_03_CompZ_fcl0_25_fch0_5_NobitNorm_VDA_NCCF' num2str(jj)];
%     print('-dpdf','-r300',filename)
    
    clear ax* ii nn iv ind hAx hLine1 hLine2 filename
    
end

%% Comparison between TSA and VDA results

tdays = (nfiles/24)/2:(nfiles/24):(size(info1,1)-(nfiles/2))/24;
figure
axes('Position',[0.1 0.15 0.7 0.3]);
plot(tdays(indtrack),cdynf(:,1),'ko','MarkerSize',4) % TSA drifts Pos lags
hold on
plot(tdays(indtrack),cdynf(:,2),'k*','MarkerSize',5) % TSA drifts Neg lags
plot(tdays(indtrack),cdynf(:,3),'k.','MarkerSize',10) % TSA drifts diff/2
plot(tdays(indtrack),cdynf(:,4),'kd','MarkerSize',4) % TSA drifts average

plot(tdays,drift,'ks','MarkerSize',4) % VDA drifts
legend('TSA^+','TSA^-','TSA^+-TSA^- /2','TSA^++TSA^- /2','VDA',...
    'Location','eastoutside')
legend boxoff

xlim([0 max(tdays)+((nfiles/24)/2)])
ylabel('\delta t (s)'); xlabel('Time (days)'); title('\delta t dyn')
clear tdays
