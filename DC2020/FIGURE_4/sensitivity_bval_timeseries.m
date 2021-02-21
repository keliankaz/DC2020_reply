% here we test the range of outcome that can be produced from using foreshock the
% traffic light system. Particularly, we test how the outcomes are
% sensitive to uniform priors for:

% - catalog start date
% - fault orientation
% - Mc correction
% 
% These are all allowable withing the strict criteria presented in Gulia and Wiemer

rng(1) = 1;
N    = 1000;

% Here are, what we consider, reasonnable parameters to vary
% variable              range               
mc_correction         = [0.1,0.3];                
catalog_start_date    = [2000,2012]; % note that later results in too few earthquakes to reach the 50 event threshold

% Strike and dip options
% (gCMT: 201907041733B CENTRAL CALIFORNIA, 201907060319A CENTRAL CALIFORNIA) 
% https://www.globalcmt.org/CMTsearch.html
foreshock_strike_dip_opt        = [227, 86; ...
                                   137, 87];
mainshock_strike_dip_opt        = [321, 81];
foreshock_blind_time  = [0.01,0.5];
mainshock_blind_time  = [0.5,5];

% Timing:

decyear = @(x) year(x) + years(x - dateshift(x, 'start', 'year'));
decyear2dt = @(x) datetime(datenum(datetime(floor(x), 1, 1) + years(x-floor(x))),'ConvertFrom','datenum');

dt_foreshock64 = datetime(2019,07,04,18,33,00);
T_foreshock64 = decyear(dt_foreshock64); % to the nearest minute

dt_mainshock71 = datetime(2019,07,06,04,19,00);
T_mainshock71 = decyear(dt_mainshock71);
                     
% Here we define the prior distribution - a uniform distribution                     
uniform_prior = @(val_range) rand(1,1) * range(val_range) + val_range(1); % random value between the start and end date
                        

figure; 
subplot(1,3,[1 2]);hold on
line_opt = {'Color',[0.8 0.8 0.8,0.1],'LineWidth',1};
[T,B] = deal(cell(1,N));
[foreshockWarning,aftershockWarning] = deal(zeros(1,N));

%%
parfor n = 1:N
    
    mci         = uniform_prior(mc_correction);
    catStart    = uniform_prior(catalog_start_date);
    fsbt        = uniform_prior(foreshock_blind_time);
    msbt        = uniform_prior(mainshock_blind_time);
    fssdi       = foreshock_strike_dip_opt(randi(2),:);
    [bi,ti] = plot_gw_timeseries(   'Earthquake','Ridgecrest', ...
                          'fsM',                   6.4, ...
                          'msM',                   7.1, ...
                          'MaxCurvatureCorrection', mci, ...
                          'CatalogStartDate',       catStart, ...
                          'ForeshockStrike',        fssdi(1), ...
                          'ForeshockDip',           fssdi(2), ...
                          'MainshockStrike',        mainshock_strike_dip_opt(1), ...
                          'MainshockDip',           mainshock_strike_dip_opt(2), ...
                          'ForeshockBlindTime',     fsbt, ...
                          'MainshockBlindTime',     msbt, ...
                          'PlotOutput',             false, ...
                          'NewFigure',              false);

     B{n} = bi;
     T{n} = ti;      
     foreshockWarning(n) = nanmedian(bi(ti<T_mainshock71));
     aftershockWarning(n) = nanmedian(bi(ti>T_mainshock71));
     
end

%%

decyear2dt = @(x) datetime(datenum(datetime(floor(x), 1, 1) + years(x-floor(x))),'ConvertFrom','datenum');
plot_ts = @(t,b) plot(decyear2dt(t),b,'-','Color',[0.5,0.5,0.5,min(10/N,1)],'Linewidth',0.1);
figure;subplot(1,4,[1 2]);hold on
tpos = [-0.15,1];
t = title('a)'); set(t,'Position',[tpos(1),tpos(2)-0.05],'Units', 'normalized')

cellfun(plot_ts,T,B);
tLimTemp = xlim;
xlim([tLimTemp(1),max(cellfun(@(t) max(decyear2dt(t)), T))])
ylabel('Relative b-value (%)')


cfs = [0.8500    0.3250    0.0980];
cms = [0.9290    0.6940    0.1250];

correctI = foreshockWarning<90 & aftershockWarning>110;
sortaCorrectI = foreshockWarning<110 & aftershockWarning>90;
correctWarningPrecent = sum(correctI)/N * 100;
sortaCorrectWarningPercent = sum(sortaCorrectI)/N * 100;

hold on
ylim([50,150])
brel_lims = ylim;
tLim = xlim;

plot(tLim,90*ones(1,2),'--', ...
                       'Color',[0.6350    0.0780    0.1840], ...
                       'Linewidth',2);
plot(tLim,110*ones(1,2),'--', ...
                        'Color',[0.4660    0.6740    0.1880], ...
                        'Linewidth',2);

brelLim = ylim;
plot([dt_mainshock71,dt_mainshock71],brelLim,'-','LineWidth',2,'Color',[cms,0.6]);
plot([dt_foreshock64,dt_foreshock64],brelLim,'-','LineWidth',2,'Color',[cfs,0.6]);

%%

% we can also visualize the distribution of warnings that may have resulted
% during the foreshock and aftershock period provided these results

subplot(1,4,3); hold on;
t = title('b)'); set(t,'Position',[tpos(1),tpos(2)-0.05],'Units', 'normalized')

E = (50:5:150);
hfs = histogram(foreshockWarning,E,'FaceColor',cfs,'EdgeColor',[1 1 1]);
hms = histogram(aftershockWarning,E,'FaceColor',cms,'EdgeColor',[1 1 1]);
% histogram(foreshockWarning(correctI),E,'faceColor',cfs)
% histogram(aftershockWarning(correctI),E,'faceColor',cms);
xlim(brel_lims)
N_lims = ylim;
plot([90,90],N_lims,'--', ...
                       'Color',[0.6350    0.0780    0.1840], ...
                       'Linewidth',2);
plot([110,110],N_lims,'--', ...
                        'Color',[0.4660    0.6740    0.1880], ...
                        'Linewidth',2);

%xlabel('Relative b-value (%)')
xticklabels([])
yticklabels([])
lgh = legend([hfs,hms],{'Foreshock','Aftershock'},'Location','southeast');
title(lgh,'Median b_{rel}')
% ylabel('Number of simulations')
view([90 -90])


%% 

subplot(1,4,4); hold on
t = title('c)'); set(t,'Position',[tpos(1),tpos(2)-0.05],'Units', 'normalized')
histogram(aftershockWarning-foreshockWarning,'EdgeColor',[1 1 1]);
xlabel('\Delta b_{rel}')
yticklabels([])

ftsz    = @(fh,fontSize) set(findall(fh,'-property','FontSize'),'FontSize',fontSize);

setsize = @(fh,dim1,dim2) set(fh,...
    'Units',        'Inches', ...
    'Position',     [0,0,dim1,dim2],...
    'PaperUnits',   'Inches',...
    'PaperSize',    [dim1,dim2]);

ftsz(gcf,10);
setsize(gcf,8,2.5);


