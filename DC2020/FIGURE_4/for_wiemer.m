% Here we test the range of outcome that can be produced from using foreshock the
% traffic light system. Particularly, we test how the outcomes are
% sensitive to uniform priors for:
% 
% - catalog start date
% - fault orientation
% - Mc correction
% 
% These are all allowable within the criteria presented in Gulia and
% Wiemer, 2019. I have left this uncommented, however, the full monte carlo
% simulation is time consuming and therefore commented out

%%

rng(1) = 1;
N    = 1000;

% Here are, what we considered reasonnable parameters to vary
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
decyear2dt = @(x) datetime(datenum(datetime(floor(x), 1, 1) + years(x-floor(x))),'ConvertFrom','datenum');

dt_foreshock64 = datetime(2019,07,04,18,33,00);
T_foreshock64 = decyear(dt_foreshock64); % to the nearest minute

dt_mainshock71 = datetime(2019,07,06,04,19,00);
T_mainshock71 = decyear(dt_mainshock71);
                     
% Here we define the prior distribution - a uniform distribution                     
uniform_prior = @(val_range) rand(1,1) * range(val_range) + val_range(1); % random value between the start and end date
     

%% plot of an example solution for Ridgecrest
% note that the results are very sensitive to the parameter choice
mc_correction = 0.1;
catalog_start_date = 2000;
strike_dip    = [227, 86];
foreshock_blind_time = 0.05;
mainshock_blind_time = 2;
line_opt = {'Color',[0 0 0],'LineWidth',2};

[bgw,tgw] = plot_gw_timeseries('Earthquake','Ridgecrest', ...
                          'fsM',                    6.4, ...
                          'msM',                    7.1, ...
                          'MaxCurvatureCorrection', mc_correction, ... 
                          'CatalogStartDate',       catalog_start_date, ...
                          'ForeshockStrike',        strike_dip(1), ...
                          'ForeshockDip',           strike_dip(2), ...
                          'MainshockStrike',        321, ...
                          'MainshockDip',           81, ...
                          'ForeshockBlindTime',     foreshock_blind_time, ...
                          'MainshockBlindTime',     mainshock_blind_time, ...
                          'PLotOptions',            line_opt, ...
                          'NewFigure',              true, ...
                          'PlotOutput',             true);
                    
                      
%% plot the verion preffered solution by gulia  
mc_correction = 0.2;
catalog_start_date = 2010;
foreshock_blind_time = 0.05;
mainshock_blind_time = 1;
line_opt = {'Color',[0 0 0],'LineWidth',2};
expertChoiceMc.pre = 'none';
strike_dip    = [137, 87];
expertChoiceMc.post= 1.5;
expertChoiceMc.post2=1.5;

[bgw,tgw] = plot_gw_timeseries('Earthquake','Ridgecrest', ...
                          'fsM',                    6.4, ...
                          'msM',                    7.1, ...
                          'MaxCurvatureCorrection', mc_correction, ...
                          'CatalogStartDate',       catalog_start_date, ...
                          'ForeshockStrike',        strike_dip(1), ...
                          'ForeshockDip',           strike_dip(2), ...
                          'MainshockStrike',        321, ...
                          'MainshockDip',           81, ...
                          'ForeshockBlindTime',     foreshock_blind_time, ...
                          'MainshockBlindTime',     mainshock_blind_time, ...
                          'McExpertChoice',         expertChoiceMc, ...
                          'PLotOptions',            line_opt, ...
                          'NewFigure',              false);

  


