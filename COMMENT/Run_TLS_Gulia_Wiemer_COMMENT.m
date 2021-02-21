% ------- Traffic Light System by Laura Gulia and Stefan Wiemer -----------
% -------------------------------------------------------------------------

% The method we propose is a real-time tool to discriminate between foreshocks and aftershocks

% This code calculates a background reference b-value for a subcatalog of events.
% This background b-value is used to calculate the variation in percentage
% from the reference, following a Magnitide M>=6. The analysis is performed using a moving-window
% approach, moved forward event by event. The window length is defined by paratemeters Npre and Npost;
% after the M>=6 quake, due to completeness issue, the removal a 'blackout
% period'(defined in the input parameter section) is required for each sequence,
% an input based on data quality analysis. 
%
% We then define a simple level of concern, expressed as a traffic light, depending on such value:

% Red  <10% change
% Yellow  +/-10% change
% Green >10% change

% The events in the subcatalog are selected according to the method by Gianfranco Vannucci, described in
% Gulia et al. (GRL, 2018)
%
% CopyRight: Laura Gulia & Stefan Wiemer, ETH Zurich
% date: 2019, tested under Matlab 2018a 
%
%% 

 %% clear
 addpath ./src
 warning off

 
% swtched the starting section to a switch/case format to avoid clumsy
% comment uncomment format
% eq_of_interest = 'Ridgecrest_Mainshock';
% eq_of_interest = 'Ridgecrest_Mainshock';



switch eq_of_interest
    case 'Ridgecrest_Foreshock'
        % for ridgecrest, M6.4 (i.e. the foreshock)%
        load ./cat/subcat/Ridgecrest/rev/Ridgecrest_foreshock_ms.mat
        load ./cat/subcat/Ridgecrest/rev/Ridgecrest_foreshock.mat
        Tms         =2.019506114559868e+03;   %time of the M6.4
        Tms2        =2.019509969332192e+03;  %time of the M7.1
        days        =0.05; % removal of days after the first M>=6.4 event
        days_post2  =2; % removal of days after the second M>=7.1 event
        LE_Tms      =2.019506112772704e+03; %last event before rc, M6.4
        LE_Tms2     =2.019509968955797e+03; %last event before , M7.1
        alt_cat_name= './cat/subcat/Ridgecrest/rev/Ridgecrest_regional.mat';
        MMS = 6.4;
        % Define the time-lenght of the catalog
        Tmin=2000;
        Tmax=2020;
        saveName = 'foreshock'; % kdc

    case 'Ridgecrest_Mainshock'
        % for ridgecrest, M7.1 (i.e. the mainshock)%      
        load ./cat/subcat/Ridgecrest/rev/Ridgecrest_mainshock_ms.mat
        load ./cat/subcat/Ridgecrest/rev/Ridgecrest_mainshock.mat
        Tms         =2.019506114559868e+03;   %time of the M6.4
        Tms2        =2.019509969332192e+03;  %time of the M7.1
        days        =0.05; % removal of days after the first M>=6.4 event
        days_post2  =2; % removal of days after the second M>=7.1 event
        LE_Tms      =2.019506112772704e+03; %last event before rc, M6.4
        LE_Tms2     =2.019509968955797e+03; %last event before , M7.1
        alt_cat_name= './cat/subcat/Ridgecrest/rev/Ridgecrest_regional.mat';
        MMS = 7.1;
        
        % Define the time-lenght of the catalog
        Tmin=2000;
        Tmax=2020;
        saveName = 'mainshock'; % kdc
        
    case 'PuertoRico_Foreshock'
        % for Puerto Rico, M5.0 (i.e. the foreshock)%
        load ./cat/subcat/PuertoRico/rev/PuertoRico_foreshock_ms.mat
        load ./cat/subcat/PuertoRico/rev/PuertoRico_foreshock.mat
        Tms         =2.019991906396816e+03;     %time of the M5
        Tms2        =2.020017350592301e+03;     %time of the M6.4
        days        =0.05;                      % removal of days after the first M>=5.0 event
        days_post2  =0.05;                      % removal of days after the second M>=6.4 event
        LE_Tms      =2.019991877045599e+03;     %last event before fs, M5
        LE_Tms2     =2.020017340344819e+03;     %last event before , M6.4
        alt_cat_name= './cat/subcat/PuertoRico/PuertoRico_regional.mat';
        MMS         = 5.0;
        
        % Define the time-lenght of the catalog
        Tmin=2003;
        Tmax=2021;
        saveName = 'foreshock'; % kdc
        
    case 'PuertoRico_Mainshock'
        % % for Puerto Rico, M6.4 (i.e. the mainshock)%
        load ./cat/subcat/PuertoRico/rev/PuertoRico_mainshock_ms.mat
        load ./cat/subcat/PuertoRico/rev/PuertoRico_mainshock.mat
        Tms         =2.019991906396816e+03;     %time of the M6.4
        Tms2        =2.020017350592301e+03;     %time of the M7.1
        days        =0.05;                      % removal of days after the first M>=5.0 event
        days_post2  =0.05;                      % removal of days after the second M>=6.4 event
        LE_Tms      =2.019991877045599e+03;     %last event before rc, M5
        LE_Tms2     =2.020017340344819e+03;     %last event before , M6.4
        alt_cat_name= './cat/subcat/PuertoRico/PuertoRico_regional.mat';
        MMS         = 6.4;
        
        % Define the time-lenght of the catalog
        Tmin=2003;
        Tmax=2021; 
        saveName = 'mainshock'; % kdc
end

%% 
%**************************************************************************
%%%%%%%%%%%%%%%%%%%%% Common input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

catS(:,6) = round(catS(:,6)*10)/10; %kdc

% define the correction to be applied to the Maximum Curvature Method
% following (Wiemer and Woessner, 2005). We allow to apply different pre-
% and post mainshock corrections but the defaults is equal = 0.2 
corr_post=0.2;
corr_pre=0.2;

%Magnitude binning
mbin    =0.1;
% Window lenght for the time-series: Npre, before the M>=6 event, Npost
% after  (due to the event abundance, we can use larger sample sizes during
% the rich aftershock sequence to reduce uncertainties. Note these sample
% sizes are before cutting at the subsequent define completeness levels. 
Npre    =250;
Npost   =400;

% in the case the subcatalog does not contain Npre events, a new subcatalog is created selecting the closest Nreg events to the hypocenter
Nreg    =250;

% Minimum number of event above the magnitude of completeness to estimate a
% the b-value. This is minimum that applies after cutting the sample size
% Npre or Npost at the determinej Mc. 
Nmin    =50;


% Define Latitude, Longitude, Magnitude, Depth and Time of the event
LATms=ms(1,2);
LONms=ms(1,1);
Mms=ms(1,6);
Monthms=ms(1,4);
Dayms=ms(1,5);
Depthms=ms(1,7);
%% 
% Some housekeeping parts and preparatory work
% 
% sorted catalog by time. 
catS=sortrows(catS,3);

% cut the catalog at Tmin to not include in the calculation of the background level:
% 1. in the case of Amatrice and Norcia, the aftershocks of L'Aquila (M6.3, 2009)
% 2. in the case of Kumamoto, the first months of the Tohoku aftershocks
ll=catS(:,3)>=Tmin;
catS=catS(ll,:);

% querry magnitude range 
Mrange=min(catS(:,6)):mbin:max(catS(:,6));

% The catalog is then divided in 3 sub-catalogs: background (catS_pre), in between (catS_post) and after the second M>=6 event (catS_post2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ll=catS(:,3)<=LE_Tms; %last event before the first M>=6 event
catS_pre=catS(ll,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %remove the first X dayfs after the mainshocks, Expert Choice (# defined above)
hh=catS(:,3)>Tms+days/365;
catS_post=catS(hh,:);
catS_post1=catS_post;   
% catS_post(1,6) % kdc

hh=catS_post(:,3)<=LE_Tms2; %%last event before the second M>=6 event
catS_post=catS_post(hh,:);
% catS_post(1,6) % kdc


% Define a raw magnitude of completeness (Mc_bulk) to pre-cut the catalog: this allow the Maximum Curvature to perform better.
% For the magnitude of completeness, conservative choices are always adopted
% calls function Mc_bulk2
[Mc_pre,cntnumb_pre]=Mc_bulk2(catS_pre, Mrange,0); % kdc
% [Mc_pre,cntnumb_pre]=Mc_bulk2(catS_pre, Mrange,corr_post); %kdc (actually correspond to the corr_post defined above)
catSpre_orig=catS_pre;
Mcpre_orig=Mc_pre;

% The same for the post period 
cat_test=catS_post;
[Mc_post,cntnumb_post]=Mc_bulk2(catS_post, Mrange,0); % kdc
% [Mc_post,cntnumb_post]=Mc_bulk2(catS_post, Mrange,corr_post); % kdc (actually correspond to the corr_post defined above)
% Mcpost_orig= 1.5; % kdc - uncomment to reproduce Gulia - ridgecrest
Mcpost_orig= Mc_post; 

% do the cutting
ll=catS_pre(:,6)>= Mc_pre; % kdc
catS_pre=catS_pre(ll,:);
catS_pre_orig=catS_pre;

ll=catS_post(:,6)>=Mcpost_orig; 
catS_post_orig=catS_post(ll);


%% 
% Step 1; Estimate the reference b-value
if length(catS_pre(:,1))>= Npre  % Do the analysis if we have enough events 
    
    disp('catS_pre: sequence')
    
    % background reference b-value estimation: the median of all the bpre estimations
    [prpre,Mcpre, bpre, apre, tpre,bprevalue,timewindowpre,result_flag_pre,bpreplus,bpreminus]...
        =abwithtime_pre(catS_pre,Npre,mbin,Mms,corr_pre);
 
    
else % Not enough events, default to regional 
    % in such case (here for the Kumamoto foreshock M6.5), since we cannot estimate a b-value using the events
    % inside the box, we upload a reduced version of the JMA catalog, to calculate a regional bvalue as the reference one
    
    disp('catS_pre=cat_Reg: regional')
    
    load(alt_cat_name)
    cat_Reg = cat_Reg(cat_Reg(:,3)<Tms,:);
    
    dists2d_Npre    = deg2km(distance(LATms,LONms,cat_Reg(:,2),cat_Reg(:,1)));
    [dists_sort,is] = sort(dists2d_Npre);
    asort_Npre      = cat_Reg(is,:);
    
    if length(cat_Reg(:,1))>=Nreg
        cat_Reg         =asort_Npre(1:Nreg,:);
        Rreg            =dists_sort(Nreg,:);
    end
    cat_Reg     = sortrows(cat_Reg,3);
    sizecat     = size(cat_Reg);
  
    [Mcpre_orig, bpre, apre,bprevalue,catSpre_orig]=breg_NLI(cat_Reg,mbin,Tms,Mms,corr_pre);
    
end

% Define the reference b-value 
breference=round(nanmean(bpre),1);

%
% b-value time series following the first M>=6 event
disp('catS_post: sequence')

% call the function abwithtime for the post- event 1 sequence data that will compute the output time series. 
[prpost,Mcpost,bpost,apost,tpost,bpostvalue,timewindow,result_flag_post, ...
    bpostplus,bpostminus]=abwithtime_post(catS_post,Npost,mbin,Mms,corr_post);

hh=catS(:,3)>Tms2+days_post2/365;
catS_post2=catS(hh,:);

% [Mc_post2,cntnumb_post2]=Mc_bulk2(catS_post2, Mrange,corr_post); % kdc
[Mc_post2,cntnumb_post2]=Mc_bulk2(catS_post2, Mrange,0); % kdc

% Mcpost_orig2=1.5; % Reproduce Gulia - Ridgecrest
Mcpost_orig2= Mc_post2; 

%ll=catS_post2(:,6)>=1.8; % kdc ERROR in original code?
ll=catS_post2(:,6)>=Mcpost_orig2; % kdc
catS_post2=catS_post2(ll,:);
catS_post_orig2=catS_post2;

% b-value time series following the second M>=6 event
% call the function abwithtime for the post- event 2 sequence data that will compute the output time series. 

[prpost2,Mcpost2,bpost2,apost2,tpost2,bpostvalue2,timewindow2,result_flag_post2, ...
    bpostplus2,bpostminus2]=abwithtime_post(catS_post2,Npost,mbin,Mms,corr_post);

[bpostmax2,ind]=max(bpost2);
tpostmax2=tpost2(ind);
magcopostmax2=Mcpost2(ind);
apostmax2=apost2(ind);

%% Now plot the summry of the results in one figure 
%plot_summary_figure   %the summary figure also plot the TLS level
plot_summary_fig_kdc

save([saveName,'_workspace.mat'])

