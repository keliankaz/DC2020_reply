%% figure 2 Ridgecrest
figure
eq_of_interest = 'Ridgecrest_Foreshock';
Run_TLS_Gulia_Wiemer_reproduce_GW
xlabel({'','M_W 6.4 source volume'},'FontSize',10,'FontWeight','bold')
 
eq_of_interest = 'Ridgecrest_Mainshock';
Run_TLS_Gulia_Wiemer_reproduce_GW
xlabel({'','M_W 7.1 source volume'},'FontSize',10,'FontWeight','bold')


%% figure 3 Puerto Rico
figure
eq_of_interest = 'PuertoRico_Foreshock';
Run_TLS_Gulia_Wiemer_reproduce_GW

% minor tweeks
xlabel({'','M_W 5.0 source volume'},'FontSize',10,'FontWeight','bold')
set(AX, 'xlim',[2,7], ...
        'ylim',[10^-2,10^5])

eq_of_interest = 'PuertoRico_Mainshock';
Run_TLS_Gulia_Wiemer_reproduce_GW

% minor tweeks:
xlabel({'','M_W 6.4 source volume'},'FontSize',10,'FontWeight','bold')
set(AX, 'xlim',[2,7], ...
        'ylim',[10^-2,10^5])
