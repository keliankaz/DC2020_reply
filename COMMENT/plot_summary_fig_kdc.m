
%% determine whether plotting  foreshock or mainshock

if strcmp(saveName,'foreshock')
   figureInfo = [];
   figureInfo(2).subplot = {2,2,1};
   figureInfo(2).abc = 'a)';
   figureInfo(3).subplot = {2,2,3};
   figureInfo(3).abc = 'c)';
elseif strcmp(saveName,'mainshock')
   figureInfo = [];
   figureInfo(2).subplot = {2,2,2};
   figureInfo(2).abc = 'b)';
   figureInfo(3).subplot = {2,2,4};
   figureInfo(3).abc = 'd)';
else
    error('Not foreshock or mainshock')
end

tpos = [-0.15,1]; % a-b-c- posistion
c= [0.7       0.7       0.7
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

%%
tLim        = [Tms-2/365, max(tpost2)];
bPlotRange  = [0.4 , 1.25];

subplot(figureInfo(2).subplot{:});

t = title(figureInfo(2).abc); set(t,'Position',tpos,'Units','normalized')

yyaxis left; hold on

pad = range(bPlotRange)/12; 

plot([Tmin max(catS_post2(:,3))], [breference breference],'Color',c(1,:))

plot([Tms Tms],     bPlotRange,'--','Color',[0.7 0.7 0.7])    % mainshock occurrence
plot([Tms2 Tms2],   bPlotRange,'--','Color',[0.7 0.7 0.7]) % mainshock occurrence
text(Tms2+0.2/365, bPlotRange(2)-pad, sprintf('M_W%0.2g',MMS),'Color',[0.7 0.7 0.7])

text(tLim(1)+0.5/365,breference+pad,'Background',...
    'Color',c(1,:))
text(tLim(1)+0.5/365,bPlotRange(1)+pad,'Foreshocks',...
    'Color',c(2,:))
text(mean(tLim),bPlotRange(1)+pad,'Aftershocks',...
    'Color',c(3,:))



xlim(tLim)
ylim(bPlotRange)

shadedErrorBar(tpost,bpost, ...
    abs([bpostplus,bpostminus]-bpost), ...
    'lineprops',{'.-','Color',c(2,:)})

shadedErrorBar(tpost2,bpost2, ...
    abs([bpostplus2,bpostminus2]-bpost2), ...
    'lineprops',{'.-','Color',c(3,:)})

ylabel('b-value')

yyaxis right
ylim(bPlotRange/breference*100)
ylabel({'Relative b-value (%)','Traffic Light System'})

yyaxis right; hold on

N = 30;
alpha = 0.02;
for n = 0:N
   t1 = tLim(1);
   DT = diff(tLim);
   t0 = t1+4.2*DT/5;
   dt = DT/10;
   W  = tLim(2)-t0+n*dt/N;
   H  = 100*bPlotRange(2)/breference - 110;
   try
   rectangle('Position',[t0+n*dt/10, 110, W, H],...
       'FaceColor', [c(5,:), alpha], ...
       'EdgeColor', 'none');
   end
   try
   rectangle('Position',[t0+n*dt/N, 90, W, 20],...
       'FaceColor', [c(3,:), alpha], ...
       'EdgeColor', 'none');
   end
   try
   rectangle('Position',[t0+n*dt/N, 100*bPlotRange(1)/breference, W, 90-100*bPlotRange(1)/breference],...
       'FaceColor', [c(7,:), alpha], ...
       'EdgeColor', 'none');
   end
end


ax = gca;
ax.YAxis(1).Color = 'k'; 
ax.YAxis(2).Color = 'k';
xticks([])
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

subplot(figureInfo(3).subplot{:})

t = title(figureInfo(3).abc); set(t,'Position',tpos,'Units','normalized')

YLIM = [floor(min(catS(:,6))),ceil(MMS)];

yyaxis left;
decyear2dt = @(x) datetime(datenum(datetime(floor(x), 1, 1) + years(x-floor(x))),'ConvertFrom','datenum');
scatter(decyear2dt(catS(:,3)),catS(:,6),(catS(:,6)-0.5).^2, ...
    'filled', ...
    'MarkerFaceColor',[0.6 0.6 0.6],...
    'MarkerFaceAlpha',0.3); 
hold on

plot(decyear2dt([Tms Tms]),     YLIM,'--','Color',[0.7 0.7 0.7])    % mainshock occurrence
plot(decyear2dt([Tms2 Tms2]),   YLIM,'--','Color',[0.7 0.7 0.7]) % mainshock occurrence


xlim(decyear2dt(tLim))
ylim(YLIM)
ylabel('Magnitude')

ftsz    = @(fh,fontSize) set(findall(fh,'-property','FontSize'),'FontSize',fontSize);
setsize = @(fh,dim1,dim2) set(fh,...
    'Units',        'Inches', ...
    'Position',     [0,0,dim1,dim2],...
    'PaperUnits',   'Inches',...
    'PaperSize',    [dim1,dim2]);


yyaxis right; hold on
plot(decyear2dt(tpost),Mcpost,  '-', 'Color',c(2,:),'LineWidth',1.5);
plot(decyear2dt(tpost2),Mcpost2,'-', 'Color',c(3,:),'LineWidth',1.5);
ylim(YLIM)
ylabel('M_c')
yticks([])
ax = gca;
ax.YAxis(1).Color = 'k'; 
ax.YAxis(2).Color = 'k';

%ftsz(gcf,10);
setsize(gcf,7,4)

xlabel({'',[saveName, ' source volume']},'FontSize',12,'FontWeight','bold')
