
%% determine whether plotting  foreshock or mainshock

if strcmp(saveName,'foreshock')
   figureInfo = [];
   figureInfo(1).subplot = {4,2,[1,3]};
   figureInfo(1).abc = 'a)';
   figureInfo(2).subplot = {4,2,5};
   figureInfo(2).abc = 'c)';
   figureInfo(3).subplot = {4,2,7};
   figureInfo(3).abc = 'e)';
elseif strcmp(saveName,'mainshock')
   figureInfo = [];
   figureInfo(1).subplot = {4,2,[2,4]};
   figureInfo(1).abc = 'b)';
   figureInfo(2).subplot = {4,2,6};
   figureInfo(2).abc = 'd)';
   figureInfo(3).subplot = {4,2,8};
   figureInfo(3).abc = 'f)';
else
    error('Not foreshock or mainshock')
end

%figure;
tpos = [-0.15,1]; % a-b-c- posistion
c= [0.7       0.7       0.7
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];


%subplot(2,2,[1,3])
%subplot(4,2,[1,3]);

subplot(figureInfo(1).subplot{:});
xlim([1, Mms+0.2])
ylabel('Cumulative Number')
xlabel('M')
t = title(figureInfo(1).abc); set(t,'Position',[tpos(1),tpos(2)-0.05],'Units', 'normalized')
hold on
%%
Mrange=min(catS(:,6)):mbin:max(catS(:,6));

cntnumb=hist(catSpre_orig(:,6),Mrange);
Numbh=cntnumb(end:-1:1);
Ncumh=cumsum(Numbh);
Ncumpre=Ncumh(end:-1:1);
Tcatpre=catSpre_orig(end,3)-Tmin;

plot(Mrange,Ncumpre./Tcatpre,'.','Color',c(1,:))


ll=catSpre_orig(:,6)>=Mcpre_orig;
catSpre_orig=catSpre_orig(ll,:);
if length(catSpre_orig(:,6))>=50
    bv=(1/(mean(catSpre_orig(:,6))-(Mcpre_orig-mbin/2)))*log10(exp(1));
    av_ann=log10(length(catSpre_orig(:,1))/Tcatpre)+ bv*Mcpre_orig;
    
    M2 = max(Mrange(Ncumpre~=0));
    
    N1=10^(av_ann-bv*Mcpre_orig);
    N2=10^(av_ann-bv*M2);
    lback = plot([Mcpre_orig,M2],[N1,N2],'Color',c(1,:));


    
    %Shi and Bold, 1982
    sig1 = (sum((catSpre_orig(:,6)-mean(catSpre_orig(:,6))).^2))/(length(catSpre_orig(:,6))*(length(catSpre_orig(:,6))-1));
        sig1 = sqrt(sig1);
        sig1_pre = 2.30*sig1*bv^2;
end
bv_pre = bv;
%%

Mcpost_orig2    =Mcpost_orig+2*corr_post; % kdc - should be twice the correction as it is applied twice
ll              =catS_post(:,6)>=Mcpost_orig2;
catS_post       =catS_post(ll,:);

cntnumb=hist(catS_post(:,6),Mrange);
Numbh=cntnumb(end:-1:1);
Ncumh=cumsum(Numbh);
Ncumpost=Ncumh(end:-1:1);
Tcatpost=catS_post(end,3)-catS_post(1,3);
plot(Mrange,Ncumpost./Tcatpost,'.','color',c(2,:));
hold on
if length(catS_post(:,6))>=50
    bv=(1/(mean(catS_post(:,6))-(Mcpost_orig2-mbin/2)))*log10(exp(1));
    av_ann=log10(length(catS_post(:,1))/Tcatpost)+ bv*Mcpost_orig2;
    N3=10^(av_ann-bv*Mcpost_orig2);
    N4=10^(av_ann-bv*max(Mrange));
    lfore = plot([Mcpost_orig2,max(Mrange)],[N3,N4],'color',c(2,:));
    
        %Shi and Bold, 1982
    sig1 = (sum((catS_post(:,6)-mean(catS_post(:,6))).^2))/(length(catS_post(:,6))*(length(catS_post(:,6))-1));
        sig1 = sqrt(sig1);
        sig1_post = 2.30*sig1*bv^2;

end

bv_post = bv;

TLS_post=(bv_post*100)/breference;

    

%%

ll          =catS_post2(:,3)<=tpostmax2;
catmax      =catS_post2(ll,:);

ll          =catS_post2(:,6)>=nanmedian(Mcpost2);
catS_post2  =catS_post2(ll,:);

if length(catmax(:,1))>Npost
    catmax=catmax(end:-1:end-Npost+1,:);
    catmax=sortrows(catmax,3);
end
size(catmax)

Tdiff3=tpostmax2-min(catmax(:,3));
cntnumb=hist(catmax(:,6),Mrange);

Numbh=cntnumb(end:-1:1);
Ncumh=cumsum(Numbh);
Ncum3=Ncumh(end:-1:1);

N5=10^(apostmax2-bpostmax2*magcopostmax2);
N6=10^(apostmax2-bpostmax2*max(Mrange));
lpost = plot([magcopostmax2,max(Mrange)],[N5,N6],'Color',c(3,:));
plot(Mrange,Ncum3/Tdiff3,'.','Color',c(3,:))
bpostmaxplot2=round(bpostmax2,2);

%Shi and Bold, 1982
sig1 = (sum((catmax(:,6)-mean(catmax(:,6))).^2))/(length(catmax(:,6))*(length(catmax(:,6))-1));
sig1 = sqrt(sig1);
sig1_post2 = 2.30*sig1*bpostmaxplot2^2;

 
TLS_post1=(bpostmaxplot2*100)/breference;
bv_post2 = bpostmax2;    
%%

leg = legend([lback,lfore,lpost], {'Background','Foreshocks','Aftershocks'});
leg.ItemTokenSize = [10,20];
leg.Box = 'off';
 
 xlabel('Magnitude')
 ylabel('Events per year, N(M>M_i)')
 
 set(gca,'yscale','log')
 AX = gca;
 pos = AX.Position;
 set(gca,'Position',[pos(1), pos(2)*1.1, pos(3:4)]);
 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%semilogy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
tLim        = [Tms-2/365, max(tpost2)];

% hack:
if any(bpostplus>1.3)
    bPlotRange  = [0.4 , 1.6];
else
    bPlotRange  = [0.4 , 1.25];
end 

subplot(figureInfo(2).subplot{:});

t = title(figureInfo(2).abc); set(t,'Position',tpos,'Units','normalized')

yyaxis left; hold on

pad = range(bPlotRange)/12; 

plot([Tmin max(catS_post2(:,3))], [breference breference],'Color',c(1,:))

plot([Tms Tms],     bPlotRange,'--','Color',[0.7 0.7 0.7])    % mainshock occurrence
plot([Tms2 Tms2],   bPlotRange,'--','Color',[0.7 0.7 0.7])    % mainshock occurrence
text(Tms2+0.2/365, bPlotRange(2), sprintf('M_W%0.2g',MMS),'Color',[0.7 0.7 0.7])

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

%subplot(2,2,4)
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

ftsz(gcf,10);
setsize(gcf,7,6)

xlabel({'',[saveName, ' source volume']},'FontSize',12,'FontWeight','bold')
