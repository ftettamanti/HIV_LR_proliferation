
%% ----- Figure 2

load('Model_RUN3.mat')

for PART=11:11
    EstPars_Cell=num2cell(EstPars72Compare(PART,11:19));
    [tzero,rN,rS,rE,g,alpha,log10beta,deltaI,log10p]=EstPars_Cell{1:9};
    
    CellsVirusTIME=MODEL_RUN{PART};
    y1fitsp=y1fits(find(y1fits(:,1)==EstPars72(PART,1)),:);
    y1observationsp=y1observations(find(y1observations(:,1)==EstPars72(PART,1)),:);
    y2fitsp=y2fits(find(y2fits(:,1)==EstPars72(PART,1)),:);
    y2observationsp=y2observations(find(y2observations(:,1)==EstPars72(PART,1)),:);
    
    
    figure
    %the figure
    w=8.7;
    h=8.7*0.75;
    u='centimeters';
    pp=0.01;
    
    set(gcf,'Units',u);
    screenpos = get(gcf,'Position');
    
    set(gcf,...
        'Position',[screenpos(1:2) w h],...
        'PaperUnits',u,...
        'PaperPosition',[pp*w pp*h w h],...
        'PaperSize',[w*(1+2*pp) h*(1+2*pp)]);

    h=tight_subplot(2,2,[.1 .15],[.1 .05],[.2 .05]);
    
    axes(h(1))%subplot(2,2,1)
    hold on
    set(gca,'fontsize',6,'TickLength',[0.05, 0.01])
    %plot(CellsVirusTIME(7,:)+tzero,log10(CellsVirusTIME(6,:)/10^3),'Color','k','LineWidth',2)%,(FRH_cols(6,:)-0.5)/255,'LineWidth',2);
    plot(y2fitsp(:,2),y2fitsp(:,6),'Color','k','LineWidth',1)%,(FRH_cols(6,:)-0.5)/255,'LineWidth',2);
    plot(y2observationsp(:,2),y2observationsp(:,3),'o','MarkerSize',2,'MarkerEdgeColor', [0.75 0.75 0.75],'MarkerFaceColor', [0.75 0.75 0.75]);
    xlim([min(y2observationsp(:,2)) max(y2observationsp(:,2))])
    ylim([0 9])
    xlabel('DOFPV (days)')
    ylabel({'Viral load','(log10 copies/ml)'})
    text(-0.4,1.02,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square
    %legend({'CD4 Cells','Ln','Li'},'Location','northeast')
    %title('Virus')
    box on
    axes(h(2))%subplot(2,2,2)
    hold on
    set(gca,'fontsize',6,'TickLength',[0.05, 0.01])
    %plot(CellsVirusTIME(7,:)+tzero,(CellsVirusTIME(2,:)+CellsVirusTIME(3,:)+CellsVirusTIME(4,:)+CellsVirusTIME(5,:))/10^6,'Color','k','LineWidth',2);%,(FRH_cols(2,:)-0.5)/255,'LineWidth',2);
    plot(y1fitsp(:,2),10.^y1fitsp(:,6),'Color','k','LineWidth',1);%,(FRH_cols(2,:)-0.5)/255,'LineWidth',2);
    plot(y1observationsp(:,2),10.^y1observationsp(:,3),'o','MarkerSize',2,'MarkerEdgeColor', [0.75 0.75 0.75],'MarkerFaceColor', [0.75 0.75 0.75])
    xlim([min(y1observationsp(:,2)) max(y2observationsp(:,2))])
    ylim([0 max(10.^y1observationsp(:,3))])
    xlabel('DOFPV (days)')
    ylabel({'Total CD4 counts','(cells/\muL)'})
    text(-0.4,1.02,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square
    %axis square
    %legend({'CD4 Cells','Ln','Li'},'Location','northeast')
    box on
    axes(h(3))%subplot(2,2,3)
    plot(CellsVirusTIME(7,:)+tzero,CellsVirusTIME(1,:)/10^6,'Color','k','LineWidth',1);%(FRH_cols(1,:)-0.5)/255,'LineWidth',2);
    set(gca,'fontsize',6,'TickLength',[0.05, 0.01])
    xlim([min(y2observationsp(:,2)) max(y2observationsp(:,2))])
    ylim([0 max(10.^y1observationsp(:,3))])
    xlabel('DOFPV (days)')
    ylabel({'CD8 counts','(cells/\muL)'})
    text(-0.4,1.02,'c)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square
    %axis square
    %ylim([0 1000])
    %legend({'CD4 Cells','Ln','Li'},'Location','northeast')
    box on
    %title('CD8')
    
    
    axes(h(4))%subplot(2,2,4)
    plot(CellsVirusTIME(7,:)+tzero,CellsVirusTIME(5,:)/10^6,'Color','k','LineWidth',1);%(FRH_cols(1,:)-0.5)/255,'LineWidth',2);
    set(gca,'fontsize',6,'TickLength',[0.05, 0.01])
    set(gca, 'YScale', 'log')
    xlim([min(y2observationsp(:,2)) max(y2observationsp(:,2))])
    xlabel('DOFPV (days)')
    ylabel({'Total HIV DNA','(cells/\muL)'})
    text(-0.4,1.02,'d)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',8)    %axis square
    
    box on
    %title('CD8')
    
    %axis square
    box on
    
    %print 
    print(sprintf('fig2',PART),'-dpng','-r600')
    print(sprintf('fig2',PART),'-depsc2','-r600')
    print(sprintf('fig2',PART),'-dsvg','-r600')
    
     close all
end



%% ------ figure 3



CORR_LR=zeros(12,4);
for PART=1:12
    CellsVirusTIME=MODEL_RUN{PART};
    idx=find(CellsVirusTIME(7,:)>=30,1);
    CORR_LR(PART,1)=CellsVirusTIME(5,idx)/10^6;
    CORR_LR(PART,2)=min(CellsVirusTIME(2,:)+CellsVirusTIME(3,:)+CellsVirusTIME(4,:)+CellsVirusTIME(5,:))/10^6;
    CellsVirusTIME=MODEL_RUN2{PART};
    idx=find(CellsVirusTIME(7,:)>=30,1);
    CORR_LR(PART,3)=CellsVirusTIME(5,idx)/10^6;
    CORR_LR(PART,4)=min(CellsVirusTIME(2,:)+CellsVirusTIME(3,:)+CellsVirusTIME(4,:)+CellsVirusTIME(5,:))/10^6;
end


figure
Cols = linspecer(12);
%print the figure
w=8.7;
h=4;
u='centimeters';
pp=0.01;

set(gcf,'Units',u);
screenpos = get(gcf,'Position');

set(gcf,...
    'Position',[screenpos(1:2) w h],...
    'PaperUnits',u,...
    'PaperPosition',[pp*w pp*h w h],...
    'PaperSize',[w*(1+2*pp) h*(1+2*pp)]);
h=tight_subplot(1,2,[.1 .05],[.2 .1],[.13 .05]);
axes(h(2))
scatter(CORR_LR(:,3),CORR_LR(:,4),[],Cols,'filled')
set(gca,'fontsize',6,'TickLength',[0.05, 0.01])
ylim([50 550])
xlim([0 0.3])
ylabel('CD4 Nadir (cells/\muL)')
title('No Proliferation')
box on
%xlabel({'Latent Reservoir Size','(cells/\muL)'})
%legend({'CD4 Cells','Ln','Li'},'Location','northeast')
%title('CD4 Nadir')
%axis square
%h(2)=subplot(1,2,2);
axes(h(1))
scatter(CORR_LR(:,1),CORR_LR(:,2),[],Cols,'filled')
box on
set(gca,'fontsize',6,'TickLength',[0.05, 0.01])
title('Proliferation')
ylim([50 550])
xlim([0 0.3])
ylabel('CD4 Nadir (cells/\muL) ')
%axis square

%ylabel({'Latent Reservoir Size','(cells/\muL)'})
%legend({'CD4 Cells','Ln','Li'},'Location','northeast')
%title('CD4 Nadir')

set(h(2),'yticklabel',[])
set(h(2),'ylabel',[])
annotation('textbox',[.3 0 0.1 0.1],'String','Total HIV DNA (cells/\muL) at 30 Days Post-infection','EdgeColor','none','FontSize',6)

print('fig3','-dpng','-r600')
print('fig3','-depsc2','-r600')
print('fig3','-dsvg','-r600')

%% ----- figure 4


for PART=11:11
    EstPars_Cell=num2cell(EstPars72(PART,11:19));
    [tzero,rN,rS,rE,g,alpha,log10beta,deltaI,log10p]=EstPars_Cell{1:9};
    
    CellsVirusTIME=MODEL_RUN{PART};
    y1fitsp=y1fits(find(y1fits(:,1)==EstPars72(PART,1)),:);
    y1observationsp=y1observations(find(y1observations(:,1)==EstPars72(PART,1)),:);
    y2fitsp=y2fits(find(y2fits(:,1)==EstPars72(PART,1)),:);
    y2observationsp=y2observations(find(y2observations(:,1)==EstPars72(PART,1)),:);
    
    
    figure
    %the figure
    w=8.7;
    h=8.7*0.75;
    u='centimeters';
    pp=0.01;
    
    set(gcf,'Units',u);
    screenpos = get(gcf,'Position');
    
    set(gcf,...
        'Position',[screenpos(1:2) w h],...
        'PaperUnits',u,...
        'PaperPosition',[pp*w pp*h w h],...
        'PaperSize',[w*(1+2*pp) h*(1+2*pp)]);

    h=tight_subplot(2,2,[.1 .15],[.1 .05],[.2 .05]);
    
    axes(h(1))%subplot(2,2,1)
    hold on
    set(gca,'fontsize',6,'TickLength',[0.05, 0.01])
    %plot(CellsVirusTIME(7,:)+tzero,log10(CellsVirusTIME(6,:)/10^3),'Color','k','LineWidth',2)%,(FRH_cols(6,:)-0.5)/255,'LineWidth',2);
    plot(y2fitsp(:,2),y2fitsp(:,6),'Color','k','LineWidth',1)%,(FRH_cols(6,:)-0.5)/255,'LineWidth',2);
    plot(y2observationsp(:,2),y2observationsp(:,3),'o','MarkerSize',2,'MarkerEdgeColor', [0.75 0.75 0.75],'MarkerFaceColor', [0.75 0.75 0.75]);
    xlim([min(y2observationsp(:,2)) max(y2observationsp(:,2))])
    ylim([0 9])
    xlabel('DOFPV (days)')
    ylabel({'Viral load','(log10 copies/ml)'})
    text(-0.4,1.02,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square
    %legend({'CD4 Cells','Ln','Li'},'Location','northeast')
    %title('Virus')
    box on
    axes(h(2))%subplot(2,2,2)
    hold on
    set(gca,'fontsize',6,'TickLength',[0.05, 0.01])
    %plot(CellsVirusTIME(7,:)+tzero,(CellsVirusTIME(2,:)+CellsVirusTIME(3,:)+CellsVirusTIME(4,:)+CellsVirusTIME(5,:))/10^6,'Color','k','LineWidth',2);%,(FRH_cols(2,:)-0.5)/255,'LineWidth',2);
    plot(y1fitsp(:,2),10.^y1fitsp(:,6),'Color','k','LineWidth',1);%,(FRH_cols(2,:)-0.5)/255,'LineWidth',2);
    plot(y1observationsp(:,2),10.^y1observationsp(:,3),'o','MarkerSize',2,'MarkerEdgeColor', [0.75 0.75 0.75],'MarkerFaceColor', [0.75 0.75 0.75])
    xlim([min(y1observationsp(:,2)) max(y2observationsp(:,2))])
    ylim([0 max(10.^y1observationsp(:,3))])
    xlabel('DOFPV (days)')
    ylabel({'Total CD4 counts','(cells/\muL)'})
    text(-0.4,1.02,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square
    %axis square
    %legend({'CD4 Cells','Ln','Li'},'Location','northeast')
    box on
    axes(h(3))%subplot(2,2,3)
    plot(CellsVirusTIME(7,:)+tzero,CellsVirusTIME(1,:)/10^6,'Color','k','LineWidth',1);%(FRH_cols(1,:)-0.5)/255,'LineWidth',2);
    set(gca,'fontsize',6,'TickLength',[0.05, 0.01])
    xlim([min(y2observationsp(:,2)) max(y2observationsp(:,2))])
    ylim([0 max(10.^y1observationsp(:,3))])
    xlabel('DOFPV (days)')
    ylabel({'CD8 counts','(cells/\muL)'})
    text(-0.4,1.02,'c)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square
    %axis square
    %ylim([0 1000])
    %legend({'CD4 Cells','Ln','Li'},'Location','northeast')
    box on
    %title('CD8')
    
    
    axes(h(4))%subplot(2,2,4)
    set(gca,'fontsize',6,'TickLength',[0.05, 0.01])
    load(sprintf('LogarithmicProliferation/PART_v2_%d/RUN_%d/f%d/rL%d/sigma%d/data_%d.mat',PART,2,1,4,5,2));
    data = LR_save;
    data( :, ~any(data,1) ) = [];  % remove all zero clones (columns)
    
    idx=datasample(1:size(data,2),1e3,'Replace',false);
    %idx(end)=find(LR_save(78,:)==max(LR_save(78,:)));
    hold on
    hold on
    for i=1:size(data,2)
        plot(data(1:78,i)/10^6)
    end
    
    plot(sum(data(1:78,:)/10^6,2),'k','LineWidth',1);
    set(gca, 'YScale', 'log')
    xlim([min(y2observationsp(:,2)) max(y2observationsp(:,2))])
    xlabel('Time Post-infection (days)')
    ylabel({'Total HIV DNA','(cells/\muL)'})
    text(-0.4,1.02,'d)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',8)    %axis square
    
    %axis square
    box on
    
    %print 
    print(sprintf('fig4',PART),'-dpng','-r600')
    print(sprintf('fig4',PART),'-depsc2','-r600')
    print(sprintf('fig4',PART),'-dsvg','-r600')
    
     close all
end



