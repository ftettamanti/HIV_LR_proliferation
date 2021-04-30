
FRH_cols=[112 35 117; 38 57 97; 62 121 125; 98 177 182;150 194 102; 237 230 77; 242 187 68];

Markers = {'+','o','*','x','v','d','^','s','>','<','p','h'};
Markers2 = {'-+','-o','-*','-x','-v','-d','-^','-s','->','-<','-p','-h'};

Cols = linspecer(12);
%%------------- Data


samLR_cell1=cell(5,12);
LR_cell1=cell(5,12);
I_cell1=cell(5,12);
trueclones1=cell(5,12);
obsclones1=cell(5,12);
HIV_DNA_conc1=cell(5,12);
sumLR_cell1=cell(5,12);
richness_cell1=cell(5,12);
singletons_cell1=cell(5,12);
clones_cell1=cell(5,12);
sumclones_cell1=cell(5,12);
TOTCD4_cell1=cell(5,12);
times=[12,17,22,27,50,100];
idx=1;
for k = 1:12
    for i=1:5
        %    for j=1:8
        load(sprintf('LogarithmicProliferation/PART_v2_%d/RUN_%d/f%d/rL%d/sigma%d/data_%d.mat',k,i,1,4,5,i));
        data = LR_save;
        samLR_cell1{i,k} = LR_save(50,:);
        %data( :, ~any(data,1) ) = [];  % remove all zero clones (columns)
        LR_cell1{i,k} = LR_save(times+1,:);%data;
        I_cell1{i,k} = CellsVirus(4,1:1000:end)';
        HIV_DNA_conc1{i,k}=(sum(LR_save,2)+I_cell1{i,k})./(sum(CellsVirus(2:4,1:1000:end))'+sum(LR_save,2));
        trueclones1{i,k} = sum(LR_save>1,2)./sum(LR_save>0,2);
        obsclones1{i,k} = sum(LR_save>1,2)./(sum(LR_save>0,2)+I_cell1{i,k});
        sumLR_cell1{i,k}=sum(LR_save,2);
        richness_cell1{i,k}=sum(LR_save~=0,2);
        singletons_cell1{i,k}=sum(LR_save==1,2);
        clones_cell1{i,k}=sum(LR_save>1,2);
        sumclones_cell1{i,k}=zeros(size(LR_save,1),1);
        for j = 1:size(LR_save,1)
            sumclones_cell1{i,k}(j,1)=sum(LR_save(j,find(LR_save(j,:)>1)),2);
        end
        TOTCD4_cell1{i,k} = sum(CellsVirus(2:4,1:1000:end))+sumLR_cell1{i,k}';
    end
    idx=idx+1;
    %    end
end


% for k = 1:12
%     for i=1:5
%         %    for j=1:8
%         load(sprintf('LogarithmicProliferation/PART_v2_%d/RUN_%d/f%d/rL%d/sigma%d/data_%d.mat',k,i,1,4,1,i));
%         sumclones_cell1{i,k}=sum(max(LR_save,1),2);
%     end
%     idx=idx+1;
%     %    end
% end


% for k = 1:12
%     for i=1:5
%         %    for j=1:8
%         load(sprintf('LogarithmicProliferation/PART_v2_%d/RUN_%d/f%d/rL%d/sigma%d/data_%d.mat',k,i,1,4,5,i));
%         sumclones_cell1{i,k}=zeros(size(LR_save,1),1);
%         for j = 1:size(LR_save,1)
%             sumclones_cell1{i,k}(j,1)=sum(LR_save(j,find(LR_save(j,:)>1)),2);
%         end
%     end
%     idx=idx+1;
%     %    end
% end

TIMES=0:100;








%% ----- figure 5


figure

%print the figure
w=8.7;
h=11;
u='centimeters';
pp=0.01;

set(gcf,'Units',u);
screenpos = get(gcf,'Position');

set(gcf,...
    'Position',[screenpos(1:2) w h],...
    'PaperUnits',u,...
    'PaperPosition',[pp*w pp*h w h],...
    'PaperSize',[w*(1+2*pp) h*(1+2*pp)]);

%h=tight_subplot(3,2,[.05 .05],[.12 .05],[.15 .05]);


%axes(h(1:2))
subplot(4,1,1)

hold on
box on
set(gca, 'YScale', 'log','TickLength',[0.02, 0.01],'FontSize',5)
for i=[1:4 9:12 5:8]
    plot(TIMES,richness_cell1{1,i},'color',Cols(i,:),'LineWidth',1.5)
    xticks([0:25:100])
    ylim([0 1e8])
    yticks([1e2 1e4 1e6 1e8])
    ylabel({'Richness'})
end
set(gca,'XTicklabel',[])
text(-0.15,1.02,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square


subplot(4,1,2)
hold on
box on
set(gca, 'YScale', 'log','TickLength',[0.02, 0.01],'FontSize',5)
for i=[1:4 9:12 5:8]
    plot(TIMES,sumclones_cell1{1,i},'color',Cols(i,:),'LineWidth',1.5)
    ylim([0 1e8])
    yticks([1e2 1e4 1e6 1e8])
    xticks([0:25:100])
    ylabel({'Number of Clonal Cells'})
    %xlabel('Time Post-infection (days)')
end
set(gca,'XTicklabel',[])
text(-0.15,1.02,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square

subplot(4,1,3);
hold on
box on
set(gca, 'TickLength',[0.02, 0.01],'FontSize',5)
for i=[1:4 9:12 5:8]
    plot(TIMES,trueclones1{1,i},'color',Cols(i,:),'LineWidth',1.75)
    ylim([0 1])
    yticks([0:0.2:1])
    xticks([0:25:100])
    ylabel({'True Clones/Population Size'})
   % xlabel('Time Post-infection (days)')

end
set(gca,'XTicklabel',[])
text(-0.15,1.02,'c)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square
%axes(h(6))


subplot(4,1,4)
hold on
box on
set(gca, 'TickLength',[0.02, 0.01],'FontSize',5)
for i=[1:4 9:12 5:8]
    plot(TIMES,richness_cell1{1,i}./sumLR_cell1{1,i},'color',Cols(i,:),'LineWidth',1.5)
    ylim([0 1])
    yticks([0:0.2:1])
    xticks([0:25:100])
    ylabel({'Richness/Population Size'})
end
text(-0.15,1.02,'d)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square

a1=annotation('textbox',[.35 0 0.1 0.1],'String','Time Post-infection (days)','EdgeColor','none','FontSize',7);



%tightfig;

print('fig5','-dpng','-r600')
print('fig5','-depsc2','-r600')
print('fig5','-dsvg','-r600')

close all






%% ----- figure 6
figure

%print the figure
w=8.7;
h=8;
u='centimeters';
pp=0.01;

set(gcf,'Units',u);
screenpos = get(gcf,'Position');

set(gcf,...
    'Position',[screenpos(1:2) w h],...
    'PaperUnits',u,...
    'PaperPosition',[pp*w pp*h w h],...
    'PaperSize',[w*(1+2*pp) h*(1+2*pp)]);

%h=tight_subplot(3,2,[.05 .05],[.12 .05],[.15 .05]);


%axes(h(1:2))

subplot(3,1,1)
hold on
box on
set(gca, 'YScale', 'log','TickLength',[0.02, 0.01],'FontSize',5)
for i=[1:4 9:12 5:8]
    plot(TIMES,richness_cell1{1,i}+I_cell1{1,i},'color',Cols(i,:),'LineWidth',1.5)
    ylim([0 1e8])
    yticks([1e2 1e4 1e6 1e8])
    xticks([0:25:100])
    ylabel({'Richness'})    
end
set(gca,'XTicklabel',[])
text(-0.15,1.02,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square

subplot(3,1,2);
hold on
box on
set(gca, 'TickLength',[0.02, 0.01],'FontSize',5)
for i=[1:4 9:12 5:8]
    plot(TIMES,obsclones1{1,i},'color',Cols(i,:),'LineWidth',1.75)
    ylim([0 1])
    yticks([0:0.2:1])
    xticks([0:25:100])
    ylabel({'True Clones/Population Size'})    
    %xlabel('Time Post-infection (days)')

end
set(gca,'XTicklabel',[])
text(-0.15,1.02,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square

%text(-0.15,1.02,'d)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square
subplot(3,1,3)
hold on
box on
set(gca, 'TickLength',[0.02, 0.01],'FontSize',5)
for i=[1:4 9:12 5:8]
    plot(TIMES,(richness_cell1{1,i}+I_cell1{1,i})./(sumLR_cell1{1,i}+I_cell1{1,i}),'color',Cols(i,:),'LineWidth',1.5)
    ylim([0 1])
    yticks([0:0.2:1])
    xticks([0:25:100])
    ylabel({'Richness/Population Size'})
end
text(-0.15,1.02,'c)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square
a1=annotation('textbox',[.35 0 0.1 0.1],'String','Time Post-infection (days)','EdgeColor','none','FontSize',7);



%tightfig;

print('fig6','-dpng','-r600')
print('fig6','-depsc2','-r600')
print('fig6','-dsvg','-r600')

close all




%% ----- figure 7


figure

%print the figure
w=8.7;
h=8.7*1.5;
u='centimeters';
pp=0.01;

set(gcf,'Units',u);
screenpos = get(gcf,'Position');

set(gcf,...
    'Position',[screenpos(1:2) w h],...
    'PaperUnits',u,...
    'PaperPosition',[pp*w pp*h w h],...
    'PaperSize',[w*(1+2*pp) h*(1+2*pp)]);

h=tight_subplot(3,2,[.05 .05],[.12 .05],[.15 .05]);


axes(h(1))
hold on
box on
set(gca, 'YScale', 'log','TickLength',[0.05, 0.01],'FontSize',7)
for i=[1:4 9:12 5:8]
    plot(TIMES,singletons_cell1{1,i},'color',Cols(i,:),'LineWidth',1.5)
    xticks([0:25:100])
    ylim([0 1e8])
    yticks([1e2 1e4 1e6 1e8])
    ylabel({'Number of Singletons (per L)'})
    title('Latently Infected Cells')
end
axes(h(2))
hold on
box on
set(gca, 'YScale', 'log','TickLength',[0.05, 0.01],'FontSize',7)
for i=[1:4 9:12 5:8]
    plot(TIMES,singletons_cell1{1,i}+I_cell1{1,i},'color',Cols(i,:),'LineWidth',1.5)
    ylim([0 1e8])
    yticks([1e2 1e4 1e6 1e8])
    xticks([0:25:100])
    ylabel({'Number of Singletons (per L)'})
    title('All Infected Cells')
end
axes(h(3))
hold on
box on
set(gca, 'YScale', 'log','TickLength',[0.05, 0.01],'FontSize',7)
for i=[1:4 9:12 5:8]
    plot(TIMES,clones_cell1{1,i},'color',Cols(i,:),'LineWidth',1.5)
    ylim([0 1e8])
    yticks([1e2 1e4 1e6 1e8])
    xticks([0:25:100])
    ylabel({'Number of True Clones (per L)'})
end
% axes(h(4))
% hold on
% box on
% set(gca, 'TickLength',[0.05, 0.01],'FontSize',8)
% for i=[1:4 9:12 5:8]
%     plot(TIMES,(richness_cell1{1,i}+I_cell1{1,i})/(sumLR_cell1{1,i}+I_cell1{1,i}),'color',Cols(i,:),'LineWidth',1.5)
%     ylim([0 1])
%     yticks([0:0.2:1])
%     xticks([0:25:100])
%     ylabel({'Richness/Population Size'})
% end
axes(h(5))
hold on
box on
set(gca, 'TickLength',[0.05, 0.01],'FontSize',7)
for i=[1:4 9:12 5:8]
    plot(TIMES,clones_cell1{1,i}./richness_cell1{1,i},'color',Cols(i,:),'LineWidth',1.75)
    ylim([0 1])
    yticks([0:0.2:1])
    xticks([0:25:100])
    ylabel({'True Clones/Richness'})
end
axes(h(6))
hold on
box on
set(gca, 'TickLength',[0.05, 0.01],'FontSize',7)
for i=[1:4 9:12 5:8]
    plot(TIMES,clones_cell1{1,i}./(richness_cell1{1,i}+I_cell1{1,i}),'color',Cols(i,:),'LineWidth',1.75)
    ylim([0 1])
    yticks([0:0.2:1])
    xticks([0:25:100])
    ylabel({'True Clones/Richness'})
end
set(h(1:4),'xticklabel',[])
set(h([2 4 6]),'yticklabel',[])
set(h([2 4 6]),'ylabel',[])

delete(h(4));

a1=annotation('textbox',[.35 0 0.1 0.1],'String','Time Post-infection (days)','EdgeColor','none','FontSize',7);



%tightfig;

print('fig16c_v3','-dpng','-r600')
print('fig16c_v3','-depsc2','-r600')
print('fig16c_v3','-dsvg','-r600')


%%----------- Singletons and clones



figure

%print the figure
w=8.7;
h=10;
u='centimeters';
pp=0.01;

set(gcf,'Units',u);
screenpos = get(gcf,'Position');

set(gcf,...
    'Position',[screenpos(1:2) w h],...
    'PaperUnits',u,...
    'PaperPosition',[pp*w pp*h w h],...
    'PaperSize',[w*(1+2*pp) h*(1+2*pp)]);

%h=tight_subplot(3,2,[.05 .05],[.12 .05],[.15 .05]);


%axes(h(1:2))
subplot(3,1,1)
hold on
box on
set(gca, 'YScale', 'log','TickLength',[0.05, 0.01],'FontSize',5)
for i=[1:4 9:12 5:8]
    plot(TIMES,singletons_cell1{1,i},'color',Cols(i,:),'LineWidth',1.5)
    xticks([0:25:100])
    ylim([0 1e8])
    yticks([1e2 1e4 1e6 1e8])
    ylabel({'Number of Singleton Cells'})
    xlabel('Time Post-infection (days)')
end
text(-0.15,1.02,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square

%axes(h(3:4))
subplot(3,1,2)
hold on
box on
set(gca, 'YScale', 'log','TickLength',[0.05, 0.01],'FontSize',5)
for i=[1:4 9:12 5:8]
    plot(TIMES,sumclones_cell1{1,i},'color',Cols(i,:),'LineWidth',1.5)
    ylim([0 1e8])
    yticks([1e2 1e4 1e6 1e8])
    xticks([0:25:100])
    ylabel({'Number of Clonal Cells'})
    xlabel('Time Post-infection (days)')
end
text(-0.15,1.02,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square
% axes(h(4))
% hold on
% box on
% set(gca, 'TickLength',[0.05, 0.01],'FontSize',8)
% for i=[1:4 9:12 5:8]
%     plot(TIMES,(richness_cell1{1,i}+I_cell1{1,i})/(sumLR_cell1{1,i}+I_cell1{1,i}),'color',Cols(i,:),'LineWidth',1.5)
%     ylim([0 1])
%     yticks([0:0.2:1])
%     xticks([0:25:100])
%     ylabel({'Richness/Population Size'})
% end
%axes(h(5))
h(1)=subplot(3,2,5);
hold on
box on
set(gca, 'TickLength',[0.05, 0.01],'FontSize',5)
for i=[1:4 9:12 5:8]
    plot(TIMES,trueclones1{1,i},'color',Cols(i,:),'LineWidth',1.75)
    ylim([0 1])
    yticks([0:0.2:1])
    xticks([0:25:100])
    ylabel({'True Clones/Population Size'})
        title('Latently Infected Cells')
   % xlabel('Time Post-infection (days)')

end
text(-0.32,1.02,'c)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square
%axes(h(6))
h(2)=subplot(3,2,6);
hold on
box on
set(gca, 'TickLength',[0.05, 0.01],'FontSize',5)
for i=[1:4 9:12 5:8]
    plot(TIMES,obsclones1{1,i},'color',Cols(i,:),'LineWidth',1.75)
    ylim([0 1])
    yticks([0:0.2:1])
    xticks([0:25:100])
    ylabel({'True Clones/Population Size'})    
    title('All Infected Cells')
    %xlabel('Time Post-infection (days)')

end
%text(-0.15,1.02,'d)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square
set(h([2]),'yticklabel',[])
set(h([2]),'ylabel',[])
a1=annotation('textbox',[.35 0 0.1 0.1],'String','Time Post-infection (days)','EdgeColor','none','FontSize',6);

print('fig16c_v4','-dpng','-r600')
print('fig16c_v4','-depsc2','-r600')
print('fig16c_v4','-dsvg','-r600')



for i=[1:4 9:12 5:8]
    fprintf('%d\n',trueclones1{1,i}(101))
end
for i=[1:4 9:12 5:8]
    fprintf('%0.4g\n',max(richness_cell1{1,i}./sumLR_cell1{1,i}))
end
for i=[1:4 9:12 5:8]
    fprintf('%d\n',TIMES(find(richness_cell1{1,i}./sumLR_cell1{1,i}<0.5,1)))
end
for i=[1:4 9:12 5:8]
    fprintf('%d\n',TIMES(find(richness_cell1{1,i}./sumLR_cell1{1,i}<0.1,1)))
end
for i=[1:4 9:12 5:8]
    fprintf('%d\n',obsclones1{1,i}(11))
end
for i=[1:4 9:12 5:8]
    fprintf('%d\n',obsclones1{1,i}(51))
end
for i=[1:4 9:12 5:8]
    argh=(richness_cell1{1,i}+I_cell1{1,i})./(sumLR_cell1{1,i}+I_cell1{1,i});
    fprintf('%d\n',argh(101))
end
