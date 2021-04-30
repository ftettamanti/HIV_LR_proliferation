
FRH_cols=[112 35 117; 38 57 97; 62 121 125; 98 177 182;150 194 102; 237 230 77; 242 187 68];

Markers = {'+','o','*','x','v','d','^','s','>','<','p','h'};
Markers2 = {'-+','-o','-*','-x','-v','-d','-^','-s','->','-<','-p','-h'};

Cols = linspecer(12);

times=[12,17,22,27,50,463];
%%------------- Data

load(sprintf('LogarithmicProliferation/PART_v3_%d/RUN_%d/f%d/rL%d/sigma%d/data_%d.mat',11,2,1,4,5,2));

figure
for t=1:6
    subplot(2,3,t)
    data=-sort(-LR_save(times(t)+1,:)); %make sure correctly ranked
    data_pa  = data/sum(data);
    data_cpa = cumsum(data_pa);
    data_r   = 1:length(data);
    %[f,x]=ecdf(LR_cell1{i,k});
    %loglog((x),(f),'.-')
    loglog((data_r),100*(data_pa),'.-','Color','k','linewidth',2,'MarkerSize',5);
    ylim([1e-6 1e2])
    xlim([1 1e4])
    ylabel('Proportional Abundance (%)')
    xlabel('Rank')
    %title_txt = sprintf('$$r_L=%0.2f \\times r_N, \\sigma  = %0.2f $$',RN2(4),Sigma(tester(k)));
    %title(title_txt,'interpreter','latex')
    %axis square
    
    
end

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

subplot(2,5,3:4)
hold on
box on
axis square
set(gca,'fontsize',6)
for t=1:6
    fprintf('%d \n ',CellsVirus(4,1+1000*times(t)+1))
    data=-sort(-[LR_save(times(t)+1,:) ones(1,round(CellsVirus(4,1+1000*times(t)+1)))]); %make sure correctly ranked
    data_pa  = data/sum(data);
    data_cpa = cumsum(data_pa);
    data_r   = 1:length(data);
    %[f,x]=ecdf(LR_cell1{i,k});
    %loglog((x),(f),'.-')
    h(t)=plot((data_r),100*(data_pa),'.-','linewidth',1,'MarkerSize',5);
    
    ylim([1e-6 1e2])
    xlim([1 1e8])
    xticks([1 1e2 1e4 1e6 1e8])
    yticks([1e-6 1e-4 1e-2 1e0 1e2])
    ylabel('Proportional Abundance (%)')
    xlabel('Rank')
    %title_txt = sprintf('$$r_L=%0.2f \\times r_N, \\sigma  = %0.2f $$',RN2(4),Sigma(tester(k)));
    %title(title_txt,'interpreter','latex')
    %axis square
    set (gca,'yscale','log');
    set (gca,'xscale','log');
end
title({'Population of  ',' Infected Cell'})
subplot(2,5,8:9)
hold on
box on
axis square
set(gca,'fontsize',6)
for t=1:6
    data=-sort(-[LR_save(times(t)+1,:) ones(1,round(CellsVirus(4,1+1000*times(t)+1)))]); %make sure correctly ranked
    population=zeros(1,sum(data));
    nums=0;
    for l=1:size(data,2)
        %     fprintf('%d\n',nums)
        %     fprintf('%d\n',data(i))
        population((nums+1):(nums+data(l)))=l;
        nums=nums+data(l);
    end
    
    data_sample=datasample(population,2e2,'Replace',false);
    tally_1=tabulate(data_sample);
    % tally_1=array2table(tally_1,'VariableNames',{'Value','Count','Percent'});
    
    tally_1(find(tally_1(:,2)==0),:)=[];
    %fprintf('%0.2f \n',size(find(tally_1(:,2)==1),1)/5e2);
    data=-sort(-tally_1(:,2)); %make sure correctly ranked
    data_pa  = data/sum(data);
    data_cpa = cumsum(data_pa);
    data_r   = 1:length(data);
    %[f,x]=ecdf(LR_cell1{i,k});
    %loglog((x),(f),'.-')
    h(2*t) = plot((data_r),100*(data_pa),'.-','linewidth',1,'MarkerSize',5);%scatter((data_r),100*(data_pa),'o','filled','MarkerEdgeColor','flat','MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.2);
    ylim([5e-1 1e2])
    xlim([1 2e2])
    xticks([1 1e1 1e2])
    %      ylim([5e-8 1e2])
    %      xlim([1 1e8])
    ylabel('Proportional Abundance (%)')
    xlabel('Rank')
    %title_txt = sprintf('$$r_L=%0.2f \\times r_N, \\sigma  = %0.2f $$',RN2(4),Sigma(tester(k)));
    %title(title_txt,'interpreter','latex')
    %axis square
    set (gca,'yscale','log');
    set (gca,'xscale','log');
end
title({'Sample of','Infected Cells'})
subplot(2,5,1:2)
hold on
box on
axis square
set(gca,'fontsize',6)
for t=1:6
    data=-sort(-LR_save(times(t)+1,:)); %make sure correctly ranked
    data_pa  = data/sum(data);
    data_cpa = cumsum(data_pa);
    data_r   = 1:length(data);
    %[f,x]=ecdf(LR_cell1{i,k});
    %loglog((x),(f),'.-')
    plot((data_r),100*(data_pa),'.-','linewidth',1,'MarkerSize',5);
    ylim([1e-6 1e2])
    xlim([1 1e8])
    xticks([1 1e2 1e4 1e6 1e8])
    yticks([1e-6 1e-4 1e-2 1e0 1e2])
    
    ylabel('Proportional Abundance (%)')
    xlabel('Rank')
    %title_txt = sprintf('$$r_L=%0.2f \\times r_N, \\sigma  = %0.2f $$',RN2(4),Sigma(tester(k)));
    %title(title_txt,'interpreter','latex')
    %axis square
    set (gca,'yscale','log');
    set (gca,'xscale','log');
end
text(-0.3,1.02,'a)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square
title({'Population of ','Latently Infected Cells'})
subplot(2,5,6:7)
hold on
box on
axis square
set(gca,'fontsize',6)
for t=1:6
    data=-sort(-LR_save(times(t)+1,:)); %make sure correctly ranked
    population=zeros(1,sum(data));
    nums=0;
    for l=1:size(data,2)
        %     fprintf('%d\n',nums)
        %     fprintf('%d\n',data(i))
        population((nums+1):(nums+data(l)))=l;
        nums=nums+data(l);
    end
    
    data_sample=datasample(population,2e2,'Replace',false);
    tally_1=tabulate(data_sample);
    % tally_1=array2table(tally_1,'VariableNames',{'Value','Count','Percent'});
    
    tally_1(find(tally_1(:,2)==0),:)=[];
    %fprintf('%0.2f \n',size(find(tally_1(:,2)==1),1)/5e2);
    data=-sort(-tally_1(:,2)); %make sure correctly ranked
    data_pa  = data/sum(data);
    data_cpa = cumsum(data_pa);
    data_r   = 1:length(data);
    %[f,x]=ecdf(LR_cell1{i,k});
    %loglog((x),(f),'.-')
    plot((data_r),100*(data_pa),'.-','linewidth',1,'MarkerSize',5);%scatter((data_r),100*(data_pa),'o','filled','MarkerEdgeColor','flat','MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.2);
    
    ylim([5e-1 1e2])
    xlim([1 2e2])
    xticks([1 1e1 1e2])
    ylabel('Proportional Abundance (%)')
    xlabel('Rank')
    %title_txt = sprintf('$$r_L=%0.2f \\times r_N, \\sigma  = %0.2f $$',RN2(4),Sigma(tester(k)));
    %title(title_txt,'interpreter','latex')
    %axis square
    set (gca,'yscale','log');
    set (gca,'xscale','log');
end
title({'Sample of ',' Latently Infected Cells'})
text(-0.3,1.02,'b)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontWeight','bold','FontSize',6)    %axis square



% hLegend = subplot(2,5,10);
% posLegend = get(hLegend,'Position');
% leg=legend(hLegend,h,'FS I','FS II','FS III','FS IV','FS V','ART 1 YR');
% axis(hLegend,'off');
% set(leg,'Position',posLegend);


%tightfig;

print('fig8','-dpng','-r600')
print('fig8','-depsc2','-r600')
print('fig8','-dsvg','-r600')