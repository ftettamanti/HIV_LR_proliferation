% 7_p1 and 7_p2 were glued together in inkscape

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

times=[12,17,22,27,50,100];
idx=1;
for k = 1:12
    for i=1:5
        %    for j=1:8
        load(sprintf('LogarithmicProliferation/PART_v2_%d/RUN_%d/f%d/rL%d/sigma%d/data_%d.mat',k,i,1,4,1,i));
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
    end
    idx=idx+1;
    %    end
end

TIMES=0:100;




%%----------- Integration Sites Sampling + Infected Cells

% Fiebig stages


Positive_10=cell(5,12);
Positive_100=cell(5,12);
Positive_1000=cell(5,12);

ylabs={{'Fiebig I','12 days'},{'Fiebig II','17 days'},{'Fiebig III','22 days'}...
    ,{'Fiebig IV','27 days'},{'Fiebig V','50 days'},{'Fiebig VI','100 days'}};

figure

%print the figure
w=8.7/2;
h=8.7*2;
u='centimeters';
pp=0.01;

set(gcf,'Units',u);
screenpos = get(gcf,'Position');

set(gcf,...
    'Position',[screenpos(1:2) w h],...
    'PaperUnits',u,...
    'PaperPosition',[pp*w pp*h w h],...
    'PaperSize',[w*(1+2*pp) h*(1+2*pp)]);

%hh=tight_subplot(6,2,[.05 .05],[.12 .05],[.2 .05]);

%tightfig;

for t = 1:6
    %axes(hh(2*t))
    hh(t) = subplot(6,1,t);
    set(gca,'TickLength',[0.05, 0.01],'FontSize',5)
    for k=1:12
        for j=1:5
            %fprintf('%d %d \n',k,j);
            data = LR_cell1{j,k}(t,:);%samLR_cell{j,k};
            data(find(data==0))=[];
            Pops=zeros(1,sum(data));
            indx=0;
            for i=1:size(data,2)
                Pops(indx+1:(indx+data(i)))=i;
                indx=indx+data(i);
            end
            Pops=[Pops (size(Pops,2)+1):(size(Pops,2)+I_cell1{j,k}(times(t)+1))];
            
            Positive_10{j,k}=0;
            Positive_100{j,k}=0;
            Positive_1000{j,k}=0;
            
            data_sample=datasample(Pops,1e1,'Replace',false);
            tally_1=tabulate(data_sample);
            Positive_10{j,k}=sum(tally_1(:,2)>1)/sum(tally_1(:,2));
            %
            data_sample=datasample(Pops,1e2,'Replace',false);
            tally_1=tabulate(data_sample);
            Positive_100{j,k}=sum(tally_1(:,2)>1)/sum(tally_1(:,2));
            %
            data_sample=datasample(Pops,1e3,'Replace',false);
            tally_1=tabulate(data_sample);
            Positive_1000{j,k}=sum(tally_1(:,2)>1)/sum(tally_1(:,2));
        end
    end
    
    
    Markers = {'+','o','*','x','v','d','^','s','>','<','p','h'};
    dat=[reshape(cell2mat(Positive_10),[],1) ...
        reshape(cell2mat(Positive_100),[],1) ...
        reshape(cell2mat(Positive_1000),[],1)];
    %dat=dat(:,[1 4 7 2 5 8 3 6 9]);
    fprintf("%0.5g %0.5g %0.5g \n",median(reshape(cell2mat(Positive_10)*100,[],1)),...
        median(reshape(cell2mat(Positive_100)*100,[],1)),median(reshape(cell2mat(Positive_1000)*100,[],1)));
    dat=dat*100;
    positions=[1:3];
    h = boxplot(dat);
    % set (gca,'yscale','log');
    %bp.XAxis.TickLabelInterpreter = 'latex';
    % axis square
    hold on
    a = get(get(gca,'children'),'children');   % Get the handles of all the objects
    set(h,{'linew'},{0.5}) % box and whiskers line width
    set(gca,'xtick',[positions(1) positions(2) positions(3)])
    set(gca,'xticklabel',{'10','100','1000'})
    set(a, 'Color', 'k');   % Set the color of the first box to black
    x1=ones(length(dat)).*(positions(1)+(rand(length(dat))-0.5)/3);
    x2=ones(length(dat)).*(positions(2)+(rand(length(dat))-0.5)/3);
    x3=ones(length(dat)).*(positions(3)+(rand(length(dat))-0.5)/3);
    % x6=ones(length(dat)).*(positions(6)+(rand(length(dat))-0.5)/3);
    % x7=ones(length(dat)).*(positions(7)+(rand(length(dat))-0.5)/3);
    % x8=ones(length(dat)).*(positions(8)+(rand(length(dat))-0.5)/3);
    % x9=ones(length(dat)).*(positions(9)+(rand(length(dat))-0.5)/3);
    for i=1:5
        for j=1:12
            f1=scatter(x1(5*(j-1)+i,1),dat(5*(j-1)+i,1),6,'k',Markers{i});f1.MarkerEdgeColor = Cols(j,:);f1.MarkerFaceColor =Cols(j,:);f1.MarkerFaceAlpha = 0.4;hold on
            f2=scatter(x2(5*(j-1)+i,1),dat(5*(j-1)+i,2),6,'k',Markers{i});f2.MarkerEdgeColor = Cols(j,:);f1.MarkerFaceColor =Cols(j,:);f1.MarkerFaceAlpha = 0.4;hold on
            f3=scatter(x3(5*(j-1)+i,1),dat(5*(j-1)+i,3),6,'k',Markers{i});f3.MarkerEdgeColor = Cols(j,:);f1.MarkerFaceColor =Cols(j,:);f1.MarkerFaceAlpha = 0.4;hold on
        end
    end
    
    set(gca,'xticklabel',{[]},'FontSize',5)

    
%     yname=ylabs(t);
%     ylabel(yname{1});
    %xlabel('Number of Infected Cells Sequenced');
    ylim([0 30])
    yticks([0:10:30])
    % axis square
    
end

set(hh(1:5),'xticklabel',[])
set(hh(1:5),'xlabel',[])
set(hh(1:5),'TickLength',[0.05, 0.01])


print('fig7_p1','-dpng','-r600')
print('fig7_p1','-depsc2','-r600')
print('fig7_p1','-dsvg','-r600')


figure

%print the figure
w=8.7/2;
h=8.7*2;
u='centimeters';
pp=0.01;

set(gcf,'Units',u);
screenpos = get(gcf,'Position');

set(gcf,...
    'Position',[screenpos(1:2) w h],...
    'PaperUnits',u,...
    'PaperPosition',[pp*w pp*h w h],...
    'PaperSize',[w*(1+2*pp) h*(1+2*pp)]);

%hh=tight_subplot(6,2,[.05 .05],[.12 .05],[.2 .05]);

%tightfig;


Positive_10=cell(5,12);
Positive_100=cell(5,12);
Positive_1000=cell(5,12);

for t = 1:6
    %axes(hh(2*t-1))
    
    hh(t) = subplot(6,1,t);
    set(gca,'TickLength',[0.05, 0.01],'FontSize',5)
    for k=1:12
        for j=1:5
            %fprintf('%d %d \n',k,j);
            data = LR_cell1{j,k}(t,:);%samLR_cell{j,k};
            data(find(data==0))=[];
            Pops=zeros(1,sum(data));
            indx=0;
            for i=1:size(data,2)
                Pops(indx+1:(indx+data(i)))=i;
                indx=indx+data(i);
            end
            
            Positive_10{j,k}=0;
            Positive_100{j,k}=0;
            Positive_1000{j,k}=0;
            
            data = Pops;
            if (size(data,2)>=10^1)
                data_sample=datasample(data,10^1,'Replace',false);
                data_sample(find(data_sample==0))=[];
                tally_1=tabulate(data_sample);
                %tally_1=tabulate(data_sample);
                Positive_10{j,k}=sum(tally_1(:,2)>1)/sum(tally_1(:,2));%sum(tally_1(:,2)>1)/sum(tally_1(:,2));%data_sample;%sum(data_sample>1)/sum(data_sample);
            else
                
                Positive_10{j,k}= NaN;
            end
            if (size(data,2)>=10^2)
                data_sample=datasample(data,10^2,'Replace',false);
                data_sample(find(data_sample==0))=[];
                tally_1=tabulate(data_sample);
                %tally_1=tabulate(data_sample);
                Positive_100{j,k}=sum(tally_1(:,2)>1)/sum(tally_1(:,2));%sum(tally_1(:,2)>1)/sum(tally_1(:,2));%data_sample;%sum(data_sample>1)/sum(data_sample);
            else
                Positive_100{j,k}=NaN;
            end
            if (size(data,2)>=10^3)
                data_sample=datasample(data,10^3,'Replace',false);
                data_sample(find(data_sample==0))=[];
                tally_1=tabulate(data_sample);
                %tally_1=tabulate(data_sample);
                Positive_1000{j,k}=sum(tally_1(:,2)>1)/sum(tally_1(:,2));%sum(tally_1(:,2)>1)/sum(tally_1(:,2));%data_sample;%sum(data_sample>1)/sum(data_sample);
            else
                Positive_1000{j,k}= NaN;
            end
        end
    end
    
    
    Markers = {'+','o','*','x','v','d','^','s','>','<','p','h'};
    dat=[reshape(cell2mat(Positive_10),[],1) ...
        reshape(cell2mat(Positive_100),[],1) ...
        reshape(cell2mat(Positive_1000),[],1)];
    dummy=reshape(cell2mat(Positive_1000)*100,[],1);
    dummy(isnan(dummy))=[];
    fprintf("%0.5g %0.5g %0.5g \n",median(reshape(cell2mat(Positive_10)*100,[],1)),...
        median(reshape(cell2mat(Positive_100)*100,[],1)),median(dummy));
    %dat=dat(:,[1 4 7 2 5 8 3 6 9]);
    dat=dat*100;
    positions=[1:3];
    h = boxplot(dat);
    % set (gca,'yscale','log');
    %bp.XAxis.TickLabelInterpreter = 'latex';
    % axis square
    hold on
    a = get(get(gca,'children'),'children');   % Get the handles of all the objects
    set(h,{'linew'},{0.5}) % box and whiskers line width
    set(gca,'xtick',[positions(1) positions(2) positions(3)])
    set(gca,'xticklabel',{'10','100','1000'})
    set(a, 'Color', 'k');   % Set the color of the first box to black
    x1=ones(length(dat)).*(positions(1)+(rand(length(dat))-0.5)/3);
    x2=ones(length(dat)).*(positions(2)+(rand(length(dat))-0.5)/3);
    x3=ones(length(dat)).*(positions(3)+(rand(length(dat))-0.5)/3);
    % x6=ones(length(dat)).*(positions(6)+(rand(length(dat))-0.5)/3);
    % x7=ones(length(dat)).*(positions(7)+(rand(length(dat))-0.5)/3);
    % x8=ones(length(dat)).*(positions(8)+(rand(length(dat))-0.5)/3);
    % x9=ones(length(dat)).*(positions(9)+(rand(length(dat))-0.5)/3);
    for i=1:5
        for j=1:12
            f1=scatter(x1(5*(j-1)+i,1),dat(5*(j-1)+i,1),6,'k',Markers{i});f1.MarkerEdgeColor = Cols(j,:);f1.MarkerFaceColor =Cols(j,:);f1.MarkerFaceAlpha = 0.4;hold on
            f2=scatter(x2(5*(j-1)+i,1),dat(5*(j-1)+i,2),6,'k',Markers{i});f2.MarkerEdgeColor = Cols(j,:);f1.MarkerFaceColor =Cols(j,:);f1.MarkerFaceAlpha = 0.4;hold on
            f3=scatter(x3(5*(j-1)+i,1),dat(5*(j-1)+i,3),6,'k',Markers{i});f3.MarkerEdgeColor = Cols(j,:);f1.MarkerFaceColor =Cols(j,:);f1.MarkerFaceAlpha = 0.4;hold on
        end
    end
%     yname=ylabs(t);
%     ylabel(yname{1});
    %xlabel('Number of Infected Cells Sequenced');
    ylim([0 30])
    yticks([0:10:30])
   set(gca,'FontSize',5)

    % axis square
        %set(gca,'xticklabel',{[]})
end
set(hh(1:6),'xticklabel',[])
set(hh(1:6),'xlabel',[])
set(hh(1:6),'TickLength',[0.05, 0.01])
% set(hh([1:10]),'xticklabel',[])
% set(hh([2 4 6 8 10 12]),'yticklabel',[])
% set(hh([2 4 6 8 10 12]),'ylabel',[])
% set(hh([1:2:12]),'ylabel','FontSize')
% set(hh([1:12]),'xlabel',[])
% a1=annotation('textbox',[.35 0 0.1 0.1],'String','Number of Cells Sampled','EdgeColor','none','FontSize',7);


print('fig7_p2','-dpng','-r600')
print('fig7_p2','-depsc2','-r600')
print('fig7_p2','-dsvg','-r600')

