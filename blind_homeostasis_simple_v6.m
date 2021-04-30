%***********************************************************************
%
% This program runs the lymphopaenia stochastic model with interpolated VL
% and CD4 from participant 12330251186 from Fresh Cohort
% with tau-leaping
% 08/04/2019 - FATB
%
% The deterministic model considered is:
% \begin{align}
% \frac{\text{d} I }{\text{d} t} &= g \beta V(t) T(t) - d_I I \\
% \frac{\text{d} L }{\text{d} t} &= f \beta V(t) T(t) + r\left( 1- T/K \right) L
% \end{align}
%
% This is application of a gillespie algorithm to simulate a birth/death
% process of different clones in the reservoir during fiebig stage 1 & 2.
% For each clone, two things
% can happen:
% 1 - the creation of a cell - with logistic growth
% 2 - the destruction of a cell
% 3 - the introduction of a new cell type - 2 introduced at start point
% %
%***********************************************************************

function blind_homeostasis_simple_v6(f,r,sigma,VIRUS,CD4,times)
%
%>>>>>>>>> parameter values from monolix fit & literature <<<<<<<<<<%
%
rC = 0.06;
beta = 10^(-4)*10^(-6);
dI = 1;
p = 10^4;
c = 23;
K = (2.8);%*10^(6);
% f = 1e-4;
a = 5.7e-5;
g = (1-f);%0.05;
instep = 1e-2;

%
%>>>>>>>>> Initialize infected cells and LR
%
% >> Infected cells
Inf = zeros(1,size(times,2));
Inf(1)=0;
% >> LR
n=1;
aa=zeros(1,2*n+1);
LR = zeros(2,10^3);
sumLR = sum(LR(1,:),2);
%rN2=0.015;
rL = zeros(1,10^6);
if r>0
    if sigma<0
        rL = zeros(1,10^6);
    else
        rL=r*exp(normrnd(0,sigma*r,1,10^6));%zeros(1,10^5);%zeros(1,10^5);%rN*ones(1,10^5);%
    end
end
LR_save = zeros(size(unique(round(times)),2),10^6);
LR_save(1,1:10^3)=LR(2,:);

%
% %>>>>>>>>> OPEN FILES
% %
% % fileID1 = fopen(sprintf('blind_LR.txt'),'w');
% % % fprintf(fileID, 'Results nn \n');
% % fprintf(fileID1,'%d ',LR(1,:));
% % fprintf(fileID1,'\n');
% 
% fileID2 = fopen(sprintf('blind_I.txt'),'w');
% % fprintf(fileID, 'Results nn \n');
% fprintf(fileID2,'%d ',Inf(1));
% fprintf(fileID2,'\t');
% 
% fileID3 = fopen(sprintf('blind_t.txt'),'w');
% % fprintf(fileID, 'Results nn \n');
% fprintf(fileID3,'%d ',times(1));
% fprintf(fileID3,'\t');
% 
% %  fileID4 = fopen(sprintf('blind_LR.txt'),'w');
% % % fprintf(fileID, 'Results nn \n');
% % fwrite(fileID4,[0 0]');
% %fprintf(fileID4,'%d\t',[0 0]');
% %fprintf(fileID4,'\n');
% % fprintf(fileID4,'%d ',LR(1,:));
% % fprintf(fileID4,'\t');

%
%>>>>>>>>> Rescale K
%
K=10^K*10^6;
CD4 = 10^6*CD4;
%
%>>>>>>>>> Initialise time index of saving
%
% CD4=CD4(find(times>=0,1):end);
% VIRUS=VIRUS(find(times>=0,1):end);
% times=times(find(times>=0,1):end);
time_index=times(1)+1;
save_index=2;
for i=2:size(times,2)
    % Infected cells
    [t,C] = ode45(@(t, C) Infected(t,C,K,r,f,g,beta,a,dI,CD4(i-1),VIRUS(i-1)*10^3,sumLR),...
        [times(i-1) times(i)],[Inf(i-1)]);
    Inf(i)=C(end);
    % LR cells
    fprintf('%d\n',n);
    aa(1:n)=rL(1:n).*LR(1,1:n);
    aa((n+1):2*n)=(a+rL(1:n)*(CD4(i-1))/K).*LR(1,1:n);
    aa(2*n+1)=(f*beta*(VIRUS(i-1))*(10^3)*CD4(i-1));
    events=poissrnd(aa*instep);
    LR(2,1:n)=LR(1,1:n)+events(1:n);
    LR(2,1:n)=LR(2,1:n)-min(LR(2,1:n),events((n+1):2*n));
    if events(2*n+1)>0
        LR(2,(n+1):(n+events(2*n+1)))=ones(1,events(2*n+1));
        n=n+events(2*n+1);
    end
    if times(i)>=time_index
        %        fprintf('%2g \n',TIME(i));
        LR_save(save_index,1:size(LR,2))=LR(2,1:size(LR,2));
%         [~, col, v] = find(sparse(LR(2,:)));
%         if isempty(col)
%             fwrite(fileID4,[0 0]');
%             %fprintf(fileID4,'%d\t',[0 0]');
%             %fprintf(fileID4,'\n');
%         else
%             fwrite(fileID4,[col v]');
%             %fprintf(fileID4,'%d\t',[col v]');
%             %fprintf(fileID4,'\n');
%         end
%         fprintf(fileID2,'%6g ',Inf(i));
%         fprintf(fileID2,'\t');
%         fprintf(fileID3,'%2g ',times(i));
%         fprintf(fileID3,'\t');
        time_index=time_index+1;
        save_index=save_index+1;
        %         fprintf('%d \n',n);
    end
    %    T4(1)=T4(2);
    LR_save(i,1:size(LR,2))=LR(2,1:size(LR,2));
    LR(1,:)=LR(2,:);
    sumLR = sum(LR(1,:),2);
    %fprintf('%d %d \n',i,sumLR)
end


% LR_save = sparse(LR_save);
% [row, col, v] = find(LR_save);
% fileID4 = fopen(sprintf('blind_LR.txt'),'w');
% fprintf(fileID4,'%d \t',[row col v]');
% fclose(fileID4);
% 
% fclose(fileID2);
% fclose(fileID3);
% 
% 
% fileID = fopen(sprintf('ProliferationRate.txt'),'w');
% % fprintf(fileID, 'Results nn \n');
% fprintf(fileID,'%d ',rL(1:n));
% fprintf(fileID,'\t');
% fclose(fileID);
% 
% [~,idx] = ismember(times(round(times)==times),times);
% times=times(idx);
% Inf=Inf(idx);
% rL=rL(1:n);

save data.mat f r n sigma rL VIRUS CD4 times Inf LR_save

end
