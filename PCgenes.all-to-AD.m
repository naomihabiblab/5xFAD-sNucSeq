% PCgenes.all-to-AD.m

% read data
[id,idx]=textread('ADmouse.ASC.7m.Clusters.txt','%s%d','headerlines',1,'delimiter','\t');
[~,d1,d2]=textread('ADmouse.ASC.7m.Diffusion-Map-Cell-Embeddings.txt','%s%f%f','headerlines',1,'delimiter','\t'); DC=[d1,d2];
[PC] = dlmread('ADmouse.ASC.7m.PCA-Cell-Embeddings.txt','\t',1,1);
[genes]=textread('ADmouse.ASC.7m.PCA-Gene-Loadings.txt','%s%*[^\n]','headerlines',1,'delimiter','\t');
[PCA]  = dlmread('ADmouse.ASC.7m.PCA-Gene-Loadings.txt','\t',1,1);
[N,K] = size(PC);
subplot = @(m,n,p) subtightplot(m,n,p,0.05,[0.05 0.05],[0.01 0.01]);
WT = cellfun(@length,regexp(upper(id),'WT'))>0; AD = ~WT; WT=find(WT);AD=find(AD);

k = max(idx);
K = 25; PC = PC(:,1:K); PCA = PCA(:,1:K);

if ~exist('PCgenes.all-to-AD.mat','file'),
    % calculate pairwise distances
    pd=squareform(pdist(PC));
    % run Munkres' Hungarian algorithm on all WT and AD cells vs each AD cluster
    parfor c=1:6,
	ii=AD(idx(AD)==c); % all AD cells in cluster c 
	jj=[AD(idx(AD)~=c);WT]; % all other cells

	% debug on smaller groups
	if 0,
	    ii=ii(randperm(length(ii),min(length(ii),1000)));
	    jj=jj(randperm(length(jj),1000));
	end

	% build reduced distance matrix
	pd0 = pd(ii,jj); size(pd0),;
	% search for optimal assignment
	tic; [p,l]=munkres(pd0); toc;
	% organize and save results
	pp{c}=p;
	iii{c}=ii;
	jjj{c}=jj;
	% compare to random sets of distances
	% L = nan(1000,1); for i=1:1000, L(i)=sum(sqrt(sum((PC(ii,:)-PC(jj(randperm(length(jj),length(ii))),:)).^2,2))); end;
	% l = sum(sqrt(sum((PC(ii,:) - PC(jj(p),:)) .^ 2,2))); 
	% [(l-mean(L))/std(L)],;
    end

    % collect results from parallel runs
    from = nan(N,1);
    for c=1:6,
	ii = iii{c};
	jj = jjj{c};
	from(ii)=jj(pp{c});
	% [histc(idx(ii),1:6)';histc(idx(from(ii)),1:6)'; 100*histc(idx(from(ii)),1:6)'./histc(idx(jj),1:6)']
    end
    clear pp iii jjj pd pd0 p fid C P J
    save('PCgenes.all-to-AD.mat');
else
    load('PCgenes.all-to-AD.mat');
end

ii=AD; jj=from(ii);
DC0 = DC; DC0(ii,:) = DC(jj,:);
PC0 = PC; PC0(ii,:) = PC(jj,:);

% transitions
for i=1:6, Z(i,:)=histc(idx(jj(idx(ii)==i)),1:6)'; end;

% or use KNN to find a weighted-average embedding of WT and AD cells onto WT cells
if 0, %knn
    [knnidx,knnD] = knnsearch(PC(WT,:),PC,'K',15);
    R3=1./max(.1,knnD).^2; W3=bsxfun(@rdivide,R3,sum(R3,2));
    DC0=nan(N,2); for i=1:N, DC0(i,:) = W3(i,:) * DC(WT(knnidx(i,:)),:); end;
    PC0=nan(N,K); for i=1:N, PC0(i,:) = W3(i,:) * PC(WT(knnidx(i,:)),:); end;
end

figure(1);
subplot(3,2,1); scatter(DC(WT,1),DC(WT,2),5,idx(WT),'filled'); set(gca,'XTick',[],'YTick',[]); ax=axis; title('WT');
subplot(3,2,2); scatter(DC(AD,1),DC(AD,2),5,idx(AD),'filled'); set(gca,'XTick',[],'YTick',[]); axis(ax); title('AD');
for i=1:max(idx),
    subplot(3,3,3+i);
    I=AD(idx(AD)==i);
    quiver(DC(I,1),DC(I,2),DC(I,1)-DC0(I,1),DC(I,2)-DC0(I,2),1,'k','ShowArrowHead','on','MaxHeadSize',.2);
    set(gca,'XTick',[],'YTick',[]); axis(ax);
    hold on;
    scatter(DC(AD,1),DC(AD,2),3,idx(AD),'filled');
    for j=[I'], plot([DC0(j,1),DC(j,1)],[DC0(j,2),DC(j,2)],'k:','LineWidth',1); end;
    hold off;
end
set(gcf,'PaperSize',[11 12],'PaperPosition',[0 0 11 12]); print(gcf,'-dpdf','-r300','astroshit-hungarian.all-to-AD.pdf');

% analyze cell-to-cell transitions
for c = 1:6 
    ii=AD(idx(AD)==c); % ii=ii(randperm(length(ii),500));
    jj = from(ii); 
    jj_WT = intersect(WT,jj);
    jj_AD = intersect(AD,jj); 

    current(c,:)=histc(idx(ii),1:6)';
    original(c,:)=histc(idx(jj),1:6)';

    fid=fopen(sprintf('Cluster%d_AD_origins.all-to-AD.stats.txt',c),'w');
    fprintf(fid, 'Origins of cluster %d (AD):\n',c);
    fprintf(fid, 'all:        '); fprintf(fid, '\t%d', histc(idx(jj),1:6)'); fprintf(fid,'\n');
    fprintf(fid, 'WT:         '); fprintf(fid, '\t%d', histc(idx(jj_WT),1:6)'); fprintf(fid,'\n');
    fprintf(fid, '(%% of orig)'); fprintf(fid, '\t%.1f%%', 100*(histc(idx(jj_WT),1:6)./histc(idx(WT),1:6))); fprintf(fid,'\n');
    fprintf(fid, 'AD:         '); fprintf(fid, '\t%d', histc(idx(jj_AD),1:6)'); fprintf(fid,'\n');
    fprintf(fid, '(%% of orig)');
    fprintf(fid, '\t%.1f%%', 100*(histc(idx(jj_AD),1:6)./(histc(idx(AD),1:6)-histc(idx(ii),1:6)))); fprintf(fid,'\n');
    fclose(fid);

    fid=fopen(sprintf('Cluster%d_AD_origins.all-to-AD.genes.txt',c),'w');
    for i=1:6;
	[C,P]=corr(mean(PC(ii(idx(jj)==i),:) - PC(jj(idx(jj)==i),:),1)',PCA');
	[~,J]=sort(C,'descend'); J=J(find(FDR(P(J))<0.001));
	fprintf(fid, 'cl %d -> AD %d [%3d]:\t',i,c,sum(idx(jj)==i));
	for j=[J], fprintf(fid,'%s(%.1f,%.1g),', genes{j},C(j),P(j)); end;
	fprintf(fid, '\n');

	[C,P]=corr(mean(PC(ii(idx(jj_AD)==i),:) - PC(jj(idx(jj_AD)==i),:),1)',PCA');
	[~,J]=sort(C,'descend'); J=J(find(FDR(P(J))<0.001));
	fprintf(fid, 'AD %d -> AD %d [%3d]:\t',i,c,sum(idx(jj_AD)==i));
	for j=[J], fprintf(fid,'%s(%.1f,%.1g),', genes{j},C(j),P(j)); end;
	fprintf(fid, '\n');

	[C,P]=corr(mean(PC(ii(idx(jj_WT)==i),:) - PC(jj(idx(jj_WT)==i),:),1)',PCA');
	[~,J]=sort(C,'descend'); J=J(find(FDR(P(J))<0.001));
	fprintf(fid, 'WT %d -> AD %d [%3d]:\t',i,c,sum(idx(jj_WT)==i));
	for j=[J], fprintf(fid,'%s(%.1f,%.1g),', genes{j},C(j),P(j)); end;
	fprintf(fid, '\n');
    end
    fclose(fid);
end

clear ii jj pd pd0 p fid C P J

trans = original - current;
origins = bsxfun(@rdivide,100*original,diag(current));

Z1 = original;
Z2 = origins;
colormap(cool);

figure(1); clf;
X=repmat([1:k],k,1); Y=repmat([1:k]',1,k); I = find(X>0);
scatter(X(I),Y(I),1+3*Z1(I),Z2(I),'filled'); axis([.5 k+.5 .5 k+.5]);
set(gca,'YDIr','reverse','XTick',[1:k],'YTick',[1:k],'CLim',[0 50]); colorbar;
hold on; for i=0:k, plot([.5 k+.5],0.5+[i i],'k-',0.5+[i i],[.5 k+.5],'k-'); end; hold off;
set(gca,'FontSize',16); ylabel('AD cluster ID'); xlabel('Cluster of nearest cell');
set(gcf,'PaperSize',[10 9],'PaperPosition',[0 0 10 9]); print(gcf,'-dpdf','-r300','astroshit-by-size.all-to-AD.pdf');

% find transitions
fid  = fopen('Transitions.all-to-AD.txt','wt');
fid1 = fopen('Transitions.AD-to-AD.txt','wt');
fid2 = fopen('Transitions.WT-to-AD.txt','wt');
for c = 1:6 
    ii=AD(idx(AD)==c);
    jj = from(ii); 
    jj_WT = intersect(WT,jj);
    jj_AD = intersect(AD,jj); 

    fprintf(fid, 'all'); fprintf(fid, '\t%d', histc(idx(jj),1:6)'); fprintf(fid,'\n');
    fprintf(fid1, 'AD'); fprintf(fid1, '\t%d', histc(idx(jj_AD),1:6)'); fprintf(fid1,'\n');
    fprintf(fid2, 'WT'); fprintf(fid2, '\t%d', histc(idx(jj_WT),1:6)'); fprintf(fid2,'\n');
end
fclose(fid); fclose(fid1); fclose(fid2);

