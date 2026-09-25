################################################################################
## Zero-fill variant of the pool PTBP1 co-target/non-target violin: uncovered
## positions in each 200-nt window are set to signal=0 (all sim positions used),
## instead of the inner-join (covered-only) default. Same sim/params otherwise.
################################################################################
suppressPackageStartupMessages({
  library(data.table); library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(Biostrings); library(BSgenome.Hsapiens.UCSC.hg19); library(RBPEqBind) })
WIN=200; MINPTS=20; FLANK=50
revDir='/Users/soonyi/Repos/BITS_Specificity/REVISION_2026-09/'
simDir=path.expand('~/Repos/BITS_Specificity/Dataset/Analysis/RBPEqBind_Simulation')
bgBedDir=file.path(simDir,'RESULTS_ANALYSIS/BEDGRAPH_PROCESSED'); sitesFile=file.path(simDir,'DATA_PROCESSED/binding_sites.csv')
ptbp2_fa=file.path(simDir,'DATA_PROCESSED/FASTA/PTBP2.fa'); genome=BSgenome.Hsapiens.UCSC.hg19; sc=tempdir()
mr=read_csv(file.path(simDir,'DATA_PROCESSED/rnacompete_affinity_scores.csv'),show_col_types=FALSE);colnames(mr)[1]='Motif';tm=file.path(sc,'m.csv');write_csv(mr,tm);RM=loadModel(tm)
so=function(fa,r,cc){rm=setModel(RM,max_affinity=setNames(1/rep(50,length(r)),r),min_affinity=1e-5);simulateBindingF(fasta_file=fa,rbp_models=rm[r],protein_concs=cc,rna_conc=6.75,k=7)}
bg_all=lapply(c(U='bg_u2af2_500nM.bedgraph',H='bg_u2af2_500nM_hnrnpc_200nM.bedgraph',P='bg_u2af2_500nM_ptbp1_200nM.bedgraph'),
  function(f)as.data.table(read_tsv(file.path(bgBedDir,f),col_names=c('chr','start','end','val'),show_col_types=FALSE)))
bg_track=function(cond,chrom,gs,ge){d=bg_all[[cond]][chr==chrom&start>=gs-1&end<=ge+1];if(!nrow(d))return(data.table(pos=integer(),signal=numeric()))
  p=unlist(lapply(seq_len(nrow(d)),function(i)seq.int(d$start[i],d$end[i]-1L)));data.table(pos=p-gs+1L,signal=rep(d$val,d$end-d$start))[!duplicated(pos)]}
sites=read_csv(sitesFile,show_col_types=FALSE)%>%mutate(chr=ifelse(grepl('^chr',chromosome),chromosome,paste0('chr',tolower(sub('Chr','',chromosome)))),start=as.numeric(start),end=as.numeric(end))
loci=sites%>%group_by(transcript)%>%summarise(chr=chr[1],smin=min(start),smax=max(end),.groups='drop')

## per window, ZERO-FILL: all sim positions in window, signal=bg or 0
cor_zf=function(simd,bg,lo,hi){s=simd[pos>=lo&pos<hi]; if(nrow(s)<MINPTS)return(NA)
  m=merge(s,bg,by='pos',all.x=TRUE); m$signal[is.na(m$signal)]=0; if(sd(m$density)==0||sd(m$signal)==0)return(NA); cor(m$density,m$signal)}

run_tx=function(tx){L=loci[loci$transcript==tx,];ch=L$chr
  if(tx=='PTBP2'){gs=97269727;ge=97272451;fa=ptbp2_fa}else{gs=round(L$smin-FLANK);ge=round(L$smax+FLANK)
    dss=DNAStringSet(getSeq(genome,ch,gs,ge));fa=file.path(sc,paste0(tx,'.fa'));writeXStringSet(setNames(dss,sprintf('%s|%s:%d-%d|+|len=%d',tx,ch,gs,ge,width(dss))),fa)}
  tn=sub('^>','',readLines(fa,n=1))
  simd=lapply(c('U','H','P'),function(cond){cc=switch(cond,U=c(U2AF2=500),H=c(U2AF2=500,HNRNPC=200),P=c(U2AF2=500,PTBP1=200));r=switch(cond,U='U2AF2',H=c('U2AF2','HNRNPC'),P=c('U2AF2','PTBP1'))
    so(fa,r,cc)[transcript==tn,.(pos,density=U2AF2_density)]});names(simd)=c('U','H','P')
  bg=lapply(c('U','H','P'),function(cond)bg_track(cond,ch,gs,ge));names(bg)=c('U','H','P')
  if(any(sapply(bg,nrow)<MINPTS))return(NULL)
  st=sites[sites$transcript==tx,];st$ctr=round((st$start+st$end)/2)-gs+1;cov=range(simd$U$pos)
  st=st[st$ctr-WIN/2>=cov[1]&st$ctr+WIN/2<=cov[2],]
  bind_rows(lapply(seq_len(nrow(st)),function(i){lo=st$ctr[i]-WIN/2;hi=st$ctr[i]+WIN/2
    data.frame(transcript=tx,U=cor_zf(simd$U,bg$U,lo,hi),H=cor_zf(simd$H,bg$H,lo,hi),P=cor_zf(simd$P,bg$P,lo,hi))}))}

pool=bind_rows(lapply(loci$transcript,function(t)tryCatch(run_tx(t),error=function(e)NULL)))%>%filter(!is.na(U),!is.na(H),!is.na(P))
cls=read_csv(paste0(revDir,'3_fig4G_cotarget_classification.csv'),show_col_types=FALSE)
mk_labels=function(gcol,pcol){sub=cls[cls$transcript%in%unique(pool$transcript),];sub=sub[order(-sub[[pcol]],sub$transcript),]
  txs=tapply(sub$transcript,sub[[gcol]],paste,collapse=', ');setNames(paste0(names(txs),'\n(',txs,')'),names(txs))}
lab=mk_labels('PTBP1_group','PTBP1_peaks')
d=pool%>%left_join(cls%>%select(transcript,PTBP1_group),by='transcript')%>%mutate(grp=lab[PTBP1_group])
write_csv(d%>%mutate(across(where(is.numeric),~round(.x,4))),paste0(revDir,'3_fig4G_cotarget_zerofill.csv'),na='')

wp=function(a,b)suppressWarnings(wilcox.test(a,b,paired=TRUE)$p.value)
cat('=== ZERO-FILL, U2AF2 vs +PTBP1 paired per PTBP1 group ===\n')
print(d%>%group_by(PTBP1_group)%>%summarise(n=n(),medU=round(median(U),3),medP=round(median(P),3),dP=round(median(P-U),3),p=signif(wp(P,U),2),.groups='drop')%>%as.data.frame())
cat('(inner-join reference: co-target dP +0.031 p=2.7e-4 ; non-target dP -0.042 p=8.8e-4)\n')

stars=function(p)ifelse(p<1e-4,'****',ifelse(p<1e-3,'***',ifelse(p<1e-2,'**',ifelse(p<0.05,'*','ns'))))
cmp=list(c('U','H'),c('H','P'),c('U','P'));xmap=c(U=1,H=2,P=3);ypos=c(0.82,0.91,1.00)
br=bind_rows(lapply(unique(d$grp),function(g){dd=d[d$grp==g,];bind_rows(lapply(seq_along(cmp),function(i){a=cmp[[i]][1];b=cmp[[i]][2]
  data.frame(grp=g,xmin=xmap[a],xmax=xmap[b],y=ypos[i],label=sprintf('%s (p=%.1g)',stars(wp(dd[[a]],dd[[b]])),wp(dd[[a]],dd[[b]])))}))}))
cols=c('U2AF2 only'='#34495E','U2AF2 + HNRNPC'='#E91E63','U2AF2 + PTBP1'='#F39C12')
long=d%>%transmute(grp,`U2AF2 only`=U,`U2AF2 + HNRNPC`=H,`U2AF2 + PTBP1`=P)%>%pivot_longer(-grp,names_to='Condition',values_to='Pearson')%>%mutate(Condition=factor(Condition,levels=names(cols)))
zerofill_violin=ggplot(long,aes(Condition,Pearson))+geom_violin(aes(fill=Condition),alpha=0.5,colour=NA)+geom_boxplot(aes(fill=Condition),width=0.22,notch=TRUE,outlier.shape=NA)+
  geom_jitter(width=0.08,size=0.4,alpha=0.3)+facet_wrap(~grp)+scale_fill_manual(values=cols)+
  ggsignif::geom_signif(data=br,aes(xmin=xmin,xmax=xmax,annotations=label,y_position=y),manual=TRUE,inherit.aes=FALSE,textsize=3,tip_length=0.01)+
  scale_y_continuous(breaks=seq(0,0.75,0.25))+coord_cartesian(ylim=c(min(long$Pearson,na.rm=TRUE),1.1))+
  labs(title='PTBP1 co-target vs non-target -- ZERO-FILL (uncovered positions = 0), 200 nt',x=NULL,y='Per-window Pearson r')+
  theme_bw()+theme(legend.position='none',plot.title=element_text(hjust=0.5,size=11),axis.text.x=element_text(size=8,angle=20,hjust=1),axis.title=element_text(size=12,face='bold'),strip.text=element_text(size=9,face='bold'))
dir.create(file.path(revDir, 'figures'), showWarnings = FALSE)
ggsave(file.path(revDir, 'figures', '3_fig4G_cotarget_violin.pdf'), zerofill_violin, width = 8.5, height = 4.8)
cat('\nbuilt + saved zerofill_violin\n')
