#!/usr/bin/env bash
#SBATCH --job-name=eQTL_SuSiE_coloc
#SBATCH --output=/home/zw529/donglab/data/target_ALS/MR/eQTL_SuSiE_coloc.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/MR/eQTL_SuSiE_coloc.err
#SBATCH --time=08:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=4

set -euo pipefail
module --force purge
module load PLINK/1.9b_7.11-x86_64
module load R

command -v plink >/dev/null || { echo "ERROR: PLINK unavailable"; exit 1; }
command -v Rscript >/dev/null || { echo "ERROR: Rscript unavailable"; exit 1; }
Rscript -e 'p<-c("data.table","susieR","coloc","ggplot2"); x<-p[!sapply(p,requireNamespace,quietly=TRUE)]; if(length(x)) stop("Missing R packages: ",paste(x,collapse=", "))'

BASE="$HOME/donglab/data/target_ALS"; PLINK_DIR="$BASE/QTL/plink"; BFILE="$PLINK_DIR/joint_all_chrs_filtered_bed"; RAW="$PLINK_DIR/joint_all_chrs_matrixEQTL.raw"
GWAS="$HOME/donglab/data/GCST90027163/GWAS/harmonised/34873335-GCST90027163-MONDO_0004976.h.tsv.gz"; GWAS_ORIG="$HOME/donglab/data/GCST90027163/GWAS/GCST90027163_buildGRCh37.tsv.gz"
TISSUES=(Cervical_Spinal_Cord Lumbar_Spinal_Cord Motor_Cortex Frontal_Cortex Cerebellum)

CANDIDATE_MODE="global_fdr"   # global_fdr | tissue_fdr | nominal | all
H4_REFERENCE=0.80; MIN_SHARED_SNPS=50; MIN_SHARED_FRAC=0.50; SUSIE_L=10

export BASE BFILE RAW GWAS GWAS_ORIG CANDIDATE_MODE H4_REFERENCE MIN_SHARED_SNPS MIN_SHARED_FRAC SUSIE_L TISSUE_LIST="${TISSUES[*]}"

Rscript --vanilla - <<'RS'
suppressPackageStartupMessages({library(data.table); library(susieR); library(coloc); library(ggplot2)})

BASE<-Sys.getenv("BASE"); BFILE<-Sys.getenv("BFILE"); RAW<-Sys.getenv("RAW"); GWAS<-Sys.getenv("GWAS"); GWAS_ORIG<-Sys.getenv("GWAS_ORIG")
TISSUES<-strsplit(Sys.getenv("TISSUE_LIST")," ")[[1]]; MODE<-Sys.getenv("CANDIDATE_MODE"); H4REF<-as.numeric(Sys.getenv("H4_REFERENCE"))
MINSHARED<-as.integer(Sys.getenv("MIN_SHARED_SNPS")); MINFRAC<-as.numeric(Sys.getenv("MIN_SHARED_FRAC")); L<-as.integer(Sys.getenv("SUSIE_L"))
comp<-c(A="T",T="A",C="G",G="C"); nch<-function(x) toupper(sub("^chr","",as.character(x),ignore.case=TRUE))

# SNP-FIRST MR FILES
candidate_file<-function(t,o) switch(MODE,
  global_fdr=file.path(o,paste0(t,"_SNP_MR_global_FDR0.05.tsv")),
  tissue_fdr=file.path(o,paste0(t,"_SNP_MR_FDR0.05.tsv")),
  nominal=file.path(o,paste0(t,"_SNP_MR_nominal_P0.05.tsv")),
  all=file.path(o,paste0(t,"_SNP_MR_with_global_FDR.tsv.gz")),
  stop("Unknown CANDIDATE_MODE: ",MODE))

empty_cs<-function(file,status) fwrite(data.table(status=status,credible_set=NA_integer_,snps=NA_character_,n_snps=0L),file,sep="\t")
write_cs<-function(S,file,snps,label){
  z<-S$sets$cs
  if(is.null(z)||!length(z)){empty_cs(file,paste0("no_",label,"_credible_set")); return(0L)}
  fwrite(rbindlist(lapply(seq_along(z),function(i)data.table(status="OK",credible_set=i,snps=paste(snps[z[[i]]],collapse=";"),n_snps=length(z[[i]])))),file,sep="\t")
  length(z)
}

cat("Reading BIM + MatrixEQTL dosage orientation...\n")
bim<-fread(paste0(BFILE,".bim"),header=FALSE,col.names=c("chr","rawid","cm","pos","a1","a2")); bim[,key:=paste0(nch(chr),":",pos)]
hdr<-scan(RAW,what="",nlines=1,quiet=TRUE); hdr<-hdr[!hdr%in%c("FID","IID","PAT","MAT","SEX","PHENOTYPE")]
rawmap<-data.table(rawcol=hdr); rawmap[,`:=`(rawid=sub("_[^_]+$","",rawcol),effect_allele=toupper(sub("^.*_","",rawcol)))]
amap<-merge(rawmap,bim[,.(rawid,key,a1=toupper(a1),a2=toupper(a2))],by="rawid",all.x=TRUE)
amap[,other_allele:=fifelse(effect_allele==a1,a2,fifelse(effect_allele==a2,a1,NA_character_))]
amap<-unique(amap[!is.na(key),.(key,rawid,effect_allele,other_allele)],by="key"); rm(rawmap,hdr,bim); gc()

cat("Reading harmonized ALS GWAS...\n")
gwas<-fread(GWAS,select=c("hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","standard_error","p_value"))
setnames(gwas,c("hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","standard_error","p_value"),
              c("snpid","gwas_chr","gwas_pos","gwas_oa","gwas_ea","gwas_beta","gwas_se","gwas_p"))
gwas[,`:=`(snpid=as.character(snpid),gwas_chr=nch(gwas_chr),gwas_pos=as.integer(gwas_pos),gwas_ea=toupper(gwas_ea),gwas_oa=toupper(gwas_oa),
           gwas_beta=as.numeric(gwas_beta),gwas_se=as.numeric(gwas_se),gwas_p=as.numeric(gwas_p))]
gwas<-unique(gwas[!is.na(snpid)],by="snpid")

gn<-fread(GWAS_ORIG,select=c("rsid","N_effective")); GWAS_N<-round(median(as.numeric(gn$N_effective),na.rm=TRUE)); rm(gn); gc()
cat("GWAS effective N =",GWAS_N,"\n")
fam<-fread(paste0(BFILE,".fam"),header=FALSE,select=1:2,col.names=c("FID","IID")); master<-list()

for(tissue in TISSUES){
  mrdir<-file.path(BASE,tissue,"MR"); candfile<-candidate_file(tissue,mrdir); outroot<-file.path(mrdir,"SuSiE_coloc")
  dir.create(outroot,recursive=TRUE,showWarnings=FALSE)

  fullfile<-file.path(BASE,tissue,"eQTL","results",paste0(tissue,"_eQTL.full_annotated.txt"))
  locfile<-file.path(BASE,tissue,"eQTL","snp_location.txt")
  covfile<-file.path(BASE,tissue,"eQTL",paste0("covariates_",tissue,"_encoded.txt"))

  if(!file.exists(candfile)||file.info(candfile)$size==0){cat("\n",tissue,": no candidate file; skipping\n",sep=""); next}
  for(f in c(fullfile,locfile,covfile)) if(!file.exists(f)||file.info(f)$size==0) stop("Missing/empty: ",f)

  # SNP-FIRST MR rows -> genes are nominated only here, downstream
  cand<-fread(candfile)
  if(!nrow(cand)){cat("\n",tissue,": zero candidates\n",sep=""); next}
  if(!all(c("snpid","geneid","p_MR")%in%names(cand))) stop(tissue,": SNP-first candidate file missing snpid/geneid/p_MR")

  cand[,p_MR:=as.numeric(p_MR)]
  if("FDR_MR_global"%in%names(cand)) cand[,FDR_MR_global:=as.numeric(FDR_MR_global)]
  if("FDR_MR_tissue"%in%names(cand)) cand[,FDR_MR_tissue:=as.numeric(FDR_MR_tissue)]

  gene_meta<-cand[,{
    z<-.SD[order(p_MR)]
    .(n_MR_SNPs=uniqueN(snpid),MR_SNPs=paste(unique(snpid),collapse=";"),
      best_MR_SNP=as.character(z$snpid[1]),best_MR_p=as.numeric(z$p_MR[1]),
      best_MR_global_FDR=if("FDR_MR_global"%in%names(z)) min(z$FDR_MR_global,na.rm=TRUE) else NA_real_,
      best_MR_tissue_FDR=if("FDR_MR_tissue"%in%names(z)) min(z$FDR_MR_tissue,na.rm=TRUE) else NA_real_)
  },by=geneid]
  gene_meta[!is.finite(best_MR_global_FDR),best_MR_global_FDR:=NA_real_]
  gene_meta[!is.finite(best_MR_tissue_FDR),best_MR_tissue_FDR:=NA_real_]

  genes<-gene_meta$geneid
  cat("\n",strrep("=",70),"\n",tissue,": ",nrow(cand)," SNP→gene MR candidates nominating ",length(genes)," genes\n",strrep("=",70),"\n",sep="")

  full<-fread(fullfile)[geneid%in%genes]
  loc<-fread(locfile); loc[,key:=paste0(nch(chr),":",pos)]; loc<-unique(loc,by="snpid")

  covhdr<-names(fread(covfile,nrows=0)); subjects<-setdiff(covhdr,covhdr[1]); keep<-fam[IID%in%subjects]
  if(!nrow(keep)) stop(tissue,": no covariate subjects matched PLINK FAM")
  N_EQTL<-nrow(keep); keepfile<-file.path(outroot,"plink_keep.txt"); fwrite(keep,keepfile,sep="\t",col.names=FALSE)
  cat("eQTL N =",N_EQTL,"\n")

  for(gene in genes){
    cat("  ",gene," ... ",sep="")
    gout<-file.path(outroot,gsub("[^A-Za-z0-9_.-]","_",gene)); dir.create(gout,showWarnings=FALSE)

    q<-copy(full[geneid==gene]); gm<-gene_meta[geneid==gene][1]
    q[,`:=`(snpid=as.character(snpid),qtl_chr=nch(chr),qtl_pos=as.integer(pos),beta_qtl=as.numeric(beta),
             se_qtl=abs(as.numeric(beta)/as.numeric(`t-stat`)),p_qtl=as.numeric(`p-value`))]
    q<-q[!is.na(snpid)&is.finite(beta_qtl)&is.finite(se_qtl)&se_qtl>0]
    q<-merge(q[,.(geneid,snpid,qtl_chr,qtl_pos,beta_qtl,se_qtl,p_qtl)],loc[,.(snpid,key)],by="snpid",all.x=TRUE)
    q<-merge(q,amap,by="key",all.x=TRUE); q<-unique(q,by="snpid"); n_locus<-nrow(q)

    if(n_locus<MINSHARED){
      cat("too few locus SNPs (",n_locus,")\n",sep="")
      master[[length(master)+1]]<-data.table(tissue,geneid=gene,status="too_few_locus_snps",n_locus=n_locus,n_shared=0,shared_fraction=0,
                                             n_MR_SNPs=gm$n_MR_SNPs,best_MR_SNP=gm$best_MR_SNP,best_MR_p=gm$best_MR_p)
      next
    }

    x<-merge(q,gwas,by="snpid",all=FALSE)
    x<-x[qtl_chr==gwas_chr&qtl_pos==gwas_pos&!is.na(effect_allele)&!is.na(other_allele)]
    x[,pal:=paste0(effect_allele,other_allele)%in%c("AT","TA","CG","GC")]
    x[,relation:=fcase(effect_allele==gwas_ea&other_allele==gwas_oa,"exact",
                       effect_allele==gwas_oa&other_allele==gwas_ea,"swap",
                       effect_allele==comp[gwas_ea]&other_allele==comp[gwas_oa],"strand",
                       effect_allele==comp[gwas_oa]&other_allele==comp[gwas_ea],"strand_swap",default="bad")]
    x<-unique(x[!pal&relation!="bad"&is.finite(gwas_beta)&is.finite(gwas_se)&gwas_se>0],by="snpid")
    x[,gwas_beta_h:=fifelse(relation%in%c("swap","strand_swap"),-gwas_beta,gwas_beta)]
    shared_frac<-nrow(x)/n_locus

    if(nrow(x)<MINSHARED||shared_frac<MINFRAC){
      status<-if(nrow(x)<MINSHARED) "too_few_shared_snps" else "low_shared_coverage"
      cat(status,", shared=",nrow(x),"/",n_locus," (",round(100*shared_frac,1),"%)\n",sep="")
      master[[length(master)+1]]<-data.table(tissue,geneid=gene,status,n_locus=n_locus,n_shared=nrow(x),shared_fraction=shared_frac,
                                             n_MR_SNPs=gm$n_MR_SNPs,best_MR_SNP=gm$best_MR_SNP,best_MR_p=gm$best_MR_p)
      next
    }

    idsfile<-file.path(gout,"raw_ids.txt"); fwrite(x[,.(rawid)],idsfile,col.names=FALSE); pref<-file.path(gout,"ldgeno")
    system2("plink",c("--bfile",BFILE,"--keep",keepfile,"--extract",idsfile,"--recode","A-transpose","--threads",Sys.getenv("SLURM_CPUS_PER_TASK","1"),"--out",pref),
            stdout=TRUE,stderr=TRUE)
    trawfile<-paste0(pref,".traw")

    if(!file.exists(trawfile)){
      cat("PLINK failed\n")
      master[[length(master)+1]]<-data.table(tissue,geneid=gene,status="PLINK_failed",n_locus=n_locus,n_shared=nrow(x),shared_fraction=shared_frac,
                                             n_MR_SNPs=gm$n_MR_SNPs,best_MR_SNP=gm$best_MR_SNP,best_MR_p=gm$best_MR_p)
      next
    }

    tr<-fread(trawfile); geno_cols<-setdiff(names(tr),c("CHR","SNP","(C)M","POS","COUNTED","ALT"))
    if(nrow(tr)<MINSHARED||!length(geno_cols)){
      cat("too few PLINK SNPs\n")
      master[[length(master)+1]]<-data.table(tissue,geneid=gene,status="too_few_PLINK_snps",n_locus=n_locus,n_shared=nrow(x),shared_fraction=shared_frac,
                                             n_MR_SNPs=gm$n_MR_SNPs,best_MR_SNP=gm$best_MR_SNP,best_MR_p=gm$best_MR_p)
      next
    }

    G<-as.matrix(tr[,..geno_cols]); storage.mode(G)<-"numeric"
    orient<-x[match(tr$SNP,rawid),effect_allele]; flip<-toupper(tr$COUNTED)!=toupper(orient); G[flip,]<-2-G[flip,]
    G<-t(apply(G,1,function(v){v[is.na(v)]<-mean(v,na.rm=TRUE); v}))
    sdv<-apply(G,1,sd); ok<-is.finite(sdv)&sdv>0; tr<-tr[ok]; G<-G[ok,,drop=FALSE]; x<-x[match(tr$SNP,rawid)]

    if(nrow(x)<MINSHARED){
      cat("too few variable SNPs\n")
      master[[length(master)+1]]<-data.table(tissue,geneid=gene,status="too_few_variable_snps",n_locus=n_locus,n_shared=nrow(x),shared_fraction=nrow(x)/n_locus,
                                             n_MR_SNPs=gm$n_MR_SNPs,best_MR_SNP=gm$best_MR_SNP,best_MR_p=gm$best_MR_p)
      next
    }

    R<-cor(t(G))
    if(!is.matrix(R)||nrow(R)!=nrow(x)) stop(gene,": malformed LD matrix")
    R[R>1]<-1; R[R< -1]<- -1; diag(R)<-1; colnames(R)<-rownames(R)<-x$snpid

    d1<-list(beta=x$beta_qtl,varbeta=x$se_qtl^2,snp=x$snpid,position=x$qtl_pos,type="quant",sdY=1,N=N_EQTL,LD=R)
    d2<-list(beta=x$gwas_beta_h,varbeta=x$gwas_se^2,snp=x$snpid,position=x$qtl_pos,type="cc",N=GWAS_N,LD=R)

    invisible(coloc::check_dataset(d1,req="LD")); invisible(coloc::check_dataset(d2,req="LD"))
    z1<-d1$beta/sqrt(d1$varbeta); z2<-d2$beta/sqrt(d2$varbeta)
    s_eqtl<-tryCatch(susieR::estimate_s_rss(z1,R,n=N_EQTL),error=function(e) NA_real_)
    s_gwas<-tryCatch(susieR::estimate_s_rss(z2,R,n=GWAS_N),error=function(e) NA_real_)

    run_coloc_results<-function(S1,S2,x,gout){
      pip1<-if(is.null(S1$pip)) rep(NA_real_,nrow(x)) else S1$pip
      pip2<-if(is.null(S2$pip)) rep(NA_real_,nrow(x)) else S2$pip
      fwrite(data.table(snp=x$snpid,position=x$qtl_pos,eqtl_pip=pip1,gwas_pip=pip2),file.path(gout,"PIP.tsv"),sep="\t")

      ncs1<-write_cs(S1,file.path(gout,"eQTL_credible_sets.tsv"),x$snpid,"eQTL")
      ncs2<-write_cs(S2,file.path(gout,"GWAS_credible_sets.tsv"),x$snpid,"GWAS")
      top1<-if(any(is.finite(pip1))) which.max(pip1) else NA_integer_
      top2<-if(any(is.finite(pip2))) which.max(pip2) else NA_integer_

      ret<-function(status,h4=NA_real_) list(
        status=status,h4=h4,ncs1=ncs1,ncs2=ncs2,
        eqtl_top=if(!is.na(top1))x$snpid[top1] else NA_character_,
        eqtl_pip=if(!is.na(top1))pip1[top1] else NA_real_,
        gwas_top=if(!is.na(top2))x$snpid[top2] else NA_character_,
        gwas_pip=if(!is.na(top2))pip2[top2] else NA_real_)

      if(ncs1==0||ncs2==0){
        status<-if(ncs1==0&&ncs2==0) "no_eQTL_or_GWAS_credible_set" else if(ncs1==0) "no_eQTL_credible_set" else "no_GWAS_credible_set"
        fwrite(data.table(status=status),file.path(gout,"coloc_susie.tsv"),sep="\t"); return(ret(status))
      }

      C<-coloc::coloc.susie(S1,S2); csum<-if(is.null(C$summary)) data.table() else as.data.table(C$summary)
      if(!nrow(csum)){fwrite(data.table(status="no_coloc_signal_pairs"),file.path(gout,"coloc_susie.tsv"),sep="\t"); return(ret("no_coloc_signal_pairs"))}

      fwrite(csum,file.path(gout,"coloc_susie.tsv"),sep="\t")
      h4<-if("PP.H4.abf"%in%names(csum)) suppressWarnings(max(csum$PP.H4.abf,na.rm=TRUE)) else NA_real_
      if(!is.finite(h4)) h4<-NA_real_
      ret("OK",h4)
    }

    ans<-tryCatch({
      S1<-coloc::runsusie(d1,maxit=200,L=L); S2<-coloc::runsusie(d2,maxit=200,L=L)
      run_coloc_results(S1,S2,x,gout)
    },error=function(e) list(status=paste0("ERROR: ",conditionMessage(e)),h4=NA_real_,ncs1=NA_integer_,ncs2=NA_integer_,
                              eqtl_top=NA_character_,eqtl_pip=NA_real_,gwas_top=NA_character_,gwas_pip=NA_real_))

    fwrite(x[,.(snpid,qtl_chr,qtl_pos,effect_allele,other_allele,beta_qtl,se_qtl,p_qtl,gwas_beta_h,gwas_se,gwas_p,relation)],
           file.path(gout,"harmonized_locus.tsv"),sep="\t")

    rec<-data.table(
      tissue,geneid=gene,status=ans$status,n_locus=n_locus,n_shared=nrow(x),shared_fraction=nrow(x)/n_locus,
      N_eQTL=N_EQTL,N_GWAS=GWAS_N,n_MR_SNPs=gm$n_MR_SNPs,MR_SNPs=gm$MR_SNPs,best_MR_SNP=gm$best_MR_SNP,
      best_MR_p=gm$best_MR_p,best_MR_tissue_FDR=gm$best_MR_tissue_FDR,best_MR_global_FDR=gm$best_MR_global_FDR,
      eqtl_LD_mismatch_s=s_eqtl,gwas_LD_mismatch_s=s_gwas,n_eQTL_CS=ans$ncs1,n_GWAS_CS=ans$ncs2,
      top_eQTL_SNP=ans$eqtl_top,top_eQTL_PIP=ans$eqtl_pip,top_GWAS_SNP=ans$gwas_top,top_GWAS_PIP=ans$gwas_pip,
      max_PP_H4=ans$h4,H4_reference_pass=!is.na(ans$h4)&ans$h4>=H4REF)

    master[[length(master)+1]]<-rec
    cat(ans$status,", shared=",nrow(x),"/",n_locus," (",round(100*nrow(x)/n_locus,1),
        "%), MR_SNPs=",gm$n_MR_SNPs,", eQTL_CS=",ans$ncs1,", GWAS_CS=",ans$ncs2,", H4=",signif(ans$h4,3),"\n",sep="")
    rm(G,R,x,q); gc()
  }
  rm(full); gc()
}

res<-if(length(master)) rbindlist(master,fill=TRUE) else data.table()
outfile<-file.path(BASE,"MR","SuSiE_coloc_summary.tsv"); fwrite(res,outfile,sep="\t")

# ==============================================================
# PP.H4 DISTRIBUTION — NO FILTERING BY THE 0.80 REFERENCE VALUE
# ==============================================================
h4<-res[is.finite(max_PP_H4),.(tissue,geneid,max_PP_H4,n_MR_SNPs,best_MR_SNP,best_MR_p,best_MR_global_FDR)]
h4file<-file.path(BASE,"MR","PPH4_values.tsv"); fwrite(h4,h4file,sep="\t")

if(nrow(h4)){
  p<-ggplot(h4,aes(x=max_PP_H4))+
    geom_histogram(binwidth=0.05,boundary=0,closed="left")+
    geom_vline(xintercept=H4REF,linetype="dashed",linewidth=0.8)+
    scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.1))+
    labs(title=paste0("SuSiE colocalization PP.H4 distribution (n=",nrow(h4),")"),
         subtitle=paste0("Dashed line = ",H4REF," reference threshold; histogram includes all finite PP.H4 values"),
         x="Maximum PP.H4 per gene–tissue locus",y="Number of loci")+
    theme_bw(base_size=12)

  plotfile<-file.path(BASE,"MR","PPH4_histogram.png")
  ggsave(plotfile,plot=p,width=8,height=5,dpi=300)

  cat("\nPP.H4 distribution:\n")
  cat("  n =",nrow(h4),"\n")
  cat("  median =",median(h4$max_PP_H4),"\n")
  cat("  IQR =",quantile(h4$max_PP_H4,.25),"to",quantile(h4$max_PP_H4,.75),"\n")
  cat("  >=0.50 =",sum(h4$max_PP_H4>=.50),"\n")
  cat("  >=0.80 =",sum(h4$max_PP_H4>=.80),"\n")
  cat("  Values:",h4file,"\n")
  cat("  Histogram:",plotfile,"\n")
} else {
  cat("\nNo finite PP.H4 values; histogram not generated.\n")
}

cat("\n",strrep("=",70),"\nSuSiE / COLOC SUMMARY\n",strrep("=",70),"\n",sep="")
print(res); cat("\nOutput: ",outfile,"\n",sep="")
RS

echo
echo "SuSiE + colocalization complete."
column -t -s $'\t' "$BASE/MR/SuSiE_coloc_summary.tsv"
