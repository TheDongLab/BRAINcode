#!/usr/bin/env bash
#SBATCH --job-name=eQTL_SuSiE_coloc
#SBATCH --output=/home/zw529/donglab/data/target_ALS/MR/eQTL_SuSiE_coloc.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/MR/eQTL_SuSiE_coloc.err
#SBATCH --time=8:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4

set -euo pipefail
module --force purge
module load PLINK/1.9b_7.11-x86_64
module load R

command -v plink >/dev/null || { echo "ERROR: PLINK unavailable"; exit 1; }
command -v Rscript >/dev/null || { echo "ERROR: Rscript unavailable"; exit 1; }
Rscript -e 'pkgs<-c("data.table","susieR","coloc"); x<-pkgs[!sapply(pkgs,requireNamespace,quietly=TRUE)]; if(length(x)) stop("Missing R packages: ",paste(x,collapse=", "))'

BASE="$HOME/donglab/data/target_ALS"; PLINK_DIR="$BASE/QTL/plink"; BFILE="$PLINK_DIR/joint_all_chrs_filtered_bed"; RAW="$PLINK_DIR/joint_all_chrs_matrixEQTL.raw"
GWAS="$HOME/donglab/data/GCST90027163/GWAS/harmonised/34873335-GCST90027163-MONDO_0004976.h.tsv.gz"; GWAS_ORIG="$HOME/donglab/data/GCST90027163/GWAS/GCST90027163_buildGRCh37.tsv.gz"
TISSUES=(Cervical_Spinal_Cord Lumbar_Spinal_Cord Motor_Cortex Frontal_Cortex Cerebellum)
CANDIDATE_MODE="global_fdr"; COLOC_H4=0.80; MIN_SHARED_SNPS=100; SUSIE_L=10

export BASE BFILE RAW GWAS GWAS_ORIG CANDIDATE_MODE COLOC_H4 MIN_SHARED_SNPS SUSIE_L TISSUE_LIST="${TISSUES[*]}"

Rscript --vanilla - <<'RS'
suppressPackageStartupMessages({library(data.table); library(susieR); library(coloc)})

BASE<-Sys.getenv("BASE"); BFILE<-Sys.getenv("BFILE"); RAW<-Sys.getenv("RAW"); GWAS<-Sys.getenv("GWAS"); GWAS_ORIG<-Sys.getenv("GWAS_ORIG")
TISSUES<-strsplit(Sys.getenv("TISSUE_LIST")," ")[[1]]; MODE<-Sys.getenv("CANDIDATE_MODE"); H4CUT<-as.numeric(Sys.getenv("COLOC_H4")); MINSHARED<-as.integer(Sys.getenv("MIN_SHARED_SNPS")); L<-as.integer(Sys.getenv("SUSIE_L"))
comp<-c(A="T",T="A",C="G",G="C"); nch<-function(x) toupper(sub("^chr","",as.character(x),ignore.case=TRUE))

candidate_file<-function(t,o) switch(MODE,
  global_fdr=file.path(o,paste0(t,"_IndependentSNP_MR_global_FDR0.05.tsv")),
  tissue_fdr=file.path(o,paste0(t,"_IndependentSNP_MR_FDR0.05.tsv")),
  nominal=file.path(o,paste0(t,"_IndependentSNP_MR_nominal_P0.05.tsv")),
  all=file.path(o,paste0(t,"_IndependentSNP_MR_with_global_FDR.tsv.gz")),
  stop("Unknown CANDIDATE_MODE: ",MODE))

cat("Reading BIM + original MatrixEQTL dosage orientation...\n")
bim<-fread(paste0(BFILE,".bim"),header=FALSE,col.names=c("chr","rawid","cm","pos","a1","a2")); bim[,key:=paste0(nch(chr),":",pos)]
hdr<-scan(RAW,what="",nlines=1,quiet=TRUE); hdr<-hdr[!hdr%in%c("FID","IID","PAT","MAT","SEX","PHENOTYPE")]
rawmap<-data.table(rawcol=hdr); rawmap[,`:=`(rawid=sub("_[^_]+$","",rawcol),effect_allele=toupper(sub("^.*_","",rawcol)))]
amap<-merge(rawmap,bim[,.(rawid,key,a1=toupper(a1),a2=toupper(a2))],by="rawid",all.x=TRUE)
amap[,other_allele:=fifelse(effect_allele==a1,a2,fifelse(effect_allele==a2,a1,NA_character_))]
amap<-unique(amap[!is.na(key),.(key,rawid,effect_allele,other_allele)],by="key"); rm(rawmap,hdr,bim); gc()

cat("Reading full harmonized ALS GWAS...\n")
gwas<-fread(GWAS,select=c("hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","standard_error","p_value"))
setnames(gwas,c("hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","standard_error","p_value"),
              c("snpid","gwas_chr","gwas_pos","gwas_oa","gwas_ea","gwas_beta","gwas_se","gwas_p"))
gwas[,`:=`(snpid=as.character(snpid),gwas_chr=nch(gwas_chr),gwas_ea=toupper(gwas_ea),gwas_oa=toupper(gwas_oa))]
gwas<-unique(gwas[!is.na(snpid)],by="snpid")

cat("Determining external GWAS effective N...\n")
gn<-fread(GWAS_ORIG,select=c("rsid","N_effective")); GWAS_N<-round(median(as.numeric(gn$N_effective),na.rm=TRUE)); rm(gn); gc()
cat("GWAS effective N =",GWAS_N,"\n")

fam<-fread(paste0(BFILE,".fam"),header=FALSE,select=1:2,col.names=c("FID","IID"))
master<-list()

for(tissue in TISSUES){
  mrdir<-file.path(BASE,tissue,"MR"); candfile<-candidate_file(tissue,mrdir); outroot<-file.path(mrdir,"SuSiE_coloc"); dir.create(outroot,recursive=TRUE,showWarnings=FALSE)
  fullfile<-file.path(BASE,tissue,"eQTL","results",paste0(tissue,"_eQTL.full_annotated.txt"))
  locfile<-file.path(BASE,tissue,"eQTL","snp_location.txt"); covfile<-file.path(BASE,tissue,"eQTL",paste0("covariates_",tissue,"_encoded.txt"))
  if(!file.exists(candfile)||file.info(candfile)$size==0){cat("\n",tissue,": no candidate file; skipping\n",sep=""); next}
  for(f in c(fullfile,locfile,covfile)) if(!file.exists(f)||file.info(f)$size==0) stop("Missing/empty: ",f)

  cand<-fread(candfile); if(nrow(cand)==0){cat("\n",tissue,": zero candidates\n",sep=""); next}
  genes<-unique(as.character(cand$geneid)); cat("\n",strrep("=",70),"\n",tissue,": ",length(genes)," candidate genes\n",strrep("=",70),"\n",sep="")

  full<-fread(fullfile); full<-full[geneid%in%genes]
  loc<-fread(locfile); loc[,key:=paste0(nch(chr),":",pos)]; loc<-unique(loc,by="snpid")
  nms<-names(fread(covfile,nrows=0)); subjects<-setdiff(nms,nms[1]); N_EQTL<-length(subjects)
  keep<-fam[IID%in%subjects]; if(nrow(keep)==0) stop(tissue,": no covariate subjects matched PLINK FAM")
  keepfile<-file.path(outroot,"plink_keep.txt"); fwrite(keep,keepfile,sep="\t",col.names=FALSE)

  for(gene in genes){
    cat("  ",gene," ... ",sep=""); gout<-file.path(outroot,gsub("[^A-Za-z0-9_.-]","_",gene)); dir.create(gout,showWarnings=FALSE)
    q<-copy(full[geneid==gene]); mr<-cand[geneid==gene][1]
    if(nrow(q)<MINSHARED){cat("too few eQTL SNPs\n"); master[[length(master)+1]]<-data.table(tissue,geneid=gene,status="too_few_eqtl_snps"); next}

    q[,`:=`(snpid=as.character(snpid),qtl_chr=nch(chr),qtl_pos=as.integer(pos),beta_qtl=as.numeric(beta),se_qtl=abs(as.numeric(beta)/as.numeric(`t-stat`)))]
    q<-merge(q[,.(geneid,snpid,qtl_chr,qtl_pos,beta_qtl,se_qtl,p_qtl=as.numeric(`p-value`))],loc[,.(snpid,key)],by="snpid",all.x=TRUE)
    q<-merge(q,amap,by="key",all.x=TRUE); x<-merge(q,gwas,by="snpid",all=FALSE)
    x<-x[qtl_chr==gwas_chr & qtl_pos==as.integer(gwas_pos) & !is.na(effect_allele) & !is.na(other_allele)]

    x[,pal:=paste0(effect_allele,other_allele)%in%c("AT","TA","CG","GC")]
    x[,relation:=fcase(effect_allele==gwas_ea & other_allele==gwas_oa,"exact",
                       effect_allele==gwas_oa & other_allele==gwas_ea,"swap",
                       effect_allele==comp[gwas_ea] & other_allele==comp[gwas_oa],"strand",
                       effect_allele==comp[gwas_oa] & other_allele==comp[gwas_ea],"strand_swap",
                       default="bad")]
    x<-x[!pal & relation!="bad" & is.finite(beta_qtl) & se_qtl>0 & is.finite(gwas_beta) & gwas_se>0]
    x[,gwas_beta_h:=fifelse(relation%in%c("swap","strand_swap"),-as.numeric(gwas_beta),as.numeric(gwas_beta))]
    x<-unique(x,by="snpid")
    if(nrow(x)<MINSHARED){cat("too few shared SNPs\n"); master[[length(master)+1]]<-data.table(tissue,geneid=gene,status="too_few_shared_snps",n_shared=nrow(x)); next}

    idsfile<-file.path(gout,"raw_ids.txt"); fwrite(x[,.(rawid)],idsfile,col.names=FALSE)
    pref<-file.path(gout,"ldgeno")
    system2("plink",c("--bfile",BFILE,"--keep",keepfile,"--extract",idsfile,"--recode","A-transpose","--out",pref),stdout=TRUE,stderr=TRUE)
    trawfile<-paste0(pref,".traw"); if(!file.exists(trawfile)){cat("PLINK failed\n"); next}
    tr<-fread(trawfile); if(nrow(tr)<MINSHARED){cat("too few PLINK SNPs\n"); next}

    geno_cols<-setdiff(names(tr),c("CHR","SNP","(C)M","POS","COUNTED","ALT")); G<-as.matrix(tr[,..geno_cols]); storage.mode(G)<-"numeric"
    orient<-x[match(tr$SNP,rawid),effect_allele]; flip<-toupper(tr$COUNTED)!=toupper(orient)
    G[flip,]<-2-G[flip,]; G<-t(apply(G,1,function(v){v[is.na(v)]<-mean(v,na.rm=TRUE); v}))
    sdv<-apply(G,1,sd); ok<-is.finite(sdv)&sdv>0; tr<-tr[ok]; G<-G[ok,,drop=FALSE]; x<-x[match(tr$SNP,rawid)]
    R<-cor(t(G)); diag(R)<-1; R[pmax(-1,pmin(1,R))!=R]<-pmax(-1,pmin(1,R))[pmax(-1,pmin(1,R))!=R]
    colnames(R)<-rownames(R)<-x$snpid

    d1<-list(beta=x$beta_qtl,varbeta=x$se_qtl^2,snp=x$snpid,position=x$qtl_pos,type="quant",sdY=1,N=N_EQTL,LD=R)
    d2<-list(beta=x$gwas_beta_h,varbeta=as.numeric(x$gwas_se)^2,snp=x$snpid,position=x$qtl_pos,type="cc",N=GWAS_N,LD=R)
    s_eqtl<-tryCatch(susieR::estimate_s_rss(d1$beta/sqrt(d1$varbeta),R),error=function(e) NA_real_)
    s_gwas<-tryCatch(susieR::estimate_s_rss(d2$beta/sqrt(d2$varbeta),R),error=function(e) NA_real_)

    ans<-tryCatch({
      S1<-coloc::runsusie(d1,maxit=200,L=L); S2<-coloc::runsusie(d2,maxit=200,L=L); C<-coloc::coloc.susie(S1,S2)
      pip<-data.table(snp=x$snpid,position=x$qtl_pos,eqtl_pip=S1$pip,gwas_pip=S2$pip); fwrite(pip,file.path(gout,"PIP.tsv"),sep="\t")
      cswrite<-function(S,file){z<-S$sets$cs; if(is.null(z)||!length(z)) return(fwrite(data.table(),file,sep="\t"))
        fwrite(rbindlist(lapply(seq_along(z),function(i)data.table(credible_set=i,snps=paste(x$snpid[z[[i]]],collapse=";"),n_snps=length(z[[i]])))),file,sep="\t")}
      cswrite(S1,file.path(gout,"eQTL_credible_sets.tsv")); cswrite(S2,file.path(gout,"GWAS_credible_sets.tsv"))
      fwrite(as.data.table(C$summary),file.path(gout,"coloc_susie.tsv"),sep="\t")
      h4<-if(nrow(C$summary)) max(C$summary$PP.H4.abf,na.rm=TRUE) else NA_real_
      best<-if(nrow(C$summary)) C$summary[which.max(C$summary$PP.H4.abf),] else NULL
      list(status="OK",h4=h4,best=best,eqtl_top=x$snpid[which.max(S1$pip)],eqtl_pip=max(S1$pip),gwas_top=x$snpid[which.max(S2$pip)],gwas_pip=max(S2$pip),
           ncs1=if(is.null(S1$sets$cs))0 else length(S1$sets$cs),ncs2=if(is.null(S2$sets$cs))0 else length(S2$sets$cs))
    },error=function(e) list(status=paste0("ERROR: ",conditionMessage(e)),h4=NA_real_,best=NULL,eqtl_top=NA,eqtl_pip=NA,gwas_top=NA,gwas_pip=NA,ncs1=NA,ncs2=NA))

    rec<-data.table(tissue,geneid=gene,status=ans$status,n_shared=nrow(x),N_eQTL=N_EQTL,N_GWAS=GWAS_N,
                    MR_method=if("method"%in%names(mr))mr$method else NA,MR_p=if("p_MR"%in%names(mr))mr$p_MR else NA,
                    MR_global_FDR=if("FDR_MR_global"%in%names(mr))mr$FDR_MR_global else NA,
                    eqtl_LD_mismatch_s=s_eqtl,gwas_LD_mismatch_s=s_gwas,n_eQTL_CS=ans$ncs1,n_GWAS_CS=ans$ncs2,
                    top_eQTL_SNP=ans$eqtl_top,top_eQTL_PIP=ans$eqtl_pip,top_GWAS_SNP=ans$gwas_top,top_GWAS_PIP=ans$gwas_pip,
                    max_PP_H4=ans$h4,coloc_pass=!is.na(ans$h4)&ans$h4>=H4CUT)
    fwrite(x[,.(snpid,qtl_chr,qtl_pos,effect_allele,other_allele,beta_qtl,se_qtl,p_qtl,gwas_beta_h,gwas_se,gwas_p,relation)],file.path(gout,"harmonized_locus.tsv"),sep="\t")
    master[[length(master)+1]]<-rec; cat(ans$status,", shared=",nrow(x),", max H4=",signif(ans$h4,3),"\n",sep="")
    rm(G,R,x,q); gc()
  }
  rm(full); gc()
}

res<-rbindlist(master,fill=TRUE); fwrite(res,file.path(BASE,"MR","SuSiE_coloc_summary.tsv"),sep="\t")
cat("\n",strrep("=",70),"\nSuSiE/COLOC SUMMARY\n",strrep("=",70),"\n",sep="")
print(res); cat("\nOutput: ",file.path(BASE,"MR","SuSiE_coloc_summary.tsv"),"\n",sep="")
RS

echo
echo "SuSiE + colocalization complete."
column -t -s $'\t' "$BASE/MR/SuSiE_coloc_summary.tsv"
