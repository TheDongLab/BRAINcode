#!/bin/bash
#SBATCH --job-name=circ_ALS_FTD_lm
#SBATCH --output=/home/zw529/donglab/data/target_ALS/circ_ALS_FTD_lm.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/circ_ALS_FTD_lm.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

module load R

export METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
export CIRC_MATRIX="/home/zw529/donglab/data/target_ALS/QTL/circ_matrix.txt"
export GENE_BED="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gene.bed6"
export OUTDIR="/home/zw529/donglab/data/target_ALS"

Rscript - <<'EOF'

library(data.table)
library(lmerTest)
library(GenomicRanges)
library(IRanges)
library(parallel)

options(width=150)

meta <- fread(Sys.getenv("METADATA"))
circ_dt <- fread(Sys.getenv("CIRC_MATRIX"))

# --- CLEAN AND MERGE CATEGORIES ---
meta$subject_group <- trimws(gsub("[\r\n\t]+"," ",meta$subject_group))
meta$subject_group[meta$subject_group=="Non Neurological Control"] <- "Non-Neurological Control"
meta$subject_group[meta$subject_group=="ALS Spectrum MND, Other Neurological Diseases"] <- "ALS Spectrum MND, Other Neurological Disorders"

# --- FIX circ IDs (IMPORTANT) ---
circ_ids <- as.character(circ_dt$circ_id)
circ_dt[,circ_id:=NULL]

clean <- function(x) gsub("-","_",gsub(" ","",x))

meta$id <- clean(meta$externalsampleid)
colnames(circ_dt) <- clean(colnames(circ_dt))

samples <- intersect(meta$id,colnames(circ_dt))
cat("Matched samples:",length(samples),"\n")

meta <- as.data.frame(meta[match(samples,meta$id)])
circ <- as.data.frame(circ_dt[,..samples])
rownames(circ) <- circ_ids

stopifnot(identical(meta$id,colnames(circ)))

keep <- rowSums(circ > 0.001,na.rm=TRUE) >= 10
circ <- circ[keep,,drop=FALSE]

cat("Remaining circRNAs:",nrow(circ),"\n")

meta$sex <- factor(tolower(trimws(meta$sex)))
meta$tissue <- factor(trimws(meta$tissue))
meta$externalsubjectid <- factor(trimws(meta$externalsubjectid))

meta$age_at_death <- suppressWarnings(as.numeric(meta$age_at_death))
meta$post_mortem_interval_in_hours <- suppressWarnings(as.numeric(meta$post_mortem_interval_in_hours))

control_group <- "Non-Neurological Control"

als_groups <- c(
"ALS Spectrum MND",
"ALS Spectrum MND, Other Neurological Disorders",
"Alzheimer\u2019s Disease, Definite: CERAD criteria,FTLD-MND/MNI,Amyotrophic Lateral Sclerosis",
"Alzheimer's Disease, Definite: CERAD criteria,FTLD-MND/MNI,Amyotrophic Lateral Sclerosis",
"FTLD-MND/MNI,Amyotrophic Lateral Sclerosis",
"Pre-fALS, Other Neurological Disorders"
)

ftd_groups <- c(
"Alzheimer\u2019s Disease, Definite: CERAD criteria,FTLD-MND/MNI,Amyotrophic Lateral Sclerosis",
"Alzheimer's Disease, Definite: CERAD criteria,FTLD-MND/MNI,Amyotrophic Lateral Sclerosis",
"FTD, TDP43 subtype",
"FTLD",
"FTLD with tau positive inclusions",
"FTLD-MND/MNI,Amyotrophic Lateral Sclerosis",
"FTLD-TDP, Cerebrovascular disease"
)

# --- ANNOTATION ---
annotate_circs <- function(ids,bed_path){
  bed <- fread(bed_path,header=FALSE,sep="\t",select=1:4)
  setnames(bed,c("chrom","start0","end","name"))

  parts <- strsplit(bed$name,"___",fixed=TRUE)
  gene_name <- vapply(parts,function(x){
    if(length(x)>=3) x[3] else tail(x,1)
  },character(1))

  gene_gr <- GRanges(
    seqnames=bed$chrom,
    ranges=IRanges(start=as.integer(bed$start0)+1L,end=as.integer(bed$end)),
    gene_name=gene_name
  )

  p <- tstrsplit(ids,":",fixed=TRUE)
  r <- tstrsplit(p[[2]],"-",fixed=TRUE)

  ann <- data.table(
    circ_id=ids,
    chrom=p[[1]],
    start=as.integer(r[[1]]),
    end=as.integer(r[[2]]),
    strand=p[[3]],
    gene_name=NA_character_
  )

  valid <- !is.na(ann$start) & !is.na(ann$end)

  if(any(valid)){
    circ_gr <- GRanges(
      seqnames=ann$chrom[valid],
      ranges=IRanges(start=ann$start[valid],end=ann$end[valid])
    )

    hits <- suppressWarnings(findOverlaps(circ_gr,gene_gr,ignore.strand=TRUE))

    if(length(hits)>0){
      g <- split(mcols(gene_gr)$gene_name[subjectHits(hits)],queryHits(hits))
      collapsed <- vapply(g,function(x){
        paste(unique(x[nzchar(x)]),collapse=";")
      },character(1))
      idx <- which(valid)
      ann$gene_name[idx[as.integer(names(collapsed))]] <- collapsed
    }
  }

  ann
}

annotation <- annotate_circs(rownames(circ),Sys.getenv("GENE_BED"))

run_analysis <- function(label,groups,beta_col){

  group <- rep(NA_character_,nrow(meta))
  group[meta$subject_group==control_group] <- "Control"
  group[meta$subject_group %in% groups] <- label

  idx <- !is.na(group)
  meta_sub <- meta[idx,,drop=FALSE]
  circ_group <- circ[,idx,drop=FALSE]

  meta_sub$status <- factor(group[idx],levels=c("Control",label))

  base <- data.frame(
    status=meta_sub$status,
    age=meta_sub$age_at_death,
    sex=meta_sub$sex,
    PMI=meta_sub$post_mortem_interval_in_hours,
    tissue=meta_sub$tissue,
    subject=meta_sub$externalsubjectid
  )

  ok_cov <- complete.cases(base)
  base <- droplevels(base[ok_cov,,drop=FALSE])
  circ_model <- circ_group[,ok_cov,drop=FALSE]

  ids <- rownames(circ_model)

  results_list <- lapply(seq_along(ids),function(i){

    if(i %% 1000==0) cat("[",label,"] Completed",i,"/",length(ids),"LMs\n")

    y <- as.numeric(circ_model[i,])
    ok <- !is.na(y)

    if(sum(ok)<10) return(NULL)

    df <- base[ok,,drop=FALSE]
    df$circ <- y[ok]
    df <- droplevels(df)

    if(length(unique(df$status))<2) return(NULL)

    fit <- suppressMessages(suppressWarnings(
      tryCatch(
        lmer(circ ~ status + age + sex + PMI + tissue + (1|subject),data=df),
        error=function(e) NULL
      )
    ))

    if(is.null(fit)) return(NULL)

    coef <- summary(fit)$coefficients
    rown <- paste0("status",label)

    if(!(rown %in% rownames(coef))) return(NULL)

    beta <- coef[rown,"Estimate"]
    pval <- coef[rown,"Pr(>|t|)"]

    data.table(
      circ_id=ids[i],
      beta=beta,
      pvalue=pval,
      n=nrow(df),
      n_subjects=length(unique(df$subject)),
      mean_control=mean(df$circ[df$status=="Control"],na.rm=TRUE),
      mean_case=mean(df$circ[df$status==label],na.rm=TRUE),
      raw_difference=mean(df$circ[df$status==label],na.rm=TRUE) -
                     mean(df$circ[df$status=="Control"],na.rm=TRUE)
    )
  })

  results <- rbindlist(results_list,fill=TRUE)
  results[,FDR:=p.adjust(pvalue,"BH")]

  results <- merge(results,annotation,by="circ_id",all.x=TRUE,sort=FALSE)
  setnames(results,"beta",beta_col)

  results[,analysis:=paste0(label,"_vs_Control")]

  results <- results[order(-abs(get(beta_col)))]

  fwrite(results,file=file.path(Sys.getenv("OUTDIR"),paste0("circRNA_",label,"_adjusted_results.csv")))

  sig <- results[pvalue<=0.05]
  sig <- sig[order(-abs(get(beta_col)))]
  top50 <- head(sig,50)

  report <- capture.output({
    for(i in seq_len(nrow(top50))){
      row <- top50[i]
      cat("\n============================================================\n")
      cat(" circRNA ",i," OF ",nrow(top50),": ",row$circ_id,"\n",sep="")
      cat("============================================================\n\n")
      cat("Gene:",row$gene_name,"\n")
      cat("Coordinates:",row$chrom,":",row$start,"-",row$end,"\n",sep="")
      cat("\nAdjusted ",label," association:\n",sep="")
      cat(" ",beta_col,":",row[[beta_col]],"\n")
      cat(" pvalue:",row$pvalue,"\n")
      cat(" FDR:",row$FDR,"\n")
      cat(" Samples:",row$n,"\n")
      cat(" Subjects:",row$n_subjects,"\n\n")

      x <- as.numeric(circ_group[row$circ_id,])
      group_df <- data.frame(subject_group=meta_sub$subject_group,circ_value=x)
      group_mean <- aggregate(circ_value ~ subject_group,data=group_df,FUN=mean,na.rm=TRUE)
      colnames(group_mean)[2] <- "mean_circ_percentage"
      print(group_mean,row.names=FALSE,right=FALSE)
      cat("\n\n")
    }
  })

  writeLines(report,file.path(Sys.getenv("OUTDIR"),paste0("circRNA_",label,"_top50_report.txt")))

  list(label=label,n_models=nrow(results),n_sig=nrow(sig))
}

configs <- list(
  list(label="ALS",groups=als_groups,beta="beta_ALS"),
  list(label="FTD",groups=ftd_groups,beta="beta_FTD")
)

summaries <- mclapply(configs,function(cfg){
  run_analysis(cfg$label,cfg$groups,cfg$beta)
},mc.cores=2)

for(s in summaries){
  cat("\n==== ",s$label," SUMMARY ====\n")
  print(s)
}

EOF
