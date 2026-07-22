#!/bin/bash
#SBATCH --job-name=circ_ALS_lm
#SBATCH --output=/home/zw529/donglab/data/target_ALS/circ_ALS_lm.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/circ_ALS_lm.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

module load R

export METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
export CIRC_MATRIX="/home/zw529/donglab/data/target_ALS/QTL/circ_matrix.txt"

Rscript - <<'EOF'

library(data.table)
library(lme4)

meta <- fread(Sys.getenv("METADATA"))
circ <- fread(Sys.getenv("CIRC_MATRIX"))

rownames(circ) <- circ$circ_id
circ$circ_id <- NULL

clean <- function(x)
    gsub("-","_",gsub(" ","",x))

meta$id <- clean(meta$externalsampleid)
colnames(circ) <- clean(colnames(circ))

samples <- intersect(meta$id,colnames(circ))

cat("Matched samples:",length(samples),"\n")

meta <- meta[match(samples,id),]
circ <- as.data.frame(circ[,..samples])

keep <- rowSums(circ > 0.001 ,na.rm=TRUE) >= 5
circ <- circ[keep,]

cat("Remaining circRNAs:",nrow(circ),"\n")

meta$is_als <- factor(
    grepl("ALS",meta$subject_group,ignore.case=TRUE),
    levels=c(FALSE,TRUE),
    labels=c("Non_ALS","ALS")
)

meta$subject_group <- factor(meta$subject_group)
meta$sex <- factor(meta$sex)
meta$tissue <- factor(meta$tissue)
meta$externalsubjectid <- factor(meta$externalsubjectid)

####################################################
# Run regression for every circRNA
####################################################
results <- lapply(
    rownames(circ),
    function(id){

        y <- as.numeric(circ[id,])

        df <- data.frame(
            circ=y,
            is_als=meta$is_als,
            age=meta$age_at_death,
            sex=meta$sex,
            PMI=meta$post_mortem_interval_in_hours,
            tissue=meta$tissue,
            subject_group=meta$subject_group,
            subject=meta$externalsubjectid
        )

        df <- df[complete.cases(df),]

        if(length(unique(df$is_als))<2)
            return(NULL)

        fit <- lmer(
            circ ~ is_als + subject_group + age + sex + PMI + tissue + (1|subject),
            data=df
        )

        coef <- summary(fit)$coefficients

        if("is_alsALS" %in% rownames(coef)){

            data.frame(
                circ_id=id,
                beta_ALS=coef["is_alsALS","Estimate"],
                pvalue=coef["is_alsALS","Pr(>|t|)"],
                n=nrow(df)
            )

        } else NULL
    }
)

results <- rbindlist(results)
results$FDR <- p.adjust(
    results$pvalue,
    method="BH"
)

results <- results[
    order(abs(beta_ALS),decreasing=TRUE)
]

top50 <- head(results,50)

####################################################
# Print detailed output
####################################################
for(i in seq_len(nrow(top50))){
    id <- top50$circ_id[i]

    cat("\n")
    cat("============================================================\n")
    cat(" circRNA ",i," OF ",nrow(top50),": ",id,"\n",sep="")
    cat("============================================================\n")
    cat("\nAdjusted ALS association:\n")
    cat(
        " beta_ALS:",
        top50$beta_ALS[i],
        "\n"
    )
    cat(
        " pvalue:",
        top50$pvalue[i],
        "\n"
    )
    cat(
        " FDR:",
        top50$FDR[i],
        "\n"
    )
    cat(
        " Samples:",
        top50$n[i],
        "\n\n"
    )
    
    ################################################
    # Raw group means
    ################################################
    x <- as.numeric(circ[id,])

    group_df <- data.frame(
        subject_group=meta$subject_group,
        circ_value=x
    )
    group_mean <- aggregate(
        circ_value ~ subject_group,
        data=group_df,
        FUN=mean,
        na.rm=TRUE
    )
    colnames(group_mean)[2] <- "mean_circ_percentage"
    options(width=150)

    print(
        group_mean,
        row.names=FALSE,
        right=FALSE
    )
    cat("\n\n")
}

####################################################
# Save ranked table
####################################################
fwrite(
    results,
    "/home/zw529/donglab/data/target_ALS/circRNA_ALS_adjusted_results.tsv",
    sep="\t"
)

EOF
