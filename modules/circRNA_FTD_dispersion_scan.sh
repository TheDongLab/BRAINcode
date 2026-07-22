#!/bin/bash
#SBATCH --job-name=circ_FTD_lm
#SBATCH --output=/home/zw529/donglab/data/target_ALS/circ_FTD_lm.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/circ_FTD_lm.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

module load R

export METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
export CIRC_MATRIX="/home/zw529/donglab/data/target_ALS/QTL/circ_matrix.txt"

meta <- fread(Sys.getenv("METADATA"))
circ <- fread(Sys.getenv("CIRC_MATRIX"))

# --- CLEAN AND MERGE CATEGORIES ---
meta$subject_group <- trimws(gsub("[\r\n\t]+"," ",meta$subject_group))
meta$subject_group[meta$subject_group=="Non Neurological Control"] <- "Non-Neurological Control"
meta$subject_group[meta$subject_group=="ALS Spectrum MND, Other Neurological Diseases"] <- "ALS Spectrum MND, Other Neurological Disorders"
meta$subject_group[meta$subject_group=="Other Neurological Disorders"] <- "Other Neurological Disorders"

Rscript - <<'EOF'

library(data.table)
library(lmerTest)

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

keep <- rowSums(circ > 0.001,na.rm=TRUE) >= 10
circ <- circ[keep,]

cat("Remaining circRNAs:",nrow(circ),"\n")

meta$subject_group <- factor(meta$subject_group)
meta$subject_group <- relevel(
    meta$subject_group,
    ref="Non-Neurological Control"
)

meta$sex <- factor(meta$sex)
meta$tissue <- factor(meta$tissue)
meta$externalsubjectid <- factor(meta$externalsubjectid)

####################################################
# Run regression for every circRNA
####################################################
results <- lapply(
    seq_along(rownames(circ)),
    function(i){

        id <- rownames(circ)[i]

        if(i %% 1000 == 0)
            cat("Completed",i,"/",nrow(circ),"LMs\n")

        y <- as.numeric(circ[id,])

        df <- data.frame(
            circ=y,
            subject_group=meta$subject_group,
            age=meta$age_at_death,
            sex=meta$sex,
            PMI=meta$post_mortem_interval_in_hours,
            tissue=meta$tissue,
            subject=meta$externalsubjectid
        )

        df <- df[complete.cases(df),]

        if(length(unique(df$subject_group))<2)
            return(NULL)

        fit <- suppressWarnings(
            tryCatch(
                lmer(
                    circ ~ subject_group + age + sex + PMI + tissue + (1|subject),
                    data=df
                ),
                error=function(e) NULL
            )
        )

        if(is.null(fit))
            return(NULL)

        coef <- summary(fit)$coefficients

        group_rows <- grep(
            "^subject_group",
            rownames(coef),
            value=TRUE
        )

        if(length(group_rows)>0){

            data.frame(
                circ_id=id,
                comparison=group_rows,
                beta=coef[group_rows,"Estimate"],
                pvalue=coef[group_rows,"Pr(>|t|)"],
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
    order(abs(beta),decreasing=TRUE)
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

    cat("\nAdjusted FTD/FTLD association:\n")
    cat(
        " comparison:",
        top50$comparison[i],
        "\n"
    )
    cat(
        " beta:",
        top50$beta[i],
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
    "/home/zw529/donglab/data/target_ALS/circRNA_FTD_FTLD_adjusted_results.tsv",
    sep="\t"
)

EOF
