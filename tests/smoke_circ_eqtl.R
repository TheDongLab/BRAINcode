# Run in an isolated temporary directory; no cohort analyses are submitted.
args <- commandArgs(TRUE)
modules <- normalizePath(args[1])
library(data.table)
root <- tempfile('circ-eqtl-smoke-'); dir.create(root)
set.seed(42)
N <- 60
ids <- paste0('D',seq_len(N))
g <- rep(0:2,20); status <- rep(0:1,each=30)
y <- 2*g+3*g*status+rnorm(N,sd=.1)
write_matrix <- function(x,id,path) {
 d<-as.data.frame(x);names(d)<-ids;d<-cbind(id=id,d);fwrite(d,path,sep='\t')
}
write_matrix(matrix(g,1), 'rs1', file.path(root,'snp.txt'))
write_matrix(matrix(qnorm((rank(y)-.5)/N),1), 'chr17:45500000-45501000:+',file.path(root,'expr.txt'))
write_matrix(rbind(rnorm(N), status), c('age','is_als'),file.path(root,'cov.txt'))
fwrite(data.table(geneid='chr17:45500000-45501000:+',chr='chr17',left=45500001,right=45501000),file.path(root,'gene_location.txt'),sep='\t')
fwrite(data.table(snpid='rs1',chr='chr17',pos=45500500),file.path(root,'snp_location.txt'),sep='\t')
fwrite(data.table(geneid='chr17:45500000-45501000:+',host_gene='HOST',host_transcript='TX'),file.path(root,'circ_annotation.tsv'),sep='\t')
writeLines('FID IID PAT MAT SEX PHENOTYPE 17:45500500:A:G_G',file.path(root,'raw.txt'))
Sys.setenv(EQTL_MODE='circ', EQTL_RAW_FILE=file.path(root,'raw.txt'))
run <- function(script, arguments) {
 code <- system2(file.path(R.home('bin'),'Rscript'),c(shQuote(file.path(modules,script)),shQuote(arguments)))
 if(code != 0) stop(script,' failed: ',code)
}
for (interaction in c(FALSE,TRUE)) {
 suffix <- if(interaction) '_LINEAR_CROSS' else ''
 prefix<-file.path(root,if(interaction) 'interaction' else 'standard')
 run(paste0('_eQTL',suffix,'.R'),c(file.path(root,c('snp.txt','expr.txt','cov.txt')),prefix,file.path(root,c('gene_location.txt','snp_location.txt'))))
 run(paste0('_eQTL_postprocess',suffix,'.R'),c(paste0(prefix,'.cis.txt'),file.path(root,c('snp_location.txt','gene_location.txt')),prefix,'0.05','1000000'))
 result<-fread(paste0(prefix,'.full_annotated.txt'));stopifnot(nrow(result)==1,result$geneid=='chr17:45500000-45501000:+',result$host_gene=='HOST')
 run(paste0('_eQTL_boxplot',suffix,'.R'),c(paste0(prefix,'.top_for_boxplot.txt'),file.path(root,c('snp.txt','expr.txt','cov.txt','snp_location.txt')),root,if(interaction)'Interaction' else 'Standard'))
 run('_eQTL_manhattan.R',c(paste0(prefix,'.full_annotated.txt'),paste0(prefix,'.lead_snps.txt'),prefix,'0.05',if(interaction)'interaction' else 'standard',paste0(file.path(root,'standard'),'.full_annotated.txt')))
}
# Optional regression against the live, unmodified association engines.
if (length(args) >= 2) {
    original <- normalizePath(args[2])
    Sys.setenv(EQTL_MODE='gene')
    for (suffix in c('', '_LINEAR_CROSS')) {
        script <- paste0('_eQTL',suffix,'.R')
        results <- list()
        for (i in 1:2) {
            prefix <- file.path(root,paste0('regression',i,suffix))
            source <- if (i==1) modules else original
            status <- system2(file.path(R.home('bin'),'Rscript'),c(shQuote(file.path(source,script)),shQuote(c(file.path(root,c('snp.txt','expr.txt','cov.txt')),prefix,file.path(root,c('gene_location.txt','snp_location.txt'))))))
            stopifnot(status==0)
            results[[i]] <- fread(paste0(prefix,'.cis.txt'))
        }
        stopifnot(identical(results[[1]],results[[2]]))
    }
    cat('DEFAULT MODE REGRESSION PASSED\n')
}
cat('SMOKE PASSED:',root,'\n')
