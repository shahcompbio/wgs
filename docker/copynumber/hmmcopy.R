#!/usr/bin/env Rscript
library(HMMcopy)

args <- commandArgs(TRUE)

tumour_copy <- args[1]
tumour_table <- args[2]
normal_copy <- args[3]
normal_table <- args[4]
segments <- args[5]
obj_file <- args[6]
sample_id <- args[7]
normal_table_out <- args[8]
tumour_table_out <- args[9]


m <- args[10] #,
mu <- args[11] #mu
kappa <- args[12]  #kappa
e <- args[13] #e
S <- args[14]
strength <- args[15]
lambda <- args[16]
nu <- args[17]
eta <- args[18]
gamma <- args[19]  

##################################################
# check the params
# the params should have the same length
# all should be set to NULL if one of them is
##################################################
checkNull = function(x) return(x == NULL)
if (m == 'NULL')
{
stopifnot(checkNull(m))
stopifnot(checkNull(mu))
stopifnot(checkNull(kappa))
stopifnot(checkNull(e))
stopifnot(checkNull(S))
stopifnot(checkNull(strength))
stopifnot(checkNull(lambda))
stopifnot(checkNull(nu))
stopifnot(checkNull(eta))
stopifnot(checkNull(gamma))
}



if (m != 'NULL')
{
m =  as.numeric(strsplit(m, ",")[[1]])
mu =  as.numeric(strsplit(mu, ",")[[1]])
kappa =  as.numeric(strsplit(kappa, ",")[[1]])
strength =  as.numeric(strsplit(strength, ",")[[1]])
lambda = as.numeric(strsplit(lambda, ",")[[1]])
e = as.numeric(strsplit(e, ",")[[1]])
S = as.numeric(strsplit(S, ",")[[1]])
nu = as.numeric(strsplit(nu, ",")[[1]])
eta = as.numeric(strsplit(eta, ",")[[1]])
gamma = as.numeric(strsplit(gamma, ",")[[1]])


a <- list(length(m), length(mu), length(kappa),
                              length(strength), length(lambda),
                              length(e), length(S), length(nu),
                              length(eta), length(gamma) )

#ensure all of these are same length
stopifnot( length(unique(list(length(m), length(mu), length(kappa),
                              length(strength), length(lambda),
                              length(e), length(S), length(nu),
                              length(eta), length(gamma) ))) == 1 )



params = data.frame(m=m, mu=mu, kappa=kappa,
                    strength=strength, lambda=lambda,
                    e=e,S=S,nu=nu, eta=eta, gamma=gamma)
} else {
params = NULL
}



load(tumour_copy)
tumour_copy <- infile_copy
remove(infile_copy)

tumour_table <- read.table(tumour_table, as.is=T, header=T, check.names=FALSE)

if (normal_copy != 'NULL')
{
    load(normal_copy)
    normal_copy <- infile_copy
    #if normal copy is not null then normalize the tumour by normal
    tumour_copy$copy <- tumour_copy$copy - normal_copy$copy

    normal_table <- read.table(normal_table, as.is=T, header=T, check.names=FALSE)
}

# call copy numbers using default parameters
tumour_segments <- HMMsegment(tumour_copy, param=params)

save(tumour_segments, file=paste(obj_file))

segs <- tumour_segments$segs
segs <- cbind(sample_id=sample_id,segs)
write.table(segs, file = segments, col.names = TRUE, row.names = FALSE, quote = FALSE, sep ="\t")

tumour_table$state <- tumour_segments$state
write.table(tumour_table, file = tumour_table_out, col.names = TRUE, row.names = FALSE, quote = FALSE, sep ="\t")

if (normal_copy != 'NULL')
{
	normal_table$state <- tumour_segments$state
	write.table(normal_table, file = normal_table_out, col.names = TRUE, row.names = FALSE, quote = FALSE, sep ="\t")
}





