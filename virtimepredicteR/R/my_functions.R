#' div function
#' This function takes an in-frame alignment and calculates Diversity at each position of the alignment
#'
#' @export
div <- function(fas, variants = c("A","C","G","T","-"), errorthreshold = 0.01){
  if (!(is(fas, "DNAStringSet"))){
    stop("fas should be a DNAStringSet. Load the fasta file in R using the readDNAStringSet function from the Biostrings package")
  }

  # check if fasta file has 2 or more sequences, excluding refs!
  fas_noRefs = fas[!grepl("HXB2|Ref", names(fas))]
  # width of the alignment
  width = unique(width(fas_noRefs))

  if (length(fas_noRefs)>=2){
    d_accross = sapply(1:width, function(j){
      # variants at position j
      posj = strex::str_elem(as.character(fas_noRefs),j)
      # frequency of each variant at position j
      freqj = sapply(variants, function(x){sum(posj == x)/length(posj)})
      # frequency of major variant at position j
      freqj_max = freqj[freqj==max(freqj)]
      # check if sum of the frequency of minor variants is above the cutoff xc
      indicator_variable = ifelse( (1-freqj_max) > errorthreshold, 1, 0)

      # diversity contribution of each variant at position j
      d = sum(sapply(freqj, function(x){ x * (1-x)}))
      # result at position j
      resultj = indicator_variable*d
    })

  } else {
    d_accross = rep(NA, width) # Do not calculate diversity for < 2 sequences
  }
  return(d_accross)
}


#' trim_start_end function
#' This function locates the first and last non-gap characters in an alignment
#'
#' @export
trim_start_end <- function(seq) {
  seq_char <- as.character(seq)

  # Find the first non-gap character positions
  start_pos <- max(regexpr("[^-]", seq_char))

  # Find the last non-gap character positions
  reversed_strings <- sapply(seq_char, function(x) paste(rev(strsplit(x, NULL)[[1]]), collapse=""))
  end_pos <- nchar(seq_char)[[1]] - max(regexpr("[^-]", reversed_strings)) + 1

  return(c(start_pos, end_pos))
}

#' get_APD function
#' This function takes an in-frame alignment and calculates APD on all nucleotides,
#' as well as APD on the 3rd codon position
#'
#' @export
get_APD <- function(fas, variants = c("A","C","G","T","-"), errorthreshold = 0.01){
  # fas should be a DNAStringSet object
  if (!(is(fas, "DNAStringSet"))){
    stop("fas should be a DNAStringSet object. Load the fasta file in R using the readDNAStringSet function from the Biostrings package")
  }

  # remove reference sequences if any
  if (grepl("ref", tolower(names(fas))) | grepl("hxb2", lower(names(fas)))){
    fas_NOREFS = fas[!grepl("ref", tolower(names(fas))) &
                       !grepl("hxb2", tolower(names(fas)))]
  } else {
    fas_NOREFS = fas
  }

  # number of sequences
  nseq = length(fas_NOREFS)

  # diversity at each position in the alignment
  d = div(fas = fas_NOREFS, variants, errorthreshold)
  d = unlist(d)

  ## remove diversity for leading and trailing terminal gaps
  d_start = trim_start_end(fas_NOREFS)[1]
  d_end = trim_start_end(fas_NOREFS)[2]
  d_trimmed = d[d_start:d_end]

  # diversity at third codon positions of the alignment
  d3_start <- d_start + (3 - (d_start %% 3))
  d3_end <- d_end - (d_end %% 3)
  d3_trimmed = d[seq(d3_start, d3_end, by = 3)]

  # width of the trimmed alignment for APD and for APD at third codon positions
  width_d = length(d_trimmed)
  width_d3 = length(d3_trimmed)

  # APD
  APD = mean(d_trimmed)
  APD3 = mean(d3_trimmed)
  return(c(unique_sequences = nseq,
           trimmed_alignment_width=width_d,
           APD=APD,
           trimmed_alignment_width_codons=width_d3,
           APD_codons=APD3))
}

#' get_viremictime function
#'
#' This function takes an "input" dataframe with the columns "id", "gene_region", "APD".
#' For each eligible gene region, it returns an estimate of viremic time for each APD value.
#'
#' @export
get_viremictime <- function(APD_value, gene_region){
  # assert that gene regions are eligible
  # i.e., among those for which we have a good predictive model with adj_R2 > 10%
  eligible_gene_regions = unique(best_models$gene_region)
  if ( !(gene_region %in% best_models$gene_region)){
    stop(paste("INELIGIBLE GENE REGION INCLUDED. CAN ONLY ESTIMATE VIREMIC TIMES FOR:",
               paste(eligibe_gene_regions, collapse = ", ")))
  } else {
    # read in the training data
    training_data$x = training_data$APD
    training_data$y = training_data$years_untreated

    # subset of training dataset for this gene region
    train.s = subset(train, gene_region == gene_region & !is.na(x) & !is.na(y))

    # check if xnew is not an oultier
    med_APD = median(train.s$APD)
    mad_APD = mad(train.s$APD)
    z_robust = (APD - med_APD) / mad_APD
    APD_outlier = abs(z_robust) > 3

    # load appropriate MCMC object
    idx = match(g,best_models$gene_region)
    mcmcfile = paste0("../data/",
                      best_models$gene_region[idx], "_",
                      best_models$APD_version[idx], "_",
                      gsub(" ","", best_models$model[idx]),
                      ".rds")

    # load dataset
    mod_sim = readRDS(mcmcfile)
    mod_csim = as.mcmc(do.call(rbind, mod_sim))

    # samples for coefficients in linear regression
    a_samples = mod_csim[,c('a')]
    b_samples = mod_csim[,c('b')]
    sig_samples = mod_csim[,c("sig")]

    # predicted mean for new x,
    mu_new = a_samples + b_samples * APD_value

    # Generate predicted y_new samples
    y_new_samples <- rnorm(1e4, mean = mu_new, sd = 1 / sqrt(sig_samples))

    # point estimate and uncertainty around it
    predicted_y <- mean(y_new_samples)
    HPDI_lwr <- HPDinterval(as.mcmc(y_new_samples), prob = 0.95)[1]
    HPDI_upr <- HPDinterval(as.mcmc(y_new_samples), prob = 0.95)[2]

    # model used
    model = best_models$model[idx]

  }
  return(c(gene_region = gene_region,
           APD=APD,
           predicted_viremic_time = predicted_y,
           HPDI_lwr = HPDI_lwr,
           HPDI_upr = HPDI_upr,
           possible_APD_outlier = as.numeric(APD_outlier),
           model = model))
}

