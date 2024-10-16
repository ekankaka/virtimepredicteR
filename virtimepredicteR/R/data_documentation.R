#' training_data
#' 
#' This is the data that was used to fit the Bayesian models.
#' Here, it is used by the get_virtime function to determine if the input APD
#' value is an outlier or not, compared to previous APD values from that gene region.
#' @format A data frane with 594 rows and 19 columns:
#' \describe{
#' \item{id}{participant identifiers.}
#' \item{gene_region}{An HIV gene or part of a gene for which a model was fitted}
#' \item{unique_sequences}{Number of non-indentical sequences used to fit the model for that gene}
#' \item{trimmed_alignment_width}{width of the fasta alignment after removing leading and trailing gaps}
#' \item{APD}{Average Pairwise Diversity calculated on all nucleotide positions}
#' \item{trimmed_alignment_width_codons}{width of the fasta alignment after removing leading and trailing gaps, the alignment only contains nucleotides at 3rd codon positions}
#' \item{APD_codons}{Average Pairwise Diversity calculated on 3rd codon positions}
#' \item{cohort}{Patient group where the HIV sequences were generated}
#' \item{sex}{Sex at Birth of the patient}
#' \item{subtype}{HIV subtype}
#' \item{on_ART_samples}{number samples from which the proviral sequences were obtained}
#' \item{nadir_cd4}{lowest CD4 ever recorded for this patient}
#' \item{logVL_before_ART}{log of the HIV viral load measured before ART}
#' \item{years_untreated}{number of years between estimated date of infection and ART initiation}
#' \item{years_on_ART}{number years between ART initiation and last sample }
#' \item{med_APD}{median APD for that gene region}
#' \item{mad_APD}{median absolute deviation for that gene region}
#' \item{z_robust}{robust Z-score for this APD value, based on all APD values for that gene region}
#' \item{APD_outlier}{1 if absolute value of Z is > 3}
#' }

#' best_models
#' 
#' This dataset includes the performance metrics of the models that explained > 10%
#' of the variability in viremic times.
#' @format A data frane with 9 rows and 8 columns:
#' \describe{
#' \item{gene_region}{An HIV gene or part of a gene for which a model was fitted}
#' \item{APD_version}{Whether Viremic Time was fitted on APD or APD calculated on third codon positions of the in-frame nucleotide alignment}
#' \item{model}{linear regression or piecewise linear regression}
#' \item{observations}{number of observations used to fit the model for that gene region}
#' \item{PD}{Penalized Deviance (performance metric)}
#' \item{MAE}{Mean Absolute Error (performance metric)}
#' \item{RSQ}{R-squared}
#' \item{RSQ_ADJ}{Adjusted R-squared}
#' }


#' example_GAG_P17.fasta
#' 
#' This fasta file can be used to demonstrate how to predict viremic time from sequences from the GAG_P17 region of HIV-1
#' The sequences from an anonymous patient are aligned to the HXB2 reference.
#' }
