
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(knitr)
library(tidyverse)
# library(purrr)


# ----------------------- Helper Functions to Implement ------------------------


#' Read the expression data "csv" file.
metadata <- read_csv("data/proj_metadata.csv")
head(metadata)[c(1:5)]


#' Function to read microarray expression data stored in a csv file. The
#' function should return a sample x gene tibble, with an extra column named
#' "subject_id" that contains the geo accession ids for each subject.
#'
#' @param filename (str): the file to read.
#'
#' @return
#' @export
#'
#' @examples expr_mat <- read_expression_table('example_intensity_data.csv')
read_expression_table <- function(filename) {
  read_data <- read.table(filename, sep = ' ', header = TRUE)
  subject_id <- colnames(read_data)
  transposed <- t(read_data) #transposing
  colnames(transposed) <- transposed[1,]
  trans_tibble <- as_tibble(transposed) #make it a tibble
  final_data <- cbind(subject_id[-1], mutate_all(trans_tibble[-1,], as.numeric))
  return (final_data)
}

expr_mat <- read_expression_table('data/example_intensity_data.csv')
tib_expr_mat <- as_tibble(expr_mat)
knitr::kable(head(tib_expr_mat[c(1:5)]))   #head(cbind(subject_id[-1], values))[c(1:5)]

#' Replaces all '.' in a string with '_'
#'
#' @param str String to operate upon.
#'
#' @return reformatted string.
#' @export
#'
#' @examples
#' period_to_underscore("foo.bar")
#' "foo_bar"
period_to_underscore <- function(str) {
  standardized <- str_replace_all(str, "[:punct:]", "_")
  return (standardized)
}
#period_to_underscore('foo.bar')

colnames(metadata) <-  lapply(colnames(metadata), period_to_underscore)
colnames(metadata)

# rename variables:
# Age_at_diagnosis to Age
# SixSubtypesClassification to Subtype
# normalizationcombatbatch to Batch

#' Rename and select specified columns.
#'
#' Function to rename Age_at_diagnosis, SixSubtypesClassification, and
#' normalizationcombatbatch columns to Age, Subtype, and Batch, respectively. A
#' subset of the data should be returned, containing only the Sex, Age, TNM_Stage,
#' Tumor_Location, geo_accession, KRAS_Mutation, Subtype, and Batch columns.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) renamed and subsetted tibble
#' @export
#'
#' @examples rename_and_select(metadata)
#' 
#' 
rename_and_select <- function(data) {
  metadata <- data %>% 
    rename(Age = Age_at_diagnosis,
           Subtype = SixSubtypesClassification,
           Batch = normalizationcombatbatch)
  #metadata %>%
   # select()
  return (metadata[,c('Sex', 'Age', 'TNM_Stage', 'Tumor_Location', 'geo_accession', 'KRAS_Mutation', 'Subtype', 
                      'Batch')])
}
selected_metadata <- rename_and_select(metadata)
colnames(selected_metadata)






#' Create new "Stage" column containing "stage " prefix.
#'
#' Creates a new column "Stage" with elements following a "stage x" format, where
#' `x` is the cancer stage data held in the existing TNM_Stage column. Stage
#' should have a factor data type.
#'
#' @param data  (tibble) metadata information for each sample
#'
#' @return (tibble) updated metadata with "Stage" column
#' @export
#'
#' @examples metadata <- stage_as_factor(metadata)
stage_as_factor <- function(data) {
  stage_data <- mutate(data,Stage = paste('stage',TNM_Stage))
  return (stage_data)
}


selected_metadata <- stage_as_factor(selected_metadata)
selected_metadata$Stage

#' Calculate age of samples from a specified sex.
#'
#' @param data (tibble) metadata information for each sample
#' @param sex (str) which sex to calculate mean age. Possible values are "M"
#' and "F"
#'
#' @return (float) mean age of specified samples
#' @export
#'
#' @examples mean_age_by_sex(metadata, "F")
#' 
#' rename_and_select <- function(data) {
#metadata <- metadata %>% 
#  rename(Age = Age_at_diagnosis,
#         Subtype = SixSubtypesClassification,
#         Batch = normalizationcombatbatch)
#metadata %>%
# select()
#return (metadata[,c('Sex', 'Age', 'TNM_Stage', 'Tumor_Location', 'geo_accession', 'KRAS_Mutation', 'Subtype', 
#selected_metadata <- rename_and_select(metadata)
#colnames(selected_metadata)
#' 
#' 
#' 
mean_age_by_sex <- function(data, sex) {
  if (sex == "M"){
    filtered <- filter(data, Sex == 'M')
    avg_age <- mean(filtered$Age)
  }else {
    filtered <- filter(data, Sex == 'F')
    avg_age <- mean(filtered$Age)
  }
  return (avg_age)
}


male_age <- mean_age_by_sex(selected_metadata, 'M') #mean male age
female_age <- mean_age_by_sex(selected_metadata, "F") #mean female age
print(paste0("The average male age is: ", round(male_age, digits = 2)))
print(paste0("The average female age is: ", round(female_age, digits = 2)))

#' Calculate average age of samples within each cancer stage. Stages should be
#' from the newly created "Stage" column.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) summarized tibble containing average age for all samples from
#' each stage.
#' @export
#'
#' @examples age_by_stage(data)
age_by_stage <- function(data) {
  grouped <- data %>% group_by(Stage) #grouped by stage column
  grouped_mean <- grouped %>% summarise(
    Age = mean(Age)
  )
  return (grouped_mean)
}

age_by_stage(selected_metadata)


#' Create a cross tabulated table for Subtype and Stage using dplyr methods.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) table where rows are the cancer stage of each sample, and the
#' columns are each cancer subtype. Elements represent the number of samples from
#' the corresponding stage and subtype. If no instances of a specific pair are
#' observed, a zero entry is expected.
#' @export
#'
#' @examples cross_tab <- dplyr_cross_tab(metadata)
subtype_stage_cross_tab <- function(data) {     #count the number of co-occurrences between cancer subtype and cancer progression
  grouped <- data %>% group_by(Stage)
  grouped_sum <- grouped %>% summarise(
    C3 = sum(Subtype == "C3"),
    C4 = sum(Subtype == "C4")
  )
  return (grouped_sum)
}


subtype_stage_cross_tab(selected_metadata)

#selected_metadata$Stage
#selected_metadata$SixSubtypesClassification
#' Summarize average expression and probe variability over expression matrix.
#'
#' @param exprs An (n x p) expression matrix, where n is the number of samples,
#' and p is the number of probes.
#'
#' @return A summarized tibble containing `main_exp`, `variance`, and `probe`
#' columns documenting average expression, probe variability, and probe ids,
#' respectively.
mean_per_col <- function(exprs) {
  vals <- exprs[,-1]
  mean_expr <- colMeans(vals)
  #  mean_exp = colMeans(exprs[,-1])
  return (mutate_all(mean_expr, as.numeric))
}
mean_exp <- as_tibble(mean_per_col(tib_expr_mat))
knitr::kable(head(mean_exp))

var_per_col <- function(exprs) {
  col_values <- exprs[,-1]
  variance <- sapply(col_values, var)
  return (mutate_all(variance, as.numeric))
}

variance <- as_tibble(var_per_col(tib_expr_mat))
knitr::kable(head(variance))

probes <- function(exprs) {
  return(colnames(exprs[1,-1]))
}

probe <- as_tibble(probes(tib_expr_mat))
names(probe)[1] <- "probe"
knitr::kable(head(probe))

total <- as_tibble(cbind(mean_exp, variance, probe)) 
head(total)

highest_exp <- dplyr::top_n(total, 10, mean_exp) %>%
  dplyr::arrange(dplyr::desc(mean_exp))

most_variable <- dplyr::top_n(total, 10, variance) %>%
  dplyr::arrange(dplyr::desc(variance))

print(paste0("The most highly expressed probes are ", paste(highest_exp$probe, collapse=', ')))
print(paste0("The most variable probes are ", paste(most_variable$probe, collapse=',')))

sub_meta <-
  dplyr::filter(selected_metadata, geo_accession %in% expr_mat$subject_id)

print(paste0("The average male age is: ", round(mean_age_by_sex(sub_meta, "M"), digits=2)))
print(paste0("The average female age is: ", round(mean_age_by_sex(sub_meta, "F"), digits=2)))
print(age_by_stage(sub_meta))
subtype_stage_cross_tab(sub_meta)
