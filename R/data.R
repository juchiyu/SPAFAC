#' USA Deaths in 2018 by Age and Cause
#'
#' This dataset contains information on the number of deaths in the USA in 2018,
#' categorized by age and cause of death. The causes of death are classified
#' according to the International Statistical Classification of Diseases and
#' Related Health Problems (ICD-10). Ages of death are grouped in 5-year intervals,
#' starting from 1 year old, with all deaths above 100 years grouped into a "100+"
#' category. Deaths before the age of 1 are excluded due to their predominance
#' in the perinatal causes category, which could overshadow other effects.
#' This exclusion is based on the observation that perinatal causes were rare in other age groups,
#' with 57 out of 118 cases occurring before 1 year old and the remaining 61 cases
#' spread across the other 21 age groups. The dataset represents a 21x19
#' contingency table of counts, showing the number of deaths for each cause
#' at each age range.
#'
#' @format A data frame with 399 rows and 3 columns:
#' \describe{
#'   \item{AgeGroup}{Age group of death, grouped in 5-year intervals, with a "100+" category.}
#'   \item{CauseOfDeath}{Cause of death, categorized according to the ICD-10.}
#'   \item{NumberOfDeaths}{The number of deaths for each cause at each age range.}
#' }
#' @source International Statistical Classification of Diseases and Related Health Problems (ICD-10).
#' @references
#' \url{https://www.who.int/classifications/icd/icdonlineversions/en/}
"example1_sCA"

#' Punctuation Marks in Literary Works from Project Gutenberg
#'
#' This dataset compiles the counts of punctuation marks in literary works
#' by authors of different times and origins. The texts were sourced from
#' the Gutenberg Project and processed using the R `gutenbergr` package.
#' The dataset includes authors from France (N = 63), the United Kingdom (N = 61),
#' and the United States (N = 36) spanning various periods (refer to Table 1 for details).
#'
#' The dataset focuses on the frequency of nine punctuation marks: comma (,),
#' period (.), question mark (?), exclamation mark (!), colon (:), semicolon (;),
#' apostrophe ('), English quotation marks (“ ”), French quotation marks (« »),
#' dashes (-), and M-dashes (—). Works that were translations were excluded to
#' ensure the dataset only contains books written in the authors' original languages.
#'
#' @format A data frame where each row represents a literary work and includes
#' the following columns for punctuation mark counts: `Comma`, `Period`, `QuestionMark`,
#' `ExclamationMark`, `Colon`, `Semicolon`, `Apostrophe`, `EnglishQuotes`, `FrenchQuotes`,
#' `Dash`, `MDash`, along with metadata about the work such as `Author`, `Country`, and `Period`.
#' @source Project Gutenberg, processed using the `gutenbergr` R package.
#' @references
#' Robinson, D. (2021). gutenbergr: Download and Process Public Domain Works from Project Gutenberg. R package version 0.2.0.
#' \url{https://CRAN.R-project.org/package=gutenbergr}
"example2_sDiSCA"

#' Chinese Version of the Independent and Interdependent Self Scale (C-IISS) Dataset
#'
#' This dataset contains responses from 130 undergraduate students (77 females and 53 males,
#' mean age = 19.49, SD = 1.52) from National Cheng Kung University to the Chinese version
#' of the Independent and Interdependent Self Scale (C-IISS). The scale, developed by
#' Lu et al. (2007), comprises forty-two 7-point Likert scale items (1 = strongly disagree,
#' 7 = strongly agree), with 21 items measuring independence (awareness and value of oneself
#' as an individual) and 21 measuring interdependence (valuation and actions based on one's
#' cohort). Participants provided written informed consent and were compensated with NTD 120.
#'
#' Prior to analysis with sMCA, responses for each item were binned into categories of
#' comparable sizes to mitigate the influence of event rarity in MCA (and sMCA) results.
#' The association patterns between these items are illustrated in Supplementary Figures
#' (refer to Figure S1 for binning and Figure S2 for association patterns).
#'
#' @format A data frame with 130 rows and 42 columns representing the Likert scale responses
#' for each item, along with demographic information of the participants.
#' @source Lu et al. (2007). Development of the Chinese version of the Independent and Interdependent
#' Self Scale. Journal reference.
#' @references
#' Lu, L. (2007). Developing the Chinese version of the Independent and Interdependent Self Scale.
#' Journal Title, Volume(Issue), Page range.
#' Supplementary Material: Figures S1 (Binning of Responses) and S2 (Association Patterns).
"example3_sMCA"

#' Math Assessment for College Students (MACS) Dataset
#'
#' This dataset comprises responses to the Math Assessment for College Students
#' (MACS) from undergraduate students enrolled in an introductory psychology
#' statistics course at an urban public college in the northeast United States.
#' The MACS, a 30-item paper-and-pencil test covering five content domains
#' (Algebra, Arithmetic, Categorization and Ranges, Decimals/Fractions/Percentages,
#' and Visual Understanding), was administered across five semesters to 460 participants.
#'
#' The MACS test measures basic mathematics skills and was completed in approximately
#' 40 minutes during the first week of the semester. Responses were binary scored
#' (0 = incorrect, 1 = correct) with no partial credit. The course was computationally
#' based, involving manual statistical tests and software programs. Academic performance
#' was assessed through three exams per semester, with scores averaged and categorized
#' into letter grades (A, B, C, D, F). For analysis purposes, grades A and B were grouped
#' together, as were D and F, resulting in three groups (AB, C, DF).
#'
#' Informed consent was obtained under an IRB-approved protocol, and no compensation
#' was provided. Demographic and academic performance data were also collected.
#' The dataset is used to examine relationships between basic mathematics skills,
#' demographic data, and academic performance.
#'
#' @format A data frame with 460 rows and variables for each MACS item, demographic
#' information, and academic performance scores. Each MACS item is labeled according
#' to its domain and number, with additional columns for overall scores and categorized
#' academic performance.
#' @source MACS test as described in Rabin (2018).
#' @references
#' Rabin, L.A. (2018). Development and Validation of the Math Assessment for College
#' Students. Journal Title, Volume(Issue), Page range.
"example4_sDiMCA"



