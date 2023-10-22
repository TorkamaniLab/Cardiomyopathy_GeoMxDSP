# Spatial_DE_UKBB_RVA_replication
## 1. Sample selection - samples for each terminology
- working directory : /gpfs/group/home/salee/ukbb/UKBB_validation_script
- input files 
  1) HF_define_terminology.tsv # terminology definition
  2) med_ukbb.RData # UKBB phenotype data
- codes
```
srun --mem=80G --cpus-per-task=16 --exclusive --pty bash -i
module load R
R --no-save < Phenotype_sample_selection.R HF_define_terminology.tsv
```
- script : Phenotype_sample_selection.R
```
#srun --mem=80G --cpus-per-task=16 --exclusive --pty bash -i
#module load R 
#R --no-save < Phenotype_sample_selection.R HF_define_terminology.tsv


library(tidyverse)
library(data.table)

options(stringsAsFactors = FALSE)

library(doParallel)
# registerDoParallel(detectCores())
registerDoParallel(16)



################################
# Get system arguments
ehr_codings_fp <- as.character(commandArgs()[3])
# ehr_codings_fp <- "UKBB_phenotype_processing_Heart_Failure.tsv"



################################
# Load input data

obj <- load("med_ukbb.RData")
med_ukbb <- med_ukbb[med_ukbb$'f.34.0' != "", ] # drop individual with high missingness
dim(med_ukbb)
# med_ukbb <- head(med_ukbb)


ukbb_phenotyping_ehr_df <- read.csv(ehr_codings_fp,
                                    sep = '\t',
                                    na.strings=c("","NA"))
age_f6150 <- c(
    'f.3894.0', # Age_heart_attack_diagnosed
    'f.3627.0', # Age_angina_diagnosed
    'f.4056.0', # Age_stroke_diagnosed
    'f.2966.0'  # Age_high_blood_pressure_diagnosed
)
med_ukbb[age_f6150][med_ukbb[age_f6150] == ""] <- NA
med_ukbb['f.3894_3627_4056_2966.0'] <- apply(med_ukbb[age_f6150], 1, 
                                             function(x) {
                                                 toString(na.omit(x))
                                                 }
                                             )
dim(med_ukbb)
head(med_ukbb)



################################
# Debugging

# ukbb_phenotyping_ehr_df <- ukbb_phenotyping_ehr_df[which(ukbb_phenotyping_ehr_df$Name == "Atherosclerotic_cardiovascular_disease"), ]
# ukbb_phenotyping_ehr_df <- ukbb_phenotyping_ehr_df[which(ukbb_phenotyping_ehr_df$Name %in% 
#                                                          c("Atherosclerotic_cardiovascular_disease", "Stroke")), ]

# tlab_transformed_ukbb <- head(tlab_transformed_ukbb, 100)
# med_ukbb <- head(med_ukbb, 100)



################################
# Define helper functions

# get the earliest matched timestamp per EHR coding system
get_med <- function(row, name, lbls, event_fids, time_fids, code_systems, types){
    

    for (idx in 1:length(code_systems)){
        type <- types[idx]
        out_colnames <- paste0('f.', event_fids[idx], '.', time_fids[idx], '.', code_systems[idx], '.', type, "_", name)
        # print(out_colnames)
        if (!is.na(lbls[idx])){
            event_col <- names(row)[grepl(event_fids[idx], names(row), perl=TRUE)]
            event_list <- unlist(str_split(row[event_col], ", "))

            time_col <- names(row)[grepl(time_fids[idx], names(row), perl=TRUE)]
            time_list <- unlist(str_split(row[time_col], ", "))

            # target events by string startswith the pattern
            # remove dot in the mearning into ukbb code label
            target_pt <- paste0('^', str_replace_all(str_replace_all(lbls[idx], '\\.', ''), ', ', '|^'))
            match_event <- str_detect(event_list, target_pt)

            # print(icd10_date_list[target_icd10])
            # print(min(as.Date(icd10_date_list[target_icd10]), na.rm = TRUE))
            match_time <- time_list[match_event]
            # print(length(target_icd10_date))
            # date_icd10 <- names(row)[grepl("41280", names(row), perl=TRUE)]
            # print(icd10_list)

            if (length(match_time) > 0){
                if (type == "EHR"){
                    min_match_time <- as.character(min(as.Date(match_time)))
                } else if (type == "INTVW"){
                    min_match_time <- min(as.integer(match_time))
                } else {
                    print("Invalid health record data type.")
                }

                row[[out_colnames]] <- min_match_time
            } else {
                row[[out_colnames]] <- NA
            }
        } else {
            row[[out_colnames]] <- NA
        }
    }
    return(row)
}



# collect medical information per disease (line of define table) by EHR, interview and flag
collect_med <- function(row){
    
    out <- data.frame()
    lbls <- list()
    event_fids <- list()
    time_fids <- list()
    code_systems <- list()
    types <- list()
    
    name <- row[['Name']]
    print(paste0("Collecting ", name))
      for (type in c("EHR", "INTVW")){
        pt <- paste0("^", type, "\\.\\.")
        in_ehr_cols <- colnames(ukbb_phenotyping_ehr_df)[grepl(pt, colnames(ukbb_phenotyping_ehr_df), perl=TRUE)]

        # if (type %in% c("EHR", "INTVW")){
        # collect set of data resources per disease per type
        for (in_ehr_pair in in_ehr_cols){
            event <- unlist(str_split(in_ehr_pair, '\\.\\.'))[2]
            event_fid <- unlist(str_split(event, '\\.'))[2]
            code_sys <- unlist(str_split(event, '\\.'))[3]

            time <- unlist(str_split(in_ehr_pair, '\\.\\.'))[3]
            time_fid <- unlist(str_split(time, '\\.'))[2]

            lbls <- append(lbls, list(row[[in_ehr_pair]]))
            event_fids <- append(event_fids, list(event_fid))
            time_fids <- append(time_fids, list(time_fid))
            code_systems <- append(code_systems, list(code_sys))
            
            types <- append(types, list(type))
        }
    }
    out <- as.data.frame(t(apply(med_ukbb, 1, get_med, 
                         name = name, 
                         lbls = lbls, 
                         event_fids = event_fids, 
                         time_fids = time_fids, 
                         code_systems = code_systems,
                         types = types
                                )))
    out <- out[colnames(out)[grepl(name, colnames(out), perl=TRUE)]]
    # out[] <- lapply(out, as.Date) # didn't work
    
    pt <- paste0("^FLAG\\.\\.")
    in_flag_cols <- colnames(ukbb_phenotyping_ehr_df)[grepl(pt, colnames(ukbb_phenotyping_ehr_df), perl=TRUE)]
    flag_fid <- row[[in_flag_cols]]
    
    if (!is.na(flag_fid)){
        flag_df <- med_ukbb[colnames(med_ukbb)[grepl(flag_fid, colnames(med_ukbb), perl=TRUE)]]
        colnames(flag_df) <- paste0("f.", flag_fid, ".FLAG_", name)
        out <- cbind(out, flag_df)
    }
    
    return(out)
}



# Collapse target labels into the earlist timestamp of each disease outcome from multiple coding system in EHR (med_ukbb)
med_dfs <- foreach(name=ukbb_phenotyping_ehr_df[['Name']]) %dopar% {
    row <- ukbb_phenotyping_ehr_df[ukbb_phenotyping_ehr_df$Name == name, ]
    collect_med_df <- collect_med(row)
    out_list <- list()
    
    ##################
    for (type in c("EHR", "INTVW", "FLAG")){

        target_ehr_df <- collect_med_df[colnames(collect_med_df)[grepl(paste0(type, "_", name, "$"), colnames(collect_med_df), perl=TRUE)]]

        if (length(colnames(target_ehr_df)) != 0) {

            fin_col <- paste0(type, "_", name)

            if (type == "EHR"){
                print(paste0("Processing EHR_", name))
                target_ehr_df[] <- lapply(target_ehr_df, as.Date) # didn't work??
                # display only the earliest timestamp among each coding system
                target_ehr_df[fin_col] <- as.Date(apply(target_ehr_df, 1, function(x) {ifelse(all(is.na(x)), NA, na.omit(min(x, na.rm=TRUE)))}))
                out_list[[type]] <- target_ehr_df[fin_col]

            } else if (type == "INTVW"){
                print(paste0("Processing INTVW_", name))
                target_ehr_df$'f.34.Year_of_birth' <- med_ukbb$'f.34.0'

                target_ehr_df$'Year_of_attending_assessment_centre' <- format(as.Date(med_ukbb$'f.53.0'), "%Y")
                yr_cols <- colnames(target_ehr_df)[!colnames(target_ehr_df) %in% c("f.34.Year_of_birth", "Year_of_attending_assessment_centre")]
                target_ehr_df[] <- lapply(target_ehr_df, as.integer)
                target_ehr_df$'f.53.Date_of_attending_assessment_centre' <- as.Date(med_ukbb$'f.53.0')
                for (yr_col in yr_cols){
                    # convert age into year
                    target_ehr_df <- target_ehr_df %>% mutate_at(yr_col, 
                                                                 ~ ifelse(!!as.name(yr_col) < 200, 
                                                                          paste0(!!as.name(yr_col) + !!as.name('f.34.Year_of_birth')),
                                                                          .
                                                                         )) # TODO: improve to eliminate for-loop
                }

                # fill timestamp to the end of year so it won't mask the earlier EHR
                target_ehr_df[fin_col] <- as.integer(
                    apply(target_ehr_df[yr_cols], 
                          1, function(x) {
                              ifelse(all(is.na(x)), 
                                     NA, 
                                     na.omit(min(x, na.rm=TRUE))
                                    )
                              }
                          )
                )

                target_ehr_df[fin_col] = ifelse(is.na(target_ehr_df[[fin_col]]),
                    NA,
                    ifelse(
                        target_ehr_df[[fin_col]] == target_ehr_df[["Year_of_attending_assessment_centre"]], 
                        target_ehr_df[['f.53.Date_of_attending_assessment_centre']],
                        paste0(target_ehr_df[[fin_col]], '-12-31')
                    )
                )

                target_ehr_df[[fin_col]] <- as.Date(target_ehr_df[[fin_col]], "%Y-%m-%d")
                out_list[[type]] <- target_ehr_df[fin_col]

            } else if (type == "FLAG"){
                print(paste0("Processing FLAG_", name))
                target_ehr_df[] <- lapply(target_ehr_df, function(x) {as.Date(x, "%Y-%m-%d")})
                target_ehr_df[fin_col] <- as.Date(apply(target_ehr_df, 1, function(x) {ifelse(all(is.na(x)), NA, na.omit(min(x, na.rm=TRUE)))}))
                out_list[[type]] <- target_ehr_df[fin_col]
            } else {
                out_list[[type]] <- NULL
                print("Invalid type.")
            }
        }
    }

    out_list <- out_list[!sapply(out_list, is.null)]
    out_df <- do.call(cbind, out_list)

    return(out_df)
}
                                                        
# # remove NULL in the vector if fields missing for all instances
med_dfs <- med_dfs[!sapply(med_dfs, is.null)]

# # combine results across all fields
med_df <- do.call(cbind, med_dfs)
HF_med_ukbb <- cbind(med_ukbb['eid'], med_df) ## Object name!!!
gc(reset = TRUE)
                                                        
# save image
out_prefix <- gsub(".tsv", "", tail(strsplit(ehr_codings_fp, "/")[[1]], 1))
save(HF_med_ukbb, file = paste0(out_prefix, ".RData")) ## filename!!!
```
- output 
HF_define_terminology.RData

## 2. Sample selection - samples for each terminology
- script
Supplf5_UKBB_rarevariantanalysis_sample_selection.ipynb
```
# Data preparation

library(stringr)
library(tidyverse)

getwd()

## 1. Loading total cohort Rdata

obj1 <- load("med_ukbb.RData")

dim(med_ukbb)

## 2. Check the last date of events for each criterion

### 1) OPCS4

head(med_ukbb$f.41282.0)

f41282_vector <-unlist((str_split(med_ukbb$f.41282.0, ",")))

max(f41282_vector)
min(f41282_vector)

head(sort(f41282_vector))
tail(sort(f41282_vector))

### 2) OPCS3

f41283_vector <-unlist((str_split(med_ukbb$f.41283.0, ",")))

max(f41283_vector)

### 2) ICD10

f41280_vector <-unlist((str_split(med_ukbb$f.41280.0, ",")))

max(f41280_vector)

### 3) ICD9

f41281_vector <-unlist((str_split(med_ukbb$f.41281.0, ",")))

max(f41281_vector)

### 4) Death 

f40000_vector <-unlist((str_split(med_ukbb$f.40000.0, ",")))

max(f40000_vector)

## 3. Remove subject with no information

### Remove subject with no information

med_ukbb2 <- med_ukbb[med_ukbb$'f.34.0' != "", ]


dim(med_ukbb2)

## 5. Loading Heart failure phenotypes

# obj2<-load("SEL_Heart_Failure_ver2.RData")
obj2<-load("HF_define_terminology.RData")

colnames(HF_med_ukbb)

### Check the id indentity

identical(med_ukbb2$eid, HF_med_ukbb$eid)

### Add Birth info to Disease cohort

HF_med_ukbb$f.34.0 <- med_ukbb2$f.34.0

# Select Control cohort

## 1. Inclusion: Died after 70 years old


all_df <- HF_med_ukbb %>% mutate(year_death=as.numeric(format(EHR_All_death,"%Y")))


all_df <- all_df %>% mutate(age_death=year_death-as.numeric(f.34.0))


# died after 70 year old 
all_df_deathafter70 <- all_df %>% filter(age_death>70)

## 2. Exclusion

### 1) Excluding patients with any cardiac conditions other than simple HTN

names(all_df_deathafter70)

names(all_df_deathafter70)[53]
names(all_df_deathafter70)[54]

all_df_deathafter70_nocardiac <- all_df_deathafter70[rowSums(is.na(all_df_deathafter70[2:54]))==53,]

dim(all_df_deathafter70_nocardiac)

all_df_deathafter70_nocardiac_wes <- all_df_deathafter70_nocardiac 


dim(all_df_deathafter70_nocardiac_wes)

write.table(all_df_deathafter70_nocardiac_wes$eid, "control_deathafter70_nocardiac_total.ver20230728", col.names=FALSE, row.names=FALSE, sep ="\t")


# Select Heart Failure cohort

## 1) Inclusion: Any Heart Failures

# Inclusion
1) DCMP_strict (EHR)
2) HCMP_strict (EHR)
3) Ischemic_CMP (EHR)
4) HF_broad (EHR, INTVW)

hf_df <- HF_med_ukbb %>% filter(is.na(EHR_DCMP_strict)==FALSE 
                                 | is.na(EHR_HCMP_strict)==FALSE 
                         | is.na(EHR_Ischemic_CMP)==FALSE
                                | is.na(EHR_HF_broad)==FALSE 
                                | is.na(INTVW_HF_broad)==FALSE)


## 2) Exclusion: Age of diagnosis, Age of death, Other known causes of heart failure

# age of diagnosis
hf_df <- hf_df %>% mutate(year_hf=pmin(as.numeric(format(EHR_DCMP_strict,"%Y")), as.numeric(format(EHR_HCMP_strict,"%Y")), as.numeric(format(EHR_Ischemic_CMP,"%Y")), as.numeric(format(EHR_HF_broad, "%Y")), as.numeric(format(INTVW_HF_broad, "%Y")), na.rm = TRUE))



# example
hf_df %>% filter(eid==1094459)

hf_df <- hf_df %>% mutate(age_hf=(year_hf-as.numeric(f.34.0)))


head(hf_df %>% select(EHR_DCMP_strict, EHR_HCMP_strict, EHR_Ischemic_CMP, EHR_HF_broad, INTVW_HF_broad, f.34.0, year_hf, age_hf))

hf_df <- hf_df %>% mutate(year_death=as.numeric(format(EHR_All_death,"%Y")))


hf_df <- hf_df %>% mutate(age_death=year_death-as.numeric(f.34.0))


#### exclusion
hf_df_excall <- hf_df %>%
filter(
    #CMP2nd
    is.na(EHR_CMP_2nd) & is.na(INTVW_CMP_2nd)
    #valvular
    & is.na(EHR_VHD_rheumatic) & is.na(INTVW_VHD_rheumatic) & is.na(EHR_VHD_others) 
    & is.na(INTVW_VHD_others) & is.na(EHR_VHD_endocarditis) & is.na(INTVW_VHD_endocarditis)
    #myocarditis
    & is.na(EHR_Myocarditis) & is.na(INTVW_Myocarditis)
    #Congenital_HD
    & is.na(EHR_Congenital_HD) & is.na(INTVW_Congenital_HD)
    #Hypertension
    & is.na(EHR_HF_hypertensive) & is.na(INTVW_HF_hypertensive)
    #Renal failure
    & is.na(EHR_Renal_failur_sig) & is.na(INTVW_Renal_failur_sig)

)

dim(hf_df_excall)

hf_df_excall_wes <- hf_df_excall 
#%>% filter(eid %in% wes_list_v) # WES set (interim only)

# Age of diagnosis distribution
dim(hf_df_excall_wes)
dim(hf_df_excall_wes[hf_df_excall_wes$age_hf<=40,])
dim(hf_df_excall_wes[hf_df_excall_wes$age_hf<=50,])
dim(hf_df_excall_wes[hf_df_excall_wes$age_hf<=60,])
dim(hf_df_excall_wes[hf_df_excall_wes$age_hf<=70,])

hf_df_excall_wes_60 <-hf_df_excall_wes[hf_df_excall_wes$age_hf<=60,]
hf_v_excall_wes_60 <- hf_df_excall_wes_60$eid
length(hf_v_excall_wes_60)


hf_df_excall_wes_70 <-hf_df_excall_wes[hf_df_excall_wes$age_hf<=70,]
hf_v_excall_wes_70 <- hf_df_excall_wes_70$eid
length(hf_v_excall_wes_70)

write.table(hf_v_excall_wes_60, "hf_v_excall_total_60_ver20230728", col.names=FALSE, row.names=FALSE, sep ="\t")
write.table(hf_v_excall_wes_70, "hf_v_excall_total_70_ver20230728", col.names=FALSE, row.names=FALSE, sep ="\t")


# Select: CMP cohort

## Inclusion: Cardiomyopathy

# Inclusion
1) DCMP_strict (EHR)
2) HCMP_strict (EHR)
3) Ischemic_CMP (EHR)

cmp_df <-tlab_med_ukbb %>% filter(is.na(EHR_DCMP_strict)==FALSE 
                                 | is.na(EHR_HCMP_strict)==FALSE 
                         | is.na(EHR_Ischemic_CMP)==FALSE)


# age of diagnosis
cmp_df <- cmp_df %>% mutate(year_cmp=pmin(as.numeric(format(EHR_DCMP_strict,"%Y")), as.numeric(format(EHR_HCMP_strict,"%Y")), as.numeric(format(EHR_Ischemic_CMP,"%Y")), na.rm = TRUE))



cmp_df %>% filter(eid==1094459)

cmp_df <- cmp_df %>% mutate(age_cmp=(year_cmp-as.numeric(f.34.0)))


head(cmp_df %>% select(EHR_DCMP_strict, EHR_HCMP_strict, EHR_Ischemic_CMP, f.34.0, year_cmp, age_cmp))

cmp_df <- cmp_df %>% mutate(year_death=as.numeric(format(EHR_All_death,"%Y")))


cmp_df <- cmp_df %>% mutate(age_death=year_death-as.numeric(f.34.0))


## Exclusion

cmp_df_excall <- cmp_df %>%
filter(
    #CMP2nd
    is.na(EHR_CMP_2nd) & is.na(INTVW_CMP_2nd)
    #valvular
    & is.na(EHR_VHD_rheumatic) & is.na(INTVW_VHD_rheumatic) & is.na(EHR_VHD_others) 
    & is.na(INTVW_VHD_others) & is.na(EHR_VHD_endocarditis) & is.na(INTVW_VHD_endocarditis)
    #myocarditis
    & is.na(EHR_Myocarditis) & is.na(INTVW_Myocarditis)
    #Congenital_HD
    & is.na(EHR_Congenital_HD) & is.na(INTVW_Congenital_HD)
    #Hypertension
    & is.na(EHR_HF_hypertensive) & is.na(INTVW_HF_hypertensive)
    #Renal failure
    & is.na(EHR_Renal_failur_sig) & is.na(INTVW_Renal_failur_sig)

)

# WES list
cmp_df_excall_wes <- cmp_df_excall # %>% filter(eid %in% wes_list_v)

dim(cmp_df_excall_wes)

dim(cmp_df_excall_wes)
dim(cmp_df_excall_wes[cmp_df_excall_wes$age_cmp<=40,])
dim(cmp_df_excall_wes[cmp_df_excall_wes$age_cmp<=50,])
dim(cmp_df_excall_wes[cmp_df_excall_wes$age_cmp<=60,])
dim(cmp_df_excall_wes[cmp_df_excall_wes$age_cmp<=70,])

cmp_df_excall_wes_70 <-cmp_df_excall_wes[cmp_df_excall_wes$age_cmp<=70,]
cmp_df_excall_wes_60 <-cmp_df_excall_wes[cmp_df_excall_wes$age_cmp<=60,]


cmp_v_excall_wes_70 <- cmp_df_excall_wes_70$eid
cmp_v_excall_wes_60 <- cmp_df_excall_wes_60$eid


length(cmp_v_excall_wes_70)
length(cmp_v_excall_wes_60)

write.table(cmp_v_excall_wes_70, "cmp_v_excall_total_70_ver20230728", col.names=FALSE, row.names=FALSE, sep ="\t")
write.table(cmp_v_excall_wes_60, "cmp_v_excall_total_60_ver20230728", col.names=FALSE, row.names=FALSE, sep ="\t")


# Total cohort for WES sample selection (control + HF cohort)

total_cohort_hfctrl_20230728<-c(all_df_deathafter70_nocardiac_wes$eid, hf_v_excall_wes_70)

length(total_cohort_hfctrl_20230728)

write.table(total_cohort_hfctrl_20230728, "total_cohort_hfctrl_ver20230728", col.names=FALSE, row.names=FALSE, sep ="\t")




# Confirm eid and add phenotype data

df_base<-read.csv("total_baseline_ethnicity_20230617_participant.csv") # data download from RAP / birth, sex, ethnicity

head(df_base)

tail(df_base)

dim(df_base)

names(df_base)[1]

names(hf_df_excall_wes_60)

# Left join
hf_df_excall_wes_60_ethnicity <- merge(x=hf_df_excall_wes_60,y=df_base, 
          by="eid", all.x=TRUE)

head(hf_df_excall_wes_60_ethnicity  %>% select(eid, f.34.0, p31, p34), 20)

tail(hf_df_excall_wes_60_ethnicity  %>% select(eid, f.34.0, p31, p34), 20)

hf_df_excall_wes_60_ethnicity <- hf_df_excall_wes_60_ethnicity %>% mutate(diff = p34-as.integer(f.34.0))

hf_df_excall_wes_60_ethnicity[hf_df_excall_wes_60_ethnicity$diff!=0,]

hf_df_excall_wes_60_ethnicity[
    hf_df_excall_wes_60_ethnicity[["p34"]]!=as.integer(hf_df_excall_wes_60_ethnicity[["f.34.0"]]), c("eid", "f.34.0",
                                                                                                     "p34") ]

dim(hf_df_excall_wes_60_ethnicity)

is.na(hf_df_excall_wes_60_ethnicity$eid)

hf_df_excall_wes_60_ethnicity_naomit <- hf_df_excall_wes_60_ethnicity %>% drop_na(eid)


dim(hf_df_excall_wes_60_ethnicity_naomit)

identical(hf_df_excall_wes_60_ethnicity_naomit$f.34.0, hf_df_excall_wes_60_ethnicity_naomit$p34)

hf_df_excall_wes_60_ethnicity_naomit[
    (hf_df_excall_wes_60_ethnicity_naomit$f.34.0!=hf_df_excall_wes_60_ethnicity_naomit$p34),]

head(hf_df_excall_wes_60_ethnicity_naomit)

head(hf_df_excall_wes_60_ethnicity_naomit[["Year.of.birth"]])

head(as.integer(hf_df_excall_wes_60_ethnicity_naomit[["f.34.0"]]))
```

- output
control_deathafter70_nocardiac_total.ver20230728
hf_v_excall_total_60_ver20230728
hf_v_excall_total_70_ver20230728
cmp_v_excall_total_60_ver20230728
cmp_v_excall_total_70_ver20230728
total_cohort_hfctrl_ver20230728

## 3. Extract UKBB whole exome seq HF cohort from final UKBB WES (UKBB RAP)
- folder
/gpfs/home/salee/knih_dcmp/script/validate_ukbb/script_final
- codes
```
dx run_dx.sh
```
- script : run_dx.sh
```
#!/bin/bash

for j in {1..22} X Y;
do sh dx_step1_bcftools_select_hfsample.chr.dx $j
done
```
- script : dx_step1_bcftools_select_hfsample.chr.dx
```
#!/bin/bash

DIRin="project-GX2g4X8JQ4BVf4QjyB00PZkB:/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release"
DIR_list="project-GX2g4X8JQ4BVf4QjyB00PZkB:/Scripts"
DIR_out="project-GX2g4X8JQ4BVf4QjyB00PZkB:/Outputs"
list=total_cohort_hfctrl_ver20230728 # all patient control included
input=ukb23157_c${1}_*.vcf.gz


for i in $(dx ls "$DIRin"/$input); do
dx run swiss-army-knife \
-iin="$DIRin"/"$i" \
-iin="$DIR_list"/$list \
-icmd='bcftools view \
--threads 8 \
--force-samples \
-S "${in_name[1]}" \
"${in_name[0]}" \
-Oz -o "${in_prefix[0]}".filtered.vcf.gz' \
--instance-type "mem1_ssd1_v2_x8" \
--destination ${DIR_out} \
-y
done
```
- Outputs : Files are downloaded at RAP folder

## 4. Manipulate vcfs 

1) indexing - step02.maketbi.bash
```
#!/bin/bash

DIR=/opt/notebooks/Files   #where files are

cd $DIR

mkdir -p ./tmp

bcfindex() {
  id="$1"
  bcftools index -f --tbi $id
}

export -f bcfindex

ls | grep "tbi" | sed "s/.tbi//g" > tbi.list # when there are some files with tbi
ls *.vcf.gz | grep -v -f tbi.list | parallel --tmpdir ./tmp --plus bcfindex {}
```

2) concatenate them by chromosomes- step03.conc.bash
```
#!/bin/bash

DIR=/opt/notebooks/Files

cd $DIR

# Concatenate files

mkdir -p ./tmp

bcfconc() {
  id="$1"
  bcftools concat --threads 8 --allow-overlaps $(ls *vcf.gz | grep "c${id}_") -Oz -o ukbb_hf_ctrl_chr${id}.vcf.gz
  #bcftools index --tbi ukbb_hf_ctrl_chr${id}.vcf.gz 
}

export -f bcfconc

echo {1..22} X Y | tr " " "\n" | parallel --tmpdir ./tmp --plus bcfconc {}

bcfindex() {
  id="$1"
  bcftools index -f --tbi $id
}

export -f bcfindex

ls | grep "tbi" | sed "s/.tbi//g" > tbi.list # when there are some files with tbi
ls *.vcf.gz | grep -v -f tbi.list | parallel --tmpdir ./tmp --plus bcfindex {}
```

3) filter in WES target lesion, filter out blacklist regions, repeat regions and filter out by DP, GQ, and AB (het) - step04.filter.bash
```
#!/bin/bash

export DIR=/opt/notebooks/Files
export REF=/opt/notebooks/Annot/GRCh38.p13.genome.fa
export BED1=/opt/notebooks/Annot/xgen_plus_spikein.GRCh38.bed
# https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=3803
export BED2=/opt/notebooks/Annot/hg38-blacklist.v2.bed
# Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z
export BED3=/opt/notebooks/Annot/hg38_RepeatMasker_genome.bed
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1680210428_Kj6cjPxvTmN9mzO4rtKyxDMq9Y0s&clade=mammal&org=&db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName=

export cpu=8
cd $DIR

export "BCFTOOLS_PLUGINS=/opt/notebooks/Tools/bcftools-1.3.1/plugins"

bcffilter() {
	id=ukbb_hf_ctrl_chr${1}.vcf.gz
bcftools view --threads $cpu -R $BED1 ${id} | bcftools view --threads $cpu -T ^$BED2| bcftools view --threads $cpu -T ^${BED3}| bcftools norm --threads $cpu -m-both -f $REF | bcftools +fill-tags --threads $cpu -Oz -o $(echo ${id} | sed 's/vcf.gz/vaf.vcf.gz/') -- -t FORMAT/VAF 

bcftools index $cpu --tbi $(echo ${id} | sed 's/vcf.gz/vaf.vcf.gz/')

bcftools filter --threads $cpu $(echo ${id} | sed 's/vcf.gz/vaf.vcf.gz/') \  
 --set-GTs . \
 -e '(GT="het" & FMT/VAF<0.30) | (TYPE="snp" & FMT/DP<7) | (TYPE="indel" & FMT/DP<10) | (FMT/GQ<30)' \
 -Oz -o $(echo ${id} | sed 's/vcf.gz/exome.excblrep.norm.vafgqfiltered.vcf.gz/g')

}

export -f bcffilter

# --plus enables {%...} to remove a substring
echo {1..22} X Y | tr " " "\n" | parallel --tmpdir ./tmp --plus bcffilter {}
```

4) indexing - step05.maketbi.bash
```
#!/bin/bash

DIR=/opt/notebooks/Files

cd $DIR

mkdir -p ./tmp

tabixindex() {
  id="$1"
  tabix -p vcf $id
}

export -f tabixindex

ls *vafgqfiltered.vcf.gz | parallel --tmpdir ./tmp --plus tabixindex {}
```

5) filter out non-variant regions - step07.gatk.bash
```
#!/bin/bash

DIR=/opt/notebooks/Files

cd $DIR

mkdir -p ./tmp

gatk_select() {
  id="$1"
  /opt/notebooks/Tools/jdk-20.0.1/bin/java \
	-Xmx64g \
	-Djava.io.tmpdir=./tmp \
	-jar /opt/notebooks/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
	SelectVariants \
    	-V ukbb_hf_ctrl_chr"$id".exome.excblrep.norm.vafgqfiltered.vcf.gz \
    	--exclude-non-variants true \
    	-O ukbb_hf_ctrl_chr"$id".exome.excblrep.norm.vafgqgatkfilt.vcf.gz

}

export -f gatk_select

echo {1..22} X Y | tr " " "\n" | parallel --tmpdir ./tmp --plus gatk_select {}
```

6) Annotation - step08.annotation.bash
```
export DIR1=/opt/notebooks/Files
export cpu=8


annotation_annovarNbcf() {
id="$1"
annot=/opt/notebooks/annovar
${annot}/table_annovar.pl ${DIR1}/$id $annot/humandb/ -buildver hg38 -out ${DIR1}/$(echo  $id | sed "s/.vcf.gz//g") -remove -protocol refGene,gnomad_exome,gnomad_genome,gnomad211_exome,gnomad312_genome,clinvar_20221231 -operation g,f,f,f,f,f -nastring . -vcfinput -polish

INPUT5=$(echo $id | sed "s/.vcf.gz/.hg38_multianno.vcf/g")
bgzip $INPUT5 && tabix -p vcf $(echo $INPUT5 | sed 's/vcf/vcf.gz/g')
INPUT6=$(echo $INPUT5 | sed 's/vcf/vcf.gz/g')

#cd $DIR1
#run for the first run
if [ -f "cadd_indel_hdr.txt" ]; then
    echo "cadd_indel_hdr.txt exists."
else
echo -e '##INFO=<ID=PHRED,Number=1,Type=Float,Description="CADD_phred annotation provided by CADD">' > cadd_indel_hdr.txt
fi

if [ -f "revel_hdr.txt" ]; then
    echo "revel_hdr.txt exists."
else
echo -e '##INFO=<ID=REVEL,Number=1,Type=Float,Description="REVEL annotation provided by REVEL">' > revel_hdr.txt
fi

annot2=/opt/notebooks/cadd_revel
bcftools annotate --threads $cpu -a $annot2/whole_genome_SNVs.chr.tsv.gz -h cadd_indel_hdr.txt -c CHROM,POS,REF,ALT,-,INFO/PHRED $INPUT6 | bcftools annotate --threads $cpu -a $annot2/gnomad.genomes.r3.0.indel.chr.tsv.gz -h cadd_indel_hdr.txt -c CHROM,POS,REF,ALT,-,INFO/PHRED | bcftools annotate --threads $cpu -a $annot2/revel.hg38.chr.sorted.nopoint.header.tsv.gz -h revel_hdr.txt -c CHROM,POS,REF,ALT,INFO/REVEL -Oz -o $(echo $INPUT6 | sed "s/.vcf.gz/_caddrevel.vcf.gz/g")
## need revel and gnomad files from original site
# - Download REVEL files from https://sites.google.com/site/revelgenomics/downloads
# - Download CADD files (indels and snps) from https://cadd.gs.washington.edu/download

tabix -p vcf $(echo $INPUT6 | sed "s/.vcf.gz/_caddrevel.vcf.gz/g")

}

export -f annotation_annovarNbcf
ls ukbb_hf_ctrl_chr*.exome.excblrep.norm.vafgqgatkfilt.vcf.gz | parallel --plus annotation_annovarNbcf {}
```
7) Concatenate each chr file - step09.conc.bash
```
#!/bin/bash

export DIR1=/opt/notebooks/Files
export cpu=20

cd $DIR1

# Concatenate files

bcftools concat --thread $cpu --allow-overlaps $(ls | grep caddrevel.vcf.gz$) -Oz -o ukbb_hf_ctrl_500k_clinvar.vcf.gz
tabix -p vcf ukbb_hf_ctrl_500k_clinvar.vcf.gz

```

8) Subset Filtered by Frequency - step10_1.filteredbyfreq.bash

```
export DIR1=/opt/notebooks/Files
export input=ukbb_hf_ctrl_500k_clinvar.vcf.gz
export DIR_out=$DIR1/filtered_vcf_gnomADall
cd $DIR1
export cpu=8

#!/bin/bash

filterbyfreq () {
        id="$1"
        mkdir -p $DIR_out
        maf=${id}

        bcftools view --threads $cpu $input  -e "(INFO/gnomAD_exome_ALL>$maf | INFO/gnomAD_genome_ALL>$maf)" -Oz -o $DIR_out/$(echo $input | sed "s/vcf.gz/MAF${maf}.vcf.gz/")
        tabix -p vcf  $DIR_out/$(echo $input | sed "s/vcf.gz/MAF${maf}.vcf.gz/")
}

export -f filterbyfreq

echo 0.01 0.001 0.0001 0.00001 | tr " " "\n" | parallel --plus filterbyfreq {}
```

9) Subset Filtered by deleteriousness (lof, revel and cadd) - step10_2.filteredbylofrevelcadd.bash
```
#!/bin/bash

export DIR=/opt/notebooks/Files
export input=ukbb_hf_ctrl_500k_clinvar.vcf.gz
export DIR_out=$DIR/filtered_vcf_gnomADall
export cpu=8
cd $DIR


filterbyfunc () {
        id="$1"

if [[ ${id} -eq 1 ]] # LOF only
then bcftools view --threads $cpu $input -i 'INFO/ExonicFunc.refGene="frameshift_deletion" | INFO/ExonicFunc.refGene="frameshift_insertion" | INFO/ExonicFunc.refGene="stopgain" | INFO/ExonicFunc.refGene="startloss" | INFO/ExonicFunc.refGene="stoploss" | INFO/Func.refGene="splicing"' -Oz -o $DIR_out/$(echo $input | sed 's/vcf.gz/lofonly.vcf.gz/')
tabix -p vcf $DIR_out/$(echo $input | sed 's/vcf.gz/lofonly.vcf.gz/')
fi

if [[ ${id} -eq 2 ]]  # Missense PHRED>20 | REVEL>0.5
 then bcftools view --threads $cpu $input  -i 'INFO/ExonicFunc.refGene="frameshift_deletion" | INFO/ExonicFunc.refGene="frameshift_insertion" | INFO/ExonicFunc.refGene="stopgain" | INFO/ExonicFunc.refGene="startloss" | INFO/ExonicFunc.refGene="stoploss" | INFO/Func.refGene="splicing" | PHRED>20 | REVEL>0.5' -Oz -o $DIR_out/$(echo $input | sed 's/vcf.gz/lofmissensecadd20orrevel5.vcf.gz/')
tabix -p vcf $DIR_out/$(echo $input | sed 's/vcf.gz/lofmissensecadd20orrevel5.vcf.gz/')
fi

if [[ ${id} -eq 3 ]] # Missense PHRED>20 & REVEL>0.5
 then bcftools view --threads $cpu $input  -i 'INFO/ExonicFunc.refGene="frameshift_deletion" | INFO/ExonicFunc.refGene="frameshift_insertion" | INFO/ExonicFunc.refGene="stopgain" | INFO/ExonicFunc.refGene="startloss" | INFO/ExonicFunc.refGene="stoploss" | INFO/Func.refGene="splicing" |(PHRED>20 & REVEL>0.5)' -Oz -o $DIR_out/$(echo $input | sed 's/vcf.gz/lofmissensecadd20andrevel5.vcf.gz/')
tabix -p vcf $DIR_out/$(echo $input | sed 's/vcf.gz/lofmissensecadd20andrevel5.vcf.gz/')
fi

if [[ ${id} -eq 4 ]]  # Missense PHRED>25 | REVEL>0.75
 then bcftools view --threads $cpu $input  -i 'INFO/ExonicFunc.refGene="frameshift_deletion" | INFO/ExonicFunc.refGene="frameshift_insertion" | INFO/ExonicFunc.refGene="stopgain" | INFO/ExonicFunc.refGene="startloss" | INFO/ExonicFunc.refGene="stoploss" | INFO/Func.refGene="splicing" | PHRED>25 | REVEL>0.75' -Oz -o $DIR_out/$(echo $input | sed 's/vcf.gz/lofmissensecadd25orrevel75.vcf.gz/')
tabix -p vcf $DIR_out/$(echo $input | sed 's/vcf.gz/lofmissensecadd25orrevel75.vcf.gz/')
fi

if [[ ${id} -eq 5 ]]  # Missense PHRED>25 & REVEL>0.75
 then bcftools view --threads $cpu $input  -i 'INFO/ExonicFunc.refGene="frameshift_deletion" | INFO/ExonicFunc.refGene="frameshift_insertion" | INFO/ExonicFunc.refGene="stopgain" | INFO/ExonicFunc.refGene="startloss" | INFO/ExonicFunc.refGene="stoploss" | INFO/Func.refGene="splicing" |(PHRED>25 & REVEL>0.75)' -Oz -o $DIR_out/$(echo $input | sed 's/vcf.gz/lofmissensecadd25andrevel75.vcf.gz/')
tabix -p vcf $DIR_out/$(echo $input | sed 's/vcf.gz/lofmissensecadd25andrevel75.vcf.gz/')
fi
}


export -f filterbyfunc

echo {1..5} | tr " " "\n" | parallel --plus filterbyfunc {}
```
10) Filtering by both frequency and deleteriousness - step10_3.bcftoolsisec.bash
```
export DIR=/opt/notebooks/Files
export input=ukbb_hf_ctrl_500k_clinvar.vcf.gz
export DIR_out=$DIR/filtered_vcf_gnomAD
export cpu=8
cd $DIR_out

isec () {
lof="$1"
for maf in $(ls *MAF*.vcf.gz);
do
DIR_out2=isec_filtered_vcf3
mkdir -p $DIR_out/$DIR_out2
bcftools isec -n=2 -w1 --threads $cpu $lof $maf -Oz -o $DIR_out2/$(echo $lof | sed "s/vcf.gz/$(echo $maf | sed 's/ukbb_hf_ctrl_500k_clinvar.targetgene.MAF0.//g')/g")
tabix -p vcf $DIR_out2/$(echo $lof | sed "s/vcf.gz/$(echo $maf | sed 's/ukbb_hf_ctrl_500k_clinvar.targetgene.MAF0.//g')/g") ;
done

}

export -f isec
ls *lof*.vcf.gz | parallel --plus isec {}

```
11) Make list of known pathogenic or likely pathogenic variants in known cardiomyopathy genes - step11.listknowncmpPLP.bash
```
#!/bin/bash

export DIR1=/opt/notebooks/Files
export cpu=20

cd $DIR1

input1=ukbb_hf_ctrl_500k_clinvar.vcf.gz

# make exclusion list - samples with variations in cmp related gene (only "pathogenic /likely pathogenic")
bcftools query -f '[%CHROM\t%POS\t%INFO/Gene.refGene\t%INFO/CLNSIG\t%GT\t%SAMPLE\n]' $input1 | grep -E "0/1|1/1" | grep -E  "Pathogenic|Likely_pathogenic" | grep -wf CMP_gene_all_exceptNoonan_addDMD.nospace.list > ukbb_hf_ctrl_500k_clinvar.list_full
cat ukbb_hf_ctrl_500k_clinvar.list_full | cut -f6 | sort | uniq > ukbb_hf_ctrl_500k_clinvar.list_sample 
 # CMP_gene_all_exceptNoonan_addDMD.nospace.list
:ACTC1,ACTN2,BAG3,CSRP3,DES,DSC2,DSG2,DSP,FHL1,FLNC,GAA,GLA,JPH2,JUP,LAMP2,LMNA,MYBPC3,MYH7,MYL2,MYL3,NEXN,PKP2,PLN,PRKAG2,RBM20,SCN5A,TCAP,TNNC1,TNNI3,TNNT2,TPM1,TTN,VCL,DMD
 # DCMP genes in Jordan E, Peterson L, Ai T, Asatryan B, Bronicki L, Brown E, Celeghin R, Edwards M, Fan J, Ingles J, James CA, Jarinova O, Johnson R, Judge DP, Lahrouchi N, Lekanne Deprez RH, Lumbers RT, Mazzarotto F, Medeiros Domingo A, Miller RL, Morales A, Murray B, Peters S, Pilichou K, Protonotarios A, Semsarian C, Shah P, Syrris P, Thaxton C, van Tintelen JP, Walsh R, Wang J, Ware J, Hershberger RE. Evidence-Based Assessment of Genes in Dilated Cardiomyopathy. Circulation. 2021 Jul 6;144(1):7-19. doi: 10.1161/CIRCULATIONAHA.120.053033. Epub 2021 May 5. PMID: 33947203; PMCID: PMC8247549.
 # Genes related to other cardiomyopathies in Whiffin N, Walsh R, Govind R, Edwards M, Ahmad M, Zhang X, Tayal U, Buchan R, Midwinter W, Wilk AE, Najgebauer H, Francis C, Wilkinson S, Monk T, Brett L, O'Regan DP, Prasad SK, Morris-Rosendahl DJ, Barton PJR, Edwards E, Ware JS, Cook SA. CardioClassifier: disease- and gene-specific computational decision support for clinical genome interpretation. Genet Med. 2018 Oct;20(10):1246-1254. doi: 10.1038/gim.2017.258. Epub 2018 Jan 25. PMID: 29369293; PMCID: PMC6558251.
 # Excluding Noonan cases, Including DMD
cat ukbb_hf_ctrl_500k_clinvar.list_full | awk '!(/0\/1/ && /GAA/)' > ukbb_hf_ctrl_500k_clinvar.noGAAhet.list_full
cat ukbb_hf_ctrl_500k_clinvar.noGAAhet.list_full | cut -f6 | sort | uniq > ukbb_hf_ctrl_500k_clinvar.noGAAhet.list_sample
 # excluding heterozygous case for GAA
```

## 5. Select the cohort and run rarevariant analysis
1) Select cohort and run rare variant analysis - step12.rarevariantanalysis.hfctrl.target.bash 
```
#!/bin/bash

export DIR=/opt/notebooks/isec_filtered_vcf
export cpu=8
cd $DIR

rarevariant () {
input="$1"

#list of samples to include
export dir_list=/opt/notebooks/List
export list_case=$dir_list/hf_v_excall_total_60_ver20230728 # cmp_v_excall_total_60_ver20230728 for cmp only
export list_control=$dir_list/control_deathafter70_nocardiac_total.ver20230728
export name_insert=hf60vsctrl_clv_noGAA_target # cmp60vsctrl_clv_noGAA_target for cmp only

export list=$dir_list/tmp_total_${name_insert}
export list_ex=$dir_list/ukbb_hf_ctrl_500k_clinvar.noGAAhet.list_sample

export dir_out=$DIR/plink_file_${name_insert}

export targetgene=$dir_list/targetgene.list
#from Spatialtranscriptomics analysis : PTPRM,CRIP3,CYP2J2,TSPAN9,CASQ1,UBE2H,ASB18,TAX1BP3,HRH2,PFKFB2,APOBEC3A_B,CALD1

if [ -f "$list" ]; then
    echo "$list exists."
else
    cat $list_case $list_control > $list
fi

mkdir -p $dir_out

vcftools --gzvcf $input \
--max-missing 0.85 \
--not-chr X \
--not-chr Y \
--not-chr M \
--keep $list \
--remove $list_ex \
--hwe 0.00001 \
--recode --recode-INFO-all \
--out $(echo $input | sed s/.vcf.gz//).${name_insert}_missing

export input1=$(echo $input | sed 's/.vcf.gz//').${name_insert}_missing.recode.vcf  #input1_2


##PLINK PCA

/opt/notebooks/Tools/plink2/1.90b3.42/bin/plink --vcf $input1 --double-id --allow-extra-chr --pca --out $dir_out/$(echo $input1 | sed 's/_missing.recode.vcf//') \
--vcf-half-call missing

## PLINK files
/opt/notebooks/Tools/plink2/1.90b3.42/bin/plink --vcf $input1 --double-id --allow-extra-chr --make-bed --out $dir_out/$(echo $input1 | sed 's/_missing.recode.vcf//') --vcf-half-call missing

cd $dir_out

## PCA plot
echo $input1 | sed 's/_missing.recode.vcf//' | parallel --plus 'Rscript /opt/notebooks/step42.R_pcaplot.R {}'

## add phenotype
input2=$(echo $input1 | sed 's/_missing.recode.vcf/.fam/')
list1=$list_case
list2=$list_control

awk 'NR==FNR{a[$0];next} $1 in a {$6=2}1' $list1 $input2 > ${input2}.tmp1
awk 'NR==FNR{a[$0];next} $1 in a {$6=1}1' $list2 ${input2}.tmp1 > ${input2}.tmp2

paste -d " " ${input2}.tmp2 $(echo $input2 | sed 's/fam/eigenvec/g') > ${input2}.tmp3
cut -d " " -f1,2,3,4,5,6,9,10 ${input2}.tmp3 > ${input2}.tmp4
mv ${input2}.tmp4 ${input2}

## make setID

cd $DIR
input3=$(echo $input1 | sed 's/_missing.recode.vcf//')
bcftools query -f '%INFO/Gene.refGene\t%CHROM:%POS:%REF:%ALT\n' ${input3}_missing.recode.vcf | grep -wf $targetgene > $dir_out/${input3}.SetID
# bcftools query -f '%INFO/Gene.refGene\t%CHROM:%POS:%REF:%ALT\n' ${input3}_missing.recode.vcf > $dir_out/${input3}.SetID
# when no target genes. 
cd $dir_out

cat $(echo $input1 | sed 's/_missing.recode.vcf/.bim/') | awk -F'\t' -vOFS='\t' '{$2="chr"$1":"$4":"$6":"$5}1' > $(echo $input1 | sed 's/_missing.recode.vcf/.bim.tmp/')
mv $(echo $input1 | sed 's/_missing.recode.vcf/.bim.tmp/') $(echo $input1 | sed 's/_missing.recode.vcf/.bim/')


# run skat with Rscript
Rscript /opt/notebooks/step48_3.skat_binaryall_PCadj_ver20230617noresample.R $input3

}


export -f rarevariant

ls *.vcf.gz | parallel --plus rarevariant {}
```
- step42.R_pcaplot.R
```
#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
path = getwd()
setwd(path)
#install.packages("ggplot2")
library(ggplot2)

eigen <- read.table(paste0(args, ".eigenvec"))
eigenvalues <- read.table(paste0(args,".eigenval"))

# proportion of variation explained by each vector
eigen_percent <- round((eigenvalues/(sum(eigenvalues))*100), 2)

# pcaplot
pcaplot_1 <-ggplot(eigen, aes(V3, V4, label = V1), size =0.1) +
geom_point() +
xlab(paste0("PC1 (",eigen_percent[1,1],"%)")) +
ylab(paste0("PC2 (",eigen_percent[2,1],"%)")) +
ggtitle(paste(args, 'PCAplot_QC_PCA1-2')) +
geom_text(aes(label=V1), hjust=0.7, vjust=2, size =3)


pcaplot_2 <- ggplot(eigen, aes(V5, V6, label = V1), size =0.1) +
geom_point() +
xlab(paste0("PC3 (",eigen_percent[3,1],"%)")) +
ylab(paste0("PC4 (",eigen_percent[4,1],"%)")) +
ggtitle(paste(args, 'PCAplot_QC_PCA3-4')) +
geom_text(aes(label=V1), hjust=0.7, vjust=2, size =3)


ggsave(paste0(args, ".pca_1-2_plot.png"), pcaplot_1, scale=0.5)
ggsave(paste0(args, ".pca_3-4_plot.png"), pcaplot_2, scale=0.5)

```

- step48_3.skat_binaryall_PCadj_ver20230617noresample.R
```
#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
path=getwd()
setwd(path)

#install.packages("SKAT")
library(SKAT)
library(ggplot2)

Input <- args[1]

#input files
File.bed <-paste0('./', Input, '.bed')
File.bim <-paste0('./', Input, '.bim') #'./dcmp_LPandP_id.bim'
File.fam <-paste0('./', Input, '.fam') #'./dcmp_LPandP_id.fam'
File.SetID <-paste0('./', Input, '.SetID') #'./dcmp.acmgLPandP.SetID'

#output files
File.SSD <-paste0('./', Input, '.SSD')
File.Info <-paste0('./', Input, '.SSD.info')

#generate SSD
Generate_SSD_SetID(File.bed, File.bim, File.fam, File.SetID, File.SSD, File.Info)

#
FAM <- Read_Plink_FAM(File.fam, Is.binary=TRUE)
y <- FAM$Phenotype
SSD.INFO<-Open_SSD(File.SSD, File.Info)
colnames(FAM)[c(7,8)] <- c("PC1", "PC2")

X1 <- FAM$PC1
X2 <- FAM$PC2

obj<-SKAT_Null_Model(y ~ X1 + X2, out_type="D")
#, n.Resampling=1000, type.Resampling="bootstrap")

# SKAT
out.skat<-SKATBinary.SSD.All(SSD.INFO, obj, method="SKAT")
write.table(out.skat$results, file=paste0('./', Input, 'skat_result.tsv'), quote=FALSE, sep='\t', col.names = NA)

# SKAT-O
out.skato<-SKATBinary.SSD.All(SSD.INFO, obj, method="SKATO")
write.table(out.skato$results, file=paste0('./', Input, 'skato_result.tsv'), quote=FALSE, sep='\t', col.names = NA)


# Burden
out.burden<-SKATBinary.SSD.All(SSD.INFO, obj, method="Burden")
write.table(out.burden$results, file=paste0('./', Input, 'burden_result.tsv'), quote=FALSE, sep='\t', col.names = NA)

##QQ plot
#QQplot_skat <-
png(paste0('./', 'QQplot_skat_', Input, '.png'))
QQPlot_Adj(out.skat$results$P.value, out.skat$results$MAP)

#QQplot_skato <-
png(paste0('./', 'QQplot_skato_', Input, '.png'))
QQPlot_Adj(out.skato$results$P.value, out.skato$results$MAP)

#QQplot_burden <-
png(paste0('./', 'QQplot_burden_', Input, '.png'))
QQPlot_Adj(out.burden$results$P.value, out.burden$results$MAP)

dev.off()

```

## 6. Summarize the results and adjust p-value with Bonferroni correction
```R
library(tidyverse)
system('ls *.tsv', intern = TRUE) # confirm the location of the files

df_hfprog_spatial <- data.frame(matrix(ncol = 12, nrow = 0))
for (i in system('ls *.tsv', intern = TRUE)){
    print(i)
    df<-read.csv(i, header=TRUE, sep="\t")
            if (dim(df)[1]>0) {

    df$filename <-i
    df$P.adj_bonf <- p.adjust(df$P.value, method = "bonferroni", n = length(df$P.value))
    df <- df %>% dplyr::select(X, SetID, P.value, P.adj_bonf, 4:10)
    df_hfprog_spatial <- rbind(df_hfprog_spatial, df)
                }
}

df_hfprog_spatial[df_hfprog_spatial$P.adj_bonf<0.05,]

write.csv(df_hfprog_spatial, "ukbb_hf_ctrl_500k_cmp60vsctrl_clv_noGAA_target.csv")   ## or ukbb_hf_ctrl_500k_hf60vsctrl_clv_noGAA_target.csv
```

