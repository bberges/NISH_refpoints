library(stockassessment)
library(msy)
library(tidyverse)
library(FLCore)
library(icesTAF)

rm(list=ls())

mkdir('report')

# ref points as at HAWG 2024 - check of consistency
name <- "NISH2024Final"
dataYear <- 2023

source('refpts_NISH_script.R')

# ref points for 2025 reissue - check the impact of SSB survey correction
name <- "NISH20242025update"
dataYear <- 2023

source('refpts_NISH_script.R')

# ref points for HAWG2025 - impact of both SSB survey correction and new 2024 data points
name <- "NISH2025"
dataYear <- 2024

source('refpts_NISH_script.R')
