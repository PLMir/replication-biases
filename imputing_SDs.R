# This is the first stage : inputing SDs that are missing to have complete data
# Collecting data stage (mostly using the metagear package)
library(metagear)

# importing csv file 
path <- file.path("MA_infants.csv")
MA <- read.csv(path)

# Subsetting into studies 

FBl <- subset(MA, seed_ID == "FBl")
NDi <- subset(MA, seed_ID == "NDis")
wWS <- subset(MA, seed_ID == "wordWSeg")
sWS <- subset(MA, seed_ID == "speechWSeg")
HTS <- subset(MA, seed_ID == "HTScram")
sIS <- subset(MA, seed_ID == "samInfS")
pIS <- subset(MA, seed_ID == "popInfS")
VST <- subset(MA, seed_ID == "VisStats")
LTS <- subset(MA, seed_ID == "LTScram")
IACT <- subset(MA, seed_ID == "IntAct")
GRS <- subset(MA, seed_ID == "GrSoc")

# NDi impute

NDi <- impute_SD(NDi, "SD_1", "x_1")
NDi <- impute_SD(NDi, "SD_2", "x_2")

# HTS impute
HTS <- impute_SD(HTS, "SD_1", "x_1")
HTS <- impute_SD(HTS, "SD_2", "x_2")


# IACT impute
IACT <- impute_SD(IACT, "SD_1", "x_1")
IACT <- impute_SD(IACT, "SD_2", "x_2")

# Binding back together
MA_imp <- rbind(FBl, NDi, wWS, sWS, HTS, sIS, pIS, VST, LTS, IACT, GRS)
write.csv(MA_imp, file = "MA_SD.csv", row.names = FALSE)
