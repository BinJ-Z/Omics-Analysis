#install ThermoRawFileParser： https://github.com/compomics/ThermoRawFileParser/releases
#install NET：https://builds.dotnet.microsoft.com/dotnet/Runtime/8.0.29/dotnet-runtime-8.0.29-osx-arm64.pkg
#setwd("/Users/bohe/Desktop/Metabolomics/Raw_data/MTBLS12932_LCMS_RAW")
#dir.create("mzML_positive", recursive = TRUE, showWarnings = FALSE)
#dir.create("mzML_negative", recursive = TRUE, showWarnings = FALSE)


setwd("/Users/bohe/Desktop/Metabolomics/Raw_data/MTBLS12932_LCMS_RAW")
dotnet <- "/usr/local/share/dotnet/dotnet"

# ThermoRawFileParser
parser <- "/Users/bohe/Desktop/Metabolomics/env/ThermoRawFileParser-v.2.0.0-dev-osx-arm64/ThermoRawFileParser.dll"


positive_raw <- list.files(
  "positive",
  pattern = "\\.raw$",
  full.names = TRUE,
  ignore.case = TRUE
)

negative_raw <- list.files(
  "negative",
  pattern = "\\.raw$",
  full.names = TRUE,
  ignore.case = TRUE
)

for (raw_file in positive_raw) {
  
  message("start：", basename(raw_file))
  
  command <- sprintf(
    '"%s" "%s" -i="%s" -o="%s" -f=2',
    dotnet,
    parser,
    normalizePath(raw_file),
    normalizePath("mzML_positive")
  )
  
  system(command)
}


for (raw_file in negative_raw) {
  
  message("start：", basename(raw_file))
  
  command <- sprintf(
    '"%s" "%s" -i="%s" -o="%s" -f=2',
    dotnet,
    parser,
    normalizePath(raw_file),
    normalizePath("mzML_negative")
  )
  
  system(command)
}

