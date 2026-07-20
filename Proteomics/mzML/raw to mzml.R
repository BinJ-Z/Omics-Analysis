
#path
setwd("/Users/bohe/Desktop/Proteomics")

project_dir <- "/Users/bohe/Desktop/Proteomics"

raw_dir <- "/Users/bohe/Desktop/Proteomics/Raw_data"

mzml_dir <- "/Users/bohe/Desktop/Proteomics/mzML"


# .NET
dotnet <- "/usr/local/share/dotnet/dotnet"


# ThermoRawFileParser
parser <- "/Users/bohe/Desktop/Metabolomics/env/ThermoRawFileParser-v.2.0.0-dev-osx-arm64/ThermoRawFileParser.dll"



raw_files <- list.files(
  raw_dir,
  pattern = "\\.raw$",
  full.names = TRUE,
  ignore.case = TRUE
)


# RAW → mzML


for (raw_file in raw_files) {
  
  message("start：", basename(raw_file))
  
  system2(
    command = dotnet,
    args = c(
      shQuote(parser),
      paste0("-i=", shQuote(raw_file)),
      paste0("-o=", shQuote(mzml_dir)),
      "-f=2"
    )
  )
  
}