library(stringr)
library(httr)

process_devco<-function(sample_id, polarity = c("Positive", "Negative"), devco_USR = "SA-PRD-TRACTION", devco_PW = "a-/g=`KS7huFw4s("){
  
  
  #sample_id = "maufiero-01-00203-3"
  #sample_id = "STAB-21-00459_S0173847"
  
  #devco_USR <- "SA-PRD-TRACTION"
  #devco_PW <-"a-/g=`KS7huFw4s(" 
  
  metadata1 = getDevCo2Sample(sample_id, devco_USR, devco_PW)
  
  # Format conversion:

  if (polarity == "Positive"){
    PEPMASS2 =  metadata1$exactMass + 1.0073
    ADDUCT2 = "M+H"
  }
  
  if (polarity == "Negative"){
    PEPMASS2 =  metadata1$exactMass - 1.0073
    ADDUCT2 = "M-H"
  }
  
  FORMULA2 = as.character(sapply(metadata1$empFormula, function(x)  gsub("[[:blank:]]", "", x)))
  
  metadata2 = cbind.data.frame(PEPMASS  = PEPMASS2, RT = metadata1$retentionTime,
                        IONMODE = polarity, ADDUCT = ADDUCT2, CHARGE = 1,
                        ID = metadata1$contAliasName, INSTRUMENT = "N/A", DATACOLLECTOR = "N/A",
                        SAMPLEID = metadata1$sampleId, ORGANISM = metadata1$projectAlias, JNJ_ID = metadata1$contAliasName,
                        SUBMITUSER = "N/A", SMILES = "N/A", INCHIKEY = "N/A", THEOMASS = PEPMASS2, FORMULA_FROM_SMILES = "N/A",
                        FORMULA = FORMULA2, COMPOUND = metadata1$contAliasName)
  return(metadata2)
}


listDevCo2Sample <- function(devco_USR, devco_PW){
  
  jSessionID = getSessionID(devco_USR, devco_PW)
  
  sampleInfo_response <- httr::GET("https://devco.jnj.com/api/sample", httr::set_cookies(JSESSIONID = jSessionID))
  #sampleInfo_response <- httr::GET("https://stg-devco.jnj.com/api/sample", httr::set_cookies(JSESSIONID = jSessionID))
  status_response <- httr::http_status(sampleInfo_response)
  if(status_response$category!= "Success") stop(status_response$message)
  sampleIfo <- jsonlite::fromJSON(httr::content(sampleInfo_response, as = "text", encoding = "UTF-8"))
  
  sampleids = unique(sampleIfo$sampleId)
  
  return(sampleids)
  
}

#########################
###Internal functions####
#########################

getSessionID <- function(devco_USR, devco_PW){
  
  response <- httr::POST(url= "https://devco.jnj.com/api/authentication", body = list(username = devco_USR, password = devco_PW, submit='Login'), encode = "form") 
  jSessionID_field <- response$headers[[1]] 
  jSessionID <- sub("JSESSIONID=", "", strsplit(jSessionID_field, split = ";")[[1]][1])
  
  if(jSessionID == "nosniff"){
    response <- httr::POST(url = "https://devco.jnj.com/api/authentication", body = list(username = devco_USR, password = devco_PW, submit='Login'), encode = "form")
    jSessionID_field <- response$headers[[1]]
    jSessionID <- sub("JSESSIONID=", "", strsplit(jSessionID_field, split = ";")[[1]][1])
  }
  return(jSessionID)
}

getDevCo2Info <- function(alias, devco_USR, devco_PW){
  
  jSessionID = getSessionID(devco_USR, devco_PW)
  aliasTypes <- c("J", "R", "T", "D")	
  
  # potential multiple aliases
  if(grepl("[[:blank:]]", alias)){
    aliases <- strsplit(alias, split = "[[:blank:]]")[[1]]
    idxAlias <- which.max(match(substr(aliases, 1, 1), aliasTypes))
    alias <- aliases[idxAlias]
  }
  
  aliasType <- substr(alias, 1, 1)
  
  if(!aliasType %in% aliasTypes){ 
    stop("The alias should start with 'J', 'R', 'T', 'D' (and of corresponding type).")
  }
  
  alias <- sub("-", "S", alias)
  compInfo_response <- httr::GET("https://devco.jnj.com/api/moleculeHistory", query = list(aliasName = alias), httr::set_cookies(JSESSIONID = jSessionID))
  #compInfo_response <- httr::GET("https://stg-devco.jnj.com/api/moleculeHistory", query = list(aliasName = alias), httr::set_cookies(JSESSIONID = jSessionID))
  status_response <- httr::http_status(compInfo_response)
  if(status_response$category!= "Success") stop(status_response$message)
  compInfo <- jsonlite::fromJSON(httr::content(compInfo_response, as = "text", encoding = "UTF-8"))
  compInfo[["statusCode"]] <- NULL
  compInfo[["moleculeId"]] <- NULL
  compInfo[["initialMoleculeId"]] <- NULL
  
  compImg_response<- httr::GET("https://devco.jnj.com/api/moleculeHistory/svgMoleculeByAliasName", httr::set_cookies(JSESSIONID = jSessionID), query = list(aliasName = alias))
  #compImg_response <- httr::GET("https://stg-devco.jnj.com/api/moleculeHistory/svgMoleculeByAliasName", httr::set_cookies(JSESSIONID = jSessionID), query = list(aliasName = alias))
  status_response2 <- httr::http_status(compImg_response)
  if(status_response2$category != "Success") stop(status_response2$message)
  compImgRaw <- compImg_response$content
  compImgBase64 <- RCurl::base64(compImgRaw)
  chemStrPath <- sprintf(paste0('<img src="%s" height="220"', "/>"), paste0("data:image/svg+xml;base64," , compImgBase64))
  
  return(list(cmpdInfo = compInfo, chemStrPath = chemStrPath))
  
}

getDevCo2Sample <- function(sample_id, devco_USR, devco_PW){

  jSessionID = getSessionID(devco_USR, devco_PW)
  
  sampleInfo_response <- httr::GET("https://devco.jnj.com/api/sample", query = list(sampleId = sample_id), httr::set_cookies(JSESSIONID = jSessionID))
  #sampleInfo_response <- httr::GET("https://stg-devco.jnj.com/api/sample", query = list(sampleId = sample_id), httr::set_cookies(JSESSIONID = jSessionID))
  status_response <- httr::http_status(sampleInfo_response)
  if(status_response$category!= "Success") stop(status_response$message)
  sampleIfo <- jsonlite::fromJSON(httr::content(sampleInfo_response, as = "text", encoding = "UTF-8"))

  return(sampleIfo)
  
}


