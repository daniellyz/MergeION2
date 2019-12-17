#' Create additional metadata for library generation
#'
#' Function used by library_generator() to create adducts and multiple charged features
#' @export
#'
metadata_adduct_editor<-function(ref, adducts = c("Default","M+H","M+Na","M+K","M+NH4","M-H","M+Cl"), max.charge = 1){

  if (!("Default" %in% adducts)){

    ref_ids = unique(ref$ID)
    new_metadata = c()
    charges = 1:max.charge

    for (id in ref_ids){
        valid = which(ref_ids==id)
        sub_ref = ref[valid,] # metadata of the id
        new_pep = c()
        new_ions = c()
        new_adducts = c()
        new_charges = c()

    # Calculate neutral masses:
        if (sub_ref$ADDUCT[1]=="M+H"){NM = sub_ref$PEPMASS[1] - 1.007276}
        if (sub_ref$ADDUCT[1]=="M+Na"){NM = sub_ref$PEPMASS[1] - 22.989221}
        if (sub_ref$ADDUCT[1]=="M+K"){NM = sub_ref$PEPMASS[1] - 38.963158}
        if (sub_ref$ADDUCT[1]=="M+NH4"){NM = sub_ref$PEPMASS[1] - 18.033826}
        if (sub_ref$ADDUCT[1]=="M-H"){NM = sub_ref$PEPMASS[1] + 1.007276}
        if (sub_ref$ADDUCT[1]=="M+Cl"){NM = sub_ref$PEPMASS[1] - 34.969401}

    # Calculate adduct masses
        for (adduct in adducts){
          for (charge in charges){
          if (adduct == "M+H"){
                  new_pep = c(new_pep, (NM + 1.007276*charge)/charge)
                  new_ions = c(new_ions, "Positive")
                  new_adducts = c(new_adducts, "M+H")
                  new_charges = c(new_charges, charge)}

          if (adduct == "M+Na"){
            new_pep = c(new_pep, (NM + 22.989221*charge)/charge)
            new_ions = c(new_ions, "Positive")
            new_adducts = c(new_adducts, "M+Na")
            new_charges = c(new_charges, charge)}

          if (adduct == "M+K"){
            new_pep = c(new_pep, (NM + 38.963158*charge)/charge)
            new_ions = c(new_ions, "Positive")
            new_adducts = c(new_adducts, "M+K")
            new_charges = c(new_charges, charge)}

          if (adduct == "M+NH4"){
            new_pep = c(new_pep, (NM + 18.033826*charge)/charge)
            new_ions = c(new_ions, "Positive")
            new_adducts = c(new_adducts, "M+NH4")
            new_charges = c(new_charges, charge)}

           if (adduct == "M-H"){
            new_pep = c(new_pep, (NM - 1.007276*charge)/charge)
            new_ions = c(new_ions, "Negative")
            new_adducts = c(new_adducts, "M+K")
            new_charges = c(new_charges, charge)}

           if (adduct == "M+Cl"){
            new_pep = c(new_pep, (NM + 34.969401*charge)/charge)
            new_ions = c(new_ions, "Negative")
            new_adducts = c(new_adducts, "M+Cl")
            new_charges = c(new_charges, charge)}
          }
        } # End of calculating adducts

    # Append new metadata for the id:
      nb_adduct = length(new_pep)
      new_sub_ref = c()
      for (i in 1:nb_adduct){new_sub_ref = rbind.data.frame(new_sub_ref, sub_ref[1,])}
      new_sub_ref$PEPMASS = round(new_pep,5)
      new_sub_ref$IONMODE = new_ions
      new_sub_ref$ADDUCT = new_adducts
      new_sub_ref$CHARGE = new_charges   # End of calculating an id

    # Add to the new metadata:
      new_metadata = rbind.data.frame(new_metadata, new_sub_ref)
    }} else {
      new_metadata = data.frame(ref)
  }

  return(new_metadata)
}
