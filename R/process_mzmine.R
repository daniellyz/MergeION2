#' Create metadata from mzmine feature quantification table 
#'
#' Function used to build molecular networks
#'
#' @importFrom stringr str_detect
#' @export
#
process_mzmine<-function(mzmine_file, polarity = c("Positive", "Negative"), 
                combine.adduct = FALSE, remove.halogene = FALSE, remove.minor.coelution = FALSE, mz_search = 0.01, rt_search = 0.2){
  
  #if (polarity == "Positive" & combine.adduct){ref_adducts = c("M+H","M+2H","M+Na","M+K","M+NH4", "M+")}
  #if (polarity == "Positive" & !add.adduct){ref_adducts = "M+H"}
  #if (polarity == "Negative" & add.adduct){ref_adducts = c("M-H","M+Cl", "M+HCOO-","M+CH3COO-")}
  #if (polarity == "Negative" & !add.adduct){ref_adducts = "M-H"}
  
  ref = read.csv(mzmine_file,sep=";",dec=".",header=T)
  if (ncol(ref)==1){ref = read.csv(mzmine_file,sep=",",dec=".",header=T)}  
  if (ncol(ref)==1){ref = read.csv(mzmine_file,sep="\t",dec=".",header=T)}   
  
  ##############
  ## Formating##
  ##############
  
  ind1 = which(colnames(ref) == "row.ID")
  ind2 = which(colnames(ref) == "row.m.z")
  ind3 = which(colnames(ref) == "row.retention.time")
  idx = which(colnames(ref) == "X")
  inds = which(str_detect(colnames(ref), "Peak")) # For feature quantification
  
  if (length(idx)==1){ref = ref[,-idx]}
  colnames(ref)[ind1] = "ID"
  colnames(ref)[ind2] = "PEPMASS"
  colnames(ref)[ind3] = "RT"
  colnames(ref)[inds] = str_remove(colnames(ref)[inds], ".Peak.area")
  
  ref = remove_duplicate_feature(ref,mz_search,rt_search) # MZMine sometimes give duplicate features
  
  avg = apply(ref[,inds,drop=FALSE], 1,  mean) # average 
  dat = cbind(ref[,c(ind1, ind2, ind3)], AVG = avg) # Feature table
  ft = cbind(ID = ref[,ind1], ref[,inds, drop=FALSE]) # Quantification table
  
  #########################
  ## Deisotoping halogene##
  #########################
  
  if (nrow(dat)>0){
    
    new_dat = c()
    new_ft = c()
    
    ft = ft[order(dat$RT),]
    dat = dat[order(dat$RT),]
    
    rt_cluster = cut_rt_list_mzmine(dat$RT, rt_search*2)
    rt_id = unique(rt_cluster)
  
    for (rt in rt_id){  
      valid = which(rt_cluster == rt)
      tmp_dat = dat[valid,,drop=FALSE]
      tmp_ft = ft[valid,,drop=FALSE]
      
      tmp_ft = tmp_ft[order(tmp_dat$PEPMASS),]
      tmp_dat = tmp_dat[order(tmp_dat$PEPMASS),]

      if (nrow(tmp_dat)>1){
        tmp_sp = cbind.data.frame(tmp_dat$PEPMASS, tmp_dat$AVG, tmp_dat$ID)
        output = simple_deisotope(tmp_sp, remove.halogene)
        filter = match(output[,3], tmp_dat$ID)
        tmp_dat = tmp_dat[filter,,drop=FALSE]
        tmp_ft = tmp_ft[filter,,drop=FALSE]
      }
      new_ft = rbind.data.frame(new_ft, tmp_ft)
      new_dat = rbind.data.frame(new_dat,tmp_dat)
    }
  } else {
    new_dat = dat
    new_ft = ft
  }
  
  #####################
  ## Grouping adducts##
  #####################
  
  if (combine.adduct){
    
    dat= new_dat
    ft = new_ft
    
    # Label & group adduct:
    dat_adducts = c()
    if (nrow(dat)>0){
      dat = dat[order(dat$RT),]
      rt_cluster = cut_rt_list_mzmine(dat$RT, rt_search*2)
      rt_id = unique(rt_cluster)

      for (rt in rt_id){  
        valid = which(rt_cluster == rt)
        tmp_dat = dat[valid,,drop=FALSE]
        tmp_dat = tmp_dat[order(tmp_dat$PEPMASS),]
        tmp_dat2 = group_adducts(tmp_dat, polarity, mz_search)
        dat_adducts = rbind(dat_adducts, tmp_dat2)
      }
    }
    
    dat_adducts = dat_adducts[order(dat_adducts$ID),]
    dat = dat[order(dat$ID),]
    ft = ft[order(ft$ID),]
    
    # Combine adduct:
    
    new_dat = c()
    new_ft = c()
    
    ID2list = unique(dat_adducts$ID2)
    NI = length(ID2list)
    ID_LIST = rep("0", NI)
    PEPMASS_LIST = rep(0, NI)
    RT_LIST = rep(0, NI)
    AVG_LIST = rep(0, NI)
    ADDUCT_TYPE = rep("0", NI)

    for (i in 1:NI){
      
      valid = which(dat_adducts$ID2 == ID2list[i])
      tmp_dat3 = dat_adducts[valid,,drop=FALSE]
      tmp_ft3 = ft[valid,,drop=FALSE]
      
      ID_LIST[i] = ID2list[i]
      PEPMASS_LIST[i] = round(min(dat_adducts$PEPMASS[valid]),4)
      RT_LIST[i] =  round(mean(dat_adducts$RT[valid]), 2)
      AVG_LIST[i] = sum(dat_adducts$AVG[valid])
      ADDUCT_TYPE[i] = paste0(unique(dat_adducts$ADDUCT[valid]), collapse = ":")
      tmp_ft = colSums(tmp_ft3)
      tmp_ft[1] =  ID2list[i]
      new_ft = rbind.data.frame(new_ft, tmp_ft)
    }
  
    colnames(new_ft) = colnames(ft)
    new_dat = cbind.data.frame(ID = ID_LIST, PEPMASS = PEPMASS_LIST, RT = RT_LIST, AVG = AVG_LIST, ADDUCT_TYPE = ADDUCT_TYPE)
  } else {
    if (polarity == "Positive"){dat$ADDUCT_TYPE ="M+H"}
    if (polarity == "Negative"){dat$ADDUCT_TYPE ="M-H"}
    new_dat = dat
    new_ft = ft
  }
  
  #############################
  ## Minor co-eluted species ##
  #############################
  
  if (remove.minor.coelution){
    
      dat= new_dat
      ft = new_ft
      
      if (nrow(dat)>0){
        
          new_dat = c()
          new_ft = c()
          
          ft = ft[order(dat$RT),]
          dat = dat[order(dat$RT),]
          
          rt_cluster = cut_rt_list_mzmine(dat$RT, rt_search*2)
          rt_id = unique(rt_cluster)
   
          for (rt in rt_id){  
            valid = which(rt_cluster == rt)
            tmp_dat = dat[valid,,drop=FALSE]
            tmp_ft = ft[valid,,drop=FALSE]
            tmp_adduct= lapply(tmp_dat$ADDUCT_TYPE, function(x) strsplit(x, ":")[[1]])
            
            valid_ind = c(which.max(tmp_dat$PEPMASS), which.max(tmp_dat$AVG), which(sapply(tmp_adduct, length)>1))
            valid_ind = unique(valid_ind)
            
            tmp_dat2 = tmp_dat[valid_ind,,drop=FALSE]
            tmp_ft2 = tmp_ft[valid_ind,,drop=FALSE]
            
            new_dat = rbind.data.frame(new_dat, tmp_dat2)
            new_ft = rbind.data.frame(new_ft, tmp_ft2)
          }
      }
  }
  
  output_metadata = cbind.data.frame(new_dat[,c(1:3,5),drop=FALSE], new_ft[,-1,drop=FALSE])
  return(output_metadata)
}

#############################
######Internal function######
#############################

cut_mz_list<-function(mzlist, mz_window){
  
  N=length(mzlist)
  
  f=1
  mz_feature=c(0, N) 
  t0 = 1 # Start index of a cluster
  
  for (k in 2:N){
    min_mz = min(mzlist[t0:(k-1)])
    avg_mz = mean(mzlist[t0:(k-1)])
    
    if (mzlist[k] - min_mz > mz_window & mzlist[k] - avg_mz > mz_window/2){
      mz_feature[t0:(k-1)] = f 
      f = f + 1
      t0 = k
    }
  }
  mz_feature[t0:N] = f
  return(mz_feature)
}

cut_rt_list_mzmine<-function(rtlist, rt_window){
  
  N=length(rtlist)
  
  f=1
  rt_feature=c(0, N)
  t0 = 1 # Start index of a cluster
  
  for (k in 2:N){
    max_rt = max(rtlist[t0:(k-1)])
    min_rt = min(rtlist[t0:(k-1)])
    avg_rt = median(rtlist[t0:(k-1)])
    if (rtlist[k] - min_rt > rt_window & rtlist[k] - avg_rt > rt_window/2){
      rt_feature[t0:(k-1)] = f
      f = f + 1
      t0 = k
    }
  }
  rt_feature[t0:N] = f
  
  return(rt_feature)
}

remove_duplicate_feature<-function(ref, mz_search, rt_search){
  
  ref0 =ref
  
  ref = ref[order(ref$PEPMASS),]
  ref$mlabel = cut_mz_list(ref$PEPMASS, mz_search*2)
  ref = ref[order(ref$RT),]
  ref$rlabel = cut_rt_list_mzmine(ref$RT, rt_search*2)
    
  ref$flabel = paste(ref$mlabel, ref$rlabel, sep = "-")
  flabels = unique(ref$flabel)

  new_ref = c()
  
  for (i in 1:length(flabels)){
    valid = which(ref$flabel == flabels[i])
    tmp_ref = ref[valid,,drop=FALSE]
    
    tmp_id = tmp_ref$ID[1]
    tmp_ft = colMeans(tmp_ref[,2:3,drop=FALSE])
    tmp_profile = colSums(tmp_ref[,4:(ncol(tmp_ref)-3),drop=FALSE])
    
    tmp_ref = c(tmp_id, tmp_ft, tmp_profile)

    new_ref = rbind.data.frame(new_ref, tmp_ref)
  }
  
  colnames(new_ref) = colnames(ref0)
  rownames(new_ref) = NULL
  
  return(new_ref)
}

simple_deisotope<-function(sp, halogene){
  
  # Monoisotope:
  
  NP = nrow(sp) # Nb of peaks
  ref_dis = 1
  charge_state = rep(0, NP)
  
  for (i in 2:NP){
    
    delta_mz = sp[i,1] - sp[1:(i-1), 1]
    errors = abs(delta_mz - ref_dis)
    valid = which(errors<0.1)
    
    if (length(valid)>0){
      delta_int = sp[i,2]/sum(sp[valid, 2])
      if(delta_int<0.5){
        charge_state[i] = 100
        #charge_state[valid] = 1   
     }
    }
  }
  
  valid = which(charge_state!=100)
  sp = sp[valid,,drop=FALSE]
  charge_state = charge_state[valid]
  
  # Chloride:
  
  if (halogene){
      NP = nrow(sp) # Nb of peaks
      ref_dis = 2

      for (i in 2:NP){
    
          delta_mz = sp[i,1] - sp[1:(i-1), 1]
          errors = abs(delta_mz - ref_dis)
          valid = which(errors<0.1)
    
          if (length(valid)>0){
            delta_int = sp[i,2]/sum(sp[valid, 2])
            if (!is.null(delta_int)){
            if(delta_int<1.2){
              charge_state[i] = 200
            }}
          }
      }
      
      valid = which(charge_state!=200)
      sp = sp[valid,,drop=FALSE]
  }
  
  return(sp)
}

group_adducts<-function(sp, polarity, mz_error){
  
 # tmp_sp = cbind.data.frame(tmp_dat$PEPMASS, tmp_dat$ID)
  
  NS = nrow(sp)
  adduct_label = rep("0", NS)
  id_label = sp$ID
  
  for (i in 1:NS){
    
    if (adduct_label[i]=="0"){
      
        mz = sp$PEPMASS[i]
        mz_dist = sp$PEPMASS- mz
    
        if (polarity == "Positive"){
        
          ind0 = which(abs(mz_dist - 0)<=mz_error)
          ind1 = which(abs(mz_dist - 17.02655)<=mz_error)
          ind2 = which(abs(mz_dist - 21.98194)<=mz_error)
          ind3 = which(abs(mz_dist - 37.95588)<=mz_error)
          
          if (length(ind1)==1 || length(ind2)==1 || length(ind3)==1){
            adduct_label[ind0] = "M+H"
            adduct_label[ind1] = "M+NH4"
            adduct_label[ind2] = "M+Na"
            adduct_label[ind3] = "M+K"
            id_label[c(i,ind1, ind2, ind3)] = sp$ID[i]
          }
        }
        
        if (polarity == "Negative"){
          
          ind0 = which(abs(mz_dist - 0)<=mz_error)
          ind1 = which(abs(mz_dist - 35.97668)<=mz_error)
          ind2 = which(abs(mz_dist - 44.998203)<=mz_error)
          ind3 = which(abs(mz_dist - 60.02113)<=mz_error)
          
          if (length(ind1)==1 || length(ind2)==1 || length(ind3)==1){
            adduct_label[ind0] = "M-H"
            adduct_label[ind1] = "M+Cl"
            adduct_label[ind2] = "M+COO"
            adduct_label[ind3] = "M+CH3COO"
            id_label[c(i,ind1, ind2, ind3)] = sp$ID[i]
          }
        }
    }
  }
  
  if (polarity == "Positive"){
    adduct_label[adduct_label=="0"] = "M+H"
  }
  if (polarity == "Negative"){
    adduct_label[adduct_label=="0"] = "M-H"
  }
  
  sp2 = cbind.data.frame(sp, ADDUCT = adduct_label, ID2 = id_label)
  return(sp2)
}

  