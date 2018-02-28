# This file computes the spatially-explicit whole-leaf energy balance model from TIR data with related meterological data
# from submitted manuscript "Spatiotemporal dynamics of leaf transpiration quantified with time-series thermal imaging"
# code written by Gerald Page and Jean Lienard

require(reshape2)
require(dplyr)
require('RColorBrewer')

system.time({

  # for the fudgeit function (colormap)
  # helper function to put a colorbar on a raster plot
  fudgeit = function(d, colramp, ...){
    require('fields')
    x = seq(0,1,le=length(d))
    y = seq(0,1,le=length(d))
    image.plot(list(x=x,y=y,z=d), col = colramp, legend.only = T, add=T,...)
  }

  colm100=rev((colorRampPalette(c("#A50026",brewer.pal(11, 'Spectral')[2:10],"#313695") ))(100))
  
  nrep=1; doplot=F
  
  SR01Up=0; IR01UpCo=0; IR01DnCo=0; LTemp=0; ATemp=0; Larea=0
 

  
  ## constants--------------------------------------------
  cp = 29.3  # heat capacity of air (J mol-1 C-1)
  eL = 0.98 # leaf emissivity
  c_pw = 4182
  bz = 5.6703 * 10^(-8) # Stefan-Boltzmann Constant
  Absorptivity = 0.70
  Atmospheric_view_factor = 0.8
  Ground_view_factor = 1 - Atmospheric_view_factor
  
  correction_missing_photos = read.table(text='0
                                       4
                                       0
                                       0
                                       3
                                       1
                                       0
                                       0
                                       2
                                       0
                                       0
                                       1
                                       0
                                       0
                                       1
                                       0
                                       1
                                       1
                                       3
                                       1
                                       0
                                       1')
  
  leave_dir = './data/leaves_wdry/'
  my_files <- list.files(leave_dir, pattern = "\\.csv$")
  my_files <- my_files[order(nchar(my_files), my_files)]
  setwd(leave_dir)
  leaves <- lapply(my_files, read.csv)
  
  # compute Rabs:
  setwd('../met/')
  for (i in 1:length(leaves)) {
    m = read.csv(paste0('l', i, '.csv'), header = T)
    # formula: Rabs = incident shortwave*absorptivity + (eL*(atmospheric view factor*atmospheric longwave + ground view factor*ground longwave))
    leaves[[i]]$Rabs[] = (m$SR01Up + SR01Up) * Absorptivity + eL * (Atmospheric_view_factor * (m$IR01UpCo + IR01UpCo) + Ground_view_factor * (m$IR01DnCo + IR01DnCo))
    leaves[[i]]$wind.int = m$wind.int
  }
  
  # read leaf data
  setwd("../leaf_params/")
  larea <- read.csv("single_leaf_areas_R.csv", header=T)
  water <- read.csv("single_leaf_water_R.csv", header=T)
  larea <- larea[order(larea$Leaf),]
  leaf_params.df <- data.frame(larea$Leaf, larea$Area_m + Larea, water$InitialKG, water$FinalKG)
  names(leaf_params.df) <- c("Leaf", "Area_m", "InitKG", "FinalKG")
  leaf_params.df$SWM_A <- leaf_params.df$InitKG/leaf_params.df$Area_m
  leaf_params.df$FWM_A <- leaf_params.df$FinalKG/leaf_params.df$Area_m
  leaf_params.df$mass_loss_g <- (leaf_params.df$InitKG - leaf_params.df$FinalKG) * 1000
  av_lwma <- mean(leaf_params.df$SWM_A) # average starting leaf water mass per area, for dry reference leaf heat capacity estimation
  
  # initialize the bootstrap variables
  med_res = NULL
  min_res = NULL
  max_res = NULL
  
  list_allT = list()
  list_allE = list()
  list_allWL = list()
  
  for (l_i in c(3, 7, 10, 12)) { #computes spatial model only for leaves that remained stationary during image capture
    
    d = leaves[[l_i]]
    leaf_area = leaf_params.df$Area_m[l_i]
    starting_water_mass = leaf_params.df$SWM_A[l_i]
    L_l = sqrt(leaf_params.df$Area_m[l_i]) # leaf length [m]
    
    # initialise new columns
    d$E = NA
    d$Loe = NA
    d$gr = NA
    d$Real_LTemp = d$LTemp + LTemp
    d$Real_DLTemp = d$DLTemp + LTemp
    d$AirT = d$AirT + ATemp
    d$missing = is.na(d$LTemp)
    ran = 2:nrow(d)
    d$gHa_recomputed = NULL
    d$gHa = NULL
    
    # compute radiative conductance term
    d$gr = (4*bz*((d$AirT + 273.15)^3))/cp # radiative conductance
    
    # initialize the bootstrap variables
    results = matrix(NA, nrow=nrep, ncol=1)
    detailed = matrix(NA, nrow=nrep, ncol=nrow(d))
    
    # create empty array and fill with TIR image data
    all_dTemp = array(NA, dim = c(480, 640, nrow(d)))
    setwd("../whole_leaf/")
    for (xp_time in 1:(1+(nrow(d))/15)) {
      all_dTemp[,,(xp_time-1)*15+1] = as.matrix(read.csv(paste0('L', l_i, '_', xp_time, '.csv'), header = F))
    }
    
    # melt array to a data frame
    all_dTemp_df <- melt(all_dTemp)
    all_E <- all_dTemp_df
    
    remove(all_dTemp)
    remove(all_dTemp_df)
    
    system.time({
      all_E_wide <- dcast(all_E, Var3 ~ Var1 * Var2, value.var = 'value')[,-1] # cast to seconds x 307200 data frame
    }) # 32 seconds
    
    # randomly generate a new temperature trace
    if (nrep == 1) {
      # only one trial => no noise
      measurement_error = 0
    } else {
      measurement_error = rnorm(sum(!d$missing), 0, 0.3)
    }
    d$DLTemp = spline(d$Seconds[!d$missing], measurement_error + d$DLTemp[!d$missing], xout = d$Seconds)$y
    
    setwd("../whole_leaf_out_vect/")
    all_E_file = paste0('./leaf_EA/leaf_EA_', l_i, '.rds')
    all_dTemp_file = paste0('./leaf_T/allT_', l_i, '.rds')
    
    
    ## Create function to fit smoothing spline to temperature data 
    SSpline <- function(x, y, n = nrow(all_E_wide), ...) {
      ## fit the spline to x, and y
      mod <- smooth.spline(x, y, ...)
      ## predict from mod for n points over range of x
      pred.dat <- seq(from = min(x), to = max(x), length.out = n)
      ## predict
      preds <- predict(mod, x = pred.dat)
      ## return
      preds
    }
    
    # determine which pixels contain leaf temperature data
    LT_missing <- apply(all_E_wide[!d$missing,], 2, anyNA)
    good_cols <- as.numeric(which(LT_missing == FALSE))
    
    if (file.exists(all_E_file) & file.exists(all_dTemp_file))
    {
      all_E = readRDS(all_E_file)
      all_dTemp = readRDS(all_dTemp_file)
    } else {
      no_temp <- which(is.na(all_E_wide[16,]))
      # fit smoothing spline to all pixels that have temperature data for all images in sequence
      system.time({
        res2 <- apply(all_E_wide[d$Seconds[!d$missing], good_cols], 2,
                      function(x, y, ...) { SSpline(x, y, ...)$y },
                      x = d$Seconds[!d$missing])
      }) # 11 seconds!
      
      # replace raw data with spline data for matching columns in all_E_wide
      all_E_wide[,colnames(res2)] = res2
      
      # create new dataframe of a change in leaf temperature for later E calcs
      res3 <- diff(res2)
      x <- rep(0, ncol(res3))
      res3 <- rbind(x, res3)
      
      #assign starting water mass for each model
      water_mass_dr = array(starting_water_mass, dim = dim(res2)) 
      colnames(water_mass_dr) <- colnames(res2)
      water_mass_simple <- water_mass_dr
      water_mass_schym <- water_mass_dr
      
      # res2 = leaf temperature data every second, by columns
      # res3 = change in leaf temperature
      # water_mass = water mass each second (to be iterated over by row)
      
      Lwatts = bz * (res2 + 273.15)^4
      Loe = (eL*bz*((d$AirT + 273.15)^4)) + (cp*d$gr*(res2 - d$AirT))
      gHa_recomputed <- rep(NA, nrow(d))
      gHa_recomputed[2:nrow(d)] = (d$Rabs[2:nrow(d)] - (c_pw*av_lwma*(d$DLTemp[2:nrow(d)] - d$DLTemp[1:nrow(d)-1])) - ((eL*bz*((d$AirT[2:nrow(d)]+ 273.15)^4)) + (cp*d$gr[2:nrow(d)]*(d$DLTemp[2:nrow(d)] - d$AirT[2:nrow(d)])))) / (cp*(d$DLTemp[2:nrow(d)] - d$AirT[2:nrow(d)]))
      gHa_recomputed[1] = gHa_recomputed[2]
      gHa = 1.4 * (0.135 * sqrt(d$wind.int/(0.72 * L_l)))
      
      E_dr <- array(NA, dim = dim(res2))
      colnames(E_dr) <- colnames(res2)
      E_simple <- E_dr
      E_schym <- E_dr
      
      system.time({
        for (i in 1:length(ran)){
          E_dr[i,] <- d$Rabs[i] - (cp*gHa_recomputed[i] * (res2[i,] - d$AirT[i])) -
            Loe[i,] - (res3[i,] * water_mass_dr[i,]*c_pw)
          water_mass_dr[i+1,] <- water_mass_dr[i,] - (E_dr[i,]/44000 * 18.01528 / 1e3)
        }
      }) # 0.5 seconds!
      
      for (i in 1:length(ran)){
        E_simple[i,] <- d$Rabs[i] - (cp*gHa[i] * (res2[i,] - d$AirT[i])) -
          Loe[i,] - (res3[i,] * water_mass_simple[i,]*c_pw)
        water_mass_simple[i+1,] <- water_mass_simple[i,] - (E_simple[i,]/44000 * 18.01528 / 1e3)
      }
      
      T_b <- rep(NA, nrow(d))
      v_a <- rep(NA, nrow(d))
      N_rel <- rep(NA, nrow(d))
      N_rec <- rep(NA, nrow(d))
      C_l <- rep(NA, nrow(d))
      N_nul <- rep(NA, nrow(d))
      k_a <- rep(NA, nrow(d))
      h_c <- rep(NA, nrow(d))
      H_schym <- array(NA, dim = dim(res2))
      
      for (i in 1:length(ran)){
        N_pr = 0.71 # dimensionless Prandtl number
        T_b[i] = mean((d$AirT[i] + res2[i,] + 273.15*2) / 2) # boundary layer temperature [K]
        v_a[i] = (9*(10^-8) * T_b[i]) - 1.13*(10^-5) # kinematic viscosity of air [m s-1]
        N_rel[i] = d$wind.int[i] * L_l / v_a[i] # leaf average Reynold's number
        N_rec1 = 3000 # cricitical reynolds number, set as per Schymanski et al. 2013
        N_rec[i] = (N_rel[i] + N_rec1 - abs(N_rec1 - N_rel[i]) ) / 2 # substituting this term for Nrec in the calculation of C_l, following Schymanski et al. 2013
        C_l[i] = 0.037 * N_rec[i]^(4/5) - 0.664 * N_rec[i]^(1/2)
        N_nul[i] = (0.037 * N_rel[i]^(4/5) - C_l[i]) * N_pr^(1/3)
        k_a[i] = 6.84 * 10^(-5) * T_b[i] + 5.62 * 10^(-3)
        h_c[i] = k_a[i] * N_nul[i] / L_l
        H_schym[i,] = 2 * h_c[i] * (res2[i,] - d$AirT[i])
        
        E_schym[i,] <- d$Rabs[i] - H_schym[i,] - Loe[i,] - (res3[i,] * water_mass_schym[i,]*c_pw)
        water_mass_schym[i+1,] <- water_mass_schym[i,] - (E_schym[i,]/44000 * 18.01528 / 1e3)
      }
      
      
      leaf_water_dr <- array(NA, dim(all_E_wide))
      colnames(leaf_water_dr) <- colnames(all_E_wide)
      leaf_water_simple <- leaf_water_dr
      leaf_water_schym <- leaf_water_dr
      leaf_E_dr <- leaf_water_dr
      leaf_E_simple <- leaf_water_dr
      leaf_E_schym <- leaf_water_dr
      
      leaf_water_dr[,colnames(water_mass_dr)] <- water_mass_dr
      leaf_water_simple[,colnames(water_mass_simple)] <- water_mass_simple
      leaf_water_schym[,colnames(water_mass_schym)] <- water_mass_schym
      
      leaf_E_dr[,colnames(E_dr)] <- E_dr
      leaf_E_simple[,colnames(E_simple)] <- E_simple
      leaf_E_schym[,colnames(E_schym)] <- E_schym
      
      ## Now, reshape it from 196x307200 to 640x480xseconds (code exists above already)
      system.time({
        leaf_E1_dr <- melt(leaf_E_dr)
        remove(leaf_E_dr); gc()
        vars <- colsplit(leaf_E1_dr$Var2, "_", c("Var1", "Var2"))
        leaf_E1_dr <- data.frame(vars$Var1, vars$Var2, leaf_E1_dr$Var1, leaf_E1_dr$value)
        names(leaf_E1_dr) <- names(all_E)
        leaf_EA_dr <- acast(leaf_E1_dr, Var1 ~ Var2 ~ Var3, value.var = 'value') # convert back to array
        remove(leaf_E1_dr); gc()
        saveRDS(leaf_EA_dr, paste0('leaf_EA_dr_', l_i, '.rds')) # save transpiration array to file
        remove(leaf_EA_dr); gc()
      }) # colsplit takes 131 seconds on the workstation
      
      leaf_E1_simple <- melt(leaf_E_simple)
      remove(leaf_E_simple); gc()
      vars <- colsplit(leaf_E1_simple$Var2, "_", c("Var1", "Var2"))
      leaf_E1_simple <- data.frame(vars$Var1, vars$Var2, leaf_E1_simple$Var1, leaf_E1_simple$value)
      names(leaf_E1_simple) <- names(all_E)
      leaf_EA_simple <- acast(leaf_E1_simple, Var1 ~ Var2 ~ Var3, value.var = 'value') # convert back to array
      # remove(leaf_E1_simple); gc()
      saveRDS(leaf_EA_simple, paste0('leaf_EA_simple_', l_i, '.rds')) # save transpiration array to file
      remove(leaf_EA_simple); gc()
      
      leaf_E1_schym <- melt(leaf_E_schym)
      vars <- colsplit(leaf_E1_schym$Var2, "_", c("Var1", "Var2"))
      leaf_E1_schym <- data.frame(vars$Var1, vars$Var2, leaf_E1_schym$Var1, leaf_E1_schym$value)
      names(leaf_E1_schym) <- names(all_E)
      remove(leaf_E_schym); gc()
      leaf_EA_schym <- acast(leaf_E1_schym, Var1 ~ Var2 ~ Var3, value.var = 'value') # convert back to array
      remove(leaf_E1_schym); gc()
      saveRDS(leaf_EA_schym, paste0('leaf_EA_schym_', l_i, '.rds')) # save transpiration array to file
      remove(leaf_EA_schym); gc()
      
      leaf_water1_dr <- melt(leaf_water_dr)
      leaf_water1_dr <- data.frame(vars$Var1, vars$Var2, leaf_water1_dr$Var1, leaf_water1_dr$value)
      names(leaf_water1_dr) <- names(all_E)
      remove(leaf_water_dr); gc()
      leaf_waterA_dr <- acast(leaf_water1_dr, Var1 ~ Var2 ~ Var3, value.var = 'value') # convert back to array
      remove(leaf_water1_dr); gc()
      saveRDS(leaf_waterA_dr, paste0('leaf_waterA_dr_', l_i, '.rds')) # save leaf water array to file
      remove(leaf_waterA_dr); gc()
      
      leaf_water1_simple <- melt(leaf_water_simple)
      leaf_water1_simple <- data.frame(vars$Var1, vars$Var2, leaf_water1_simple$Var1, leaf_water1_simple$value)
      names(leaf_water1_simple) <- names(all_E)
      remove(leaf_water_simple); gc()
      leaf_waterA_simple <- acast(leaf_water1_simple, Var1 ~ Var2 ~ Var3, value.var = 'value') # convert back to array
      remove(leaf_water1_simple); gc()
      saveRDS(leaf_waterA_simple, paste0('leaf_waterA_simple_', l_i, '.rds')) # save leaf water array to file
      remove(leaf_waterA_simple); gc
      
      leaf_water1_schym <- melt(leaf_water_schym)
      leaf_water1_schym <- data.frame(vars$Var1, vars$Var2, leaf_water1_schym$Var1, leaf_water1_schym$value)
      names(leaf_water1_schym) <- names(all_E)
      remove(leaf_water_schym); gc()
      leaf_waterA_schym <- acast(leaf_water1_schym, Var1 ~ Var2 ~ Var3, value.var = 'value') # convert back to array
      remove(leaf_water1_schym); gc()
      saveRDS(leaf_waterA_schym, paste0('leaf_waterA_schym_', l_i, '.rds')) # save leaf water array to file
      remove(leaf_waterA_schym); gc()
      
      all_temp <- melt(all_E_wide)
      all_temp <- bind_cols(vars, all_temp)
      all_temp$variable <- rep(1:nrow(d), times = 307200)
      names(all_temp) <- names(all_E)
      remove(all_E_wide); gc()
      leaf_temp <- acast(all_temp, Var1 ~ Var2 ~ Var3, value.var = 'value') # convert back to array
      remove(all_temp); gc()
      saveRDS(leaf_temp, paste0('allT_', l_i, '.rds')) # save interpolated temperature to file.
      remove(leaf_temp); gc()
        
      print(paste('Finished Leaf', l_i))
      
    }
    # 
    # # plot spatial data
    # x_min = min(apply(leaf_EA_schym[,,2], 1, function(...)min(which(!is.na(c(...))), na.rm=T), na.rm=T))/640
    # x_max = max(apply(leaf_EA_schym[,,2], 1, function(...)max(which(!is.na(c(...))), na.rm=T), na.rm=T))/640
    # y_min = min(apply(leaf_EA_schym[,,2], 2, function(...)min(which(!is.na(c(...))), na.rm=T), na.rm=T))/640
    # y_max = max(apply(leaf_EA_schym[,,2], 2, function(...)max(which(!is.na(c(...))), na.rm=T), na.rm=T))/640
    # 
    # 
    # system.time({
    #   summed_E = apply(leaf_EA_schym, c(1,2), sum, na.rm=T)
    #   summed_E[summed_E == 0] = NA
    #   plot(1,1,type='n', xlim=c(x_min, x_max), ylim=c(y_min, y_max), xaxt='n', yaxt='n', bty='n', xlab='', ylab='')
    #   image(summed_E, col=colm100, add=T)
    #   mtext(side=3, 'Total water loss (mg)', font=2, cex=1.25)
    #   fudgeit(range(c(summed_E), na.rm=T) / 44000 * 18.01528 * leaf_area / sum(!is.na(summed_E[])) * 1e3, colm100)
    #   # dev.off()
    # })
    # 
  }

})

