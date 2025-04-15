figures_path <- "C:\\Users\\itayta.WISMAIN\\Desktop\\allcells\\figures_final_layout"







supp_figure_jsd <- function(paths) 
{
  write_path <- sprintf("%s\\supp_figure_jsd\\",figures_path)
  dir.create(write_path)
  
  true_data_paths <- sprintf("%s\\%s", paths, "equalized")
  
  sessions_to_use=c(5:8,13:16)
  
  unique_paths_both <- 
    unique(unlist(lapply(str_split(true_data_paths, "Matrices"), function(sr) {return(sr[[1]])})))
  
  percent_df_list <- lapply(unique_paths_both, 
                            function(p) {
                              load(sprintf("%s\\percent_df.R", p))
                              return(final_df_list)
                            })
  
  all_df <- c()
  
  for(i in 1:len(percent_df_list)) {
    measured_percent_df <- percent_df_list[[i]] 
    measurements_df <- c()
    
    for (col_idx in 1:ncol(measured_percent_df[[1]])) {
      
      
      col_sessions <- 
        do.call(rbind,
                lapply(measured_percent_df,
                       function(df_rep) {
                         return(df_rep[sessions_to_use, col_idx])
                       }))
      
      measurements_df <- cbind(measurements_df,
                               colMeans(col_sessions))
    }
    
    measurements_df <- cbind(measurements_df,
                             sessions_to_use)
    
    measurements_df <- as.data.frame(measurements_df)
    
    measurements_df <- cbind(measurements_df,
                             rep(str_split(str_split(unique_paths_both[[i]], "CA1\\\\")[[1]][2], "\\\\")[[1]][1],
                                 times=len(sessions_to_use)))
    
    all_df <- rbind(all_df,
                    measurements_df)
  }
  
  colnames(all_df) <- c(colnames(measured_percent_df[[1]]),
                        "Session", "Mice")
  
  estimation_list <- list()
  likelihood_list <- list()
  measured_percent_list <- list()
  estimated_per_percent_list <- list()
  session_metavar_list <- list()
  sessions_to_use <- c(5:8,13:16)
  
  for (working_path in true_data_paths) {
    
    # extract all paths
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
    
    # Run through all paths
    for (tpath in session_paths) {
      
      idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
      if (sum(idx == sessions_to_use) == 0){
        next
      }
      
      if (len((grep("pct_estimate_df.R", list.files(sprintf("%s\\simulations_JSD_likelihood", tpath))))) == 0) {
        print(tpath)
        print("Missing")
        next
      }
      
      
      
      final_df <- get_fit_params_figures_from_path(tpath, just_likelihood = T, use_JSD = T, edge_free = F)
      
      likelihood_idx <- which(colnames(final_df) == "likelihood")
      likelihood <- final_df[,likelihood_idx] # Put last column on different vector
      final_df <- final_df[,-likelihood_idx] # remove last column
      
      # Mistakely, dataframes were saved reversed with the JSD estimation, meaning that instea of 0.5, 0.6 ... .1 
      # They were saved as 1, 0.9 ... 0.5
      likelihood_t <- likelihood
      likelihood <- rev(likelihood)
      names(likelihood) <- names(likelihood_t)
      means_vec <- rowMeans(final_df)
      means_vec_t <- means_vec
      means_vec <- rev(means_vec)
      names(means_vec) <- names(means_vec_t) 
      
      if ("measured_pct" %in% colnames(final_df)) {
        col_idx <- which(colnames(final_df) == "measured_pct")
        measured_p <- final_df[1,col_idx]
        final_df <- final_df[,-col_idx]
      } else {
        
        #### IN older versions, after likelihood estimation, mistakenly the dataframe was saved without
        # the measured percentage, however as the likelihood, reflects the probability to receive MEASURED PERCENTAGE
        # under the distribution of simulations, via integration, we can re-integrate them and find out EXACTLY numerically
        # what was the measured percentage of that session (as it is random - equalizing the prior is a random process)
        # however, this is not an approximation, it's the exact percentage, just a numerical trick to get it
      
        percents_measured <- c()
        
        for (i in 1:nrow(final_df)) {
          
          kdf_func <- approxfun(density(final_df[i,], width=0.10, from=0, to=1))
          measured_p1 <- sapply(seq(0,0.999, by=0.001), function(v) {euclidean(integrate(kdf_func, v, 1, stop.on.error = F)$value, likelihood[i])})
          measured_p2 <- sapply(seq(0.001,1,  by=0.001), function(v) {euclidean(integrate(kdf_func, 0, v, stop.on.error = F)$value, likelihood[i])})
          
          measured_p1 <- seq(0,0.999, by=0.001)[which.min(measured_p1)]
          measured_p2 <- seq(0.001,1,  by=0.001)[which.min(measured_p2)]
          
          percents_measured <- rbind(percents_measured, c(measured_p1, measured_p2))
        }
        
        percents_tab <- table(percents_measured)
        percents_tab <- percents_tab[as.numeric(names(percents_tab)) > 0.01]
        measured_p <- as.numeric(names(which.max(percents_tab)))
      
      }
      
      measured_percent_list <- append(measured_percent_list, list(measured_p))
      estimation_list <- append(estimation_list, list(as.numeric(names(which.max(likelihood)))))
      likelihood_list <- append(likelihood_list, list(likelihood))
      estimated_per_percent_list <- append(estimated_per_percent_list, list(means_vec))
      
      
      
      mice_str <- str_split(str_split(tpath, "CA1\\\\")[[1]][2], "\\\\Matrices")[[1]]
      
      session_metavar_list <- append(session_metavar_list, list(c(idx, mice_str[1], ifelse(grepl("Left", mice_str[2]), "Left", "Right"))))
      
    }
  }
  
  likelihood_df <- do.call(rbind, likelihood_list)
  estimated_per_percent_df <- do.call(rbind, estimated_per_percent_list)
  session_df <- do.call(rbind, session_metavar_list)
  
  estimation_vec <- unlist(estimation_list)
  measured_vec <- unlist(measured_percent_list)
  
  estimated_df <- data.frame(percent_estimated=estimation_vec)
  estimated_df <- cbind(estimated_df, session_df)
  colnames(estimated_df) <- c("Estimated", "Session", "Mice", "Direction")
  
  # Create weighted estimations of % of place cells for both running directions
  both_dir_df <- ddply(estimated_df, .(Mice), 
                       function(sub_df) {
                         left <- sub_df[which(sub_df[,"Direction"] == "Left"),]
                         right <- sub_df[which(sub_df[,"Direction"] == "Right"),]
                         
                         intersecting_sessions <- as.numeric(intersect(left[,"Session"], right[,"Session"]))
                         
                         mice = sub_df[1,"Mice"]
                         
                         for (ses in intersecting_sessions) {
                           left_ses <- left[which(as.numeric(left[,"Session"]) == ses),]
                           right_ses <- right[which(as.numeric(right[,"Session"]) == ses),]
                           
                           # Get active percentage for weighted averages
                           tmp_all_df_row <- all_df[which(all_df[,"Session"] == ses & all_df[,"Mice"] == mice),]
                           
                           right_p <- tmp_all_df_row[,"Right active"]
                           left_p <- tmp_all_df_row[,"Left active"]
                           
                           
                           weighted_estimate_both_dir <- 
                             (right_ses[,"Estimated"] * right_p + left_ses[,"Estimated"] * left_p) / (left_p + right_p)
                           
                           sub_df <- rbind(sub_df,
                                           c(weighted_estimate_both_dir, ses, mice, "Both"))
                         }
                         
                         return(sub_df)
                       })
  
  estimated_per_percent_df <- cbind(estimated_per_percent_df, measured_vec)
  estimated_per_percent_df <- as.data.frame(estimated_per_percent_df)
  estimated_per_percent_df <- cbind(estimated_per_percent_df, session_df)
  
  
  
  melted_likelihood <- melt(t(apply(likelihood_df[,-1], 1, function(r) {return(r/sum(r))})))
  colnames(melted_likelihood) <- c("#", "Simulated percent", "Normalized likelihood")
  
  both_dir_df <- as.data.frame(both_dir_df)
  both_dir_df[,"Estimated"] <- as.numeric(both_dir_df[,"Estimated"])
  mean_sd_df <- ddply(both_dir_df, .(Direction), 
                      function(sub_df) {
                        return(c(mean(sub_df[,"Estimated"]),
                                 sd(sub_df[,"Estimated"])))
                      })
  
  colnames(mean_sd_df) <- c("Direction", "Estimated", "Sd")
  gest <- 
    ggplot(mean_sd_df, aes(x=Direction, y=Estimated)) + 
    geom_boxplot(data=both_dir_df, aes(group=Direction), width=0.5, fill=adjustcolor("gray80", alpha=0.8), color="gray65", size=1) + 
    geom_jitter(data=both_dir_df,
                aes(x=Direction, y=Estimated, group=Direction), position=position_jitter(0.2), color="gray20", size=0.75, alpha=0.1) + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank()) + 
    xlab("Running direction") +
    ylab("Fraction of place cells (estimated)") + 
    ylim(0,1.05)
  
  
  likelihood_mean_sd_df <- ddply(melted_likelihood, .(`Simulated percent`), 
                                 function(sub_df) {
                                   return(c(mean(sub_df[,"Normalized likelihood"]),
                                            sd(sub_df[,"Normalized likelihood"])))
                                 })
  
  colnames(likelihood_mean_sd_df) <- c("Simulated percent", "Mean", "Sd")
  
  melted_likelihood$`Simulated percent` <- factor(melted_likelihood$`Simulated percent` * 100,
                                                  levels=as.character(likelihood_mean_sd_df$`Simulated percent`  * 100))
  likelihood_mean_sd_df$`Simulated percent` <- factor(likelihood_mean_sd_df$`Simulated percent` * 100, 
                                                      levels=as.character(likelihood_mean_sd_df$`Simulated percent` * 100))
  
  glikelihood <- 
    ggplot(likelihood_mean_sd_df, aes(x=`Simulated percent`, y=Mean)) + 
    geom_bar(stat="summary", aes(group=`Simulated percent`), 
             width=0.5, fill="gray50", color="black", position = "dodge") + 
    geom_errorbar(aes(ymin=Mean, ymax=Mean + Sd, group=`Simulated percent`),
                  size=2, width=0.3) +
    # geom_jitter(data=melted_likelihood,
    #             aes(x=`Simulated percent`, y=`Normalized likelihood`, group=`Simulated percent`), 
    #             position=position_jitter(0.05), color="gray40", size=3, alpha=0.5) + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank()) + 
    xlab("Simulated percent (%)") +
    ylab("Normalized likelihood") + 
    ylim(0,1.05)
  
  
  scatter_list <- list()
  for (est_p in as.character(seq(0.5, 1, by=0.1))) {
    
    tmp_df <- data.frame(Measured=estimated_per_percent_df[,"measured_vec"],
                         Estimated=estimated_per_percent_df[,est_p],
                         Session=as.numeric(estimated_per_percent_df[,ncol(estimated_per_percent_df) - 2]))
    
    tmp_df$Session[tmp_df$Session > 8] <- tmp_df$Session[tmp_df$Session > 8] - 8
    tmp_df$Session <- factor(sprintf("%dth", tmp_df$Session), c("5th", "6th", "7th", "8th"))
    
    linear_line_df <- data.frame(x=c(0.3,1),
                                 y=c(0.3,1))
    gs <- 
      ggplot(linear_line_df) +
      geom_line(aes(x=x,y=y), linetype="longdash", color="#c0c1c3",
                size=1.5, alpha=0.5) +    
      geom_point(data=tmp_df, aes(x=Measured, y=Estimated, color=Session), size=1.5) +
      theme_light() +     
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position="NA",
            plot.margin=unit(c(0,0,0,0), "cm")) +
      scale_color_brewer(palette="RdYlBu") + 
      ylim(c(0.3,1)) + 
      xlim(c(0.3,1)) + 
      ylab("") +
      xlab("")
    print(sprintf("MSE = %f",
                  mean((tmp_df$Estimated - tmp_df$Measured) ** 2)))
    
    scatter_list <- append(scatter_list, list(gs))
  }
  
  scatter_list$nrow = 2
  #scatter_list$left="Estimated place cell percentage (%)"
  #scatter_list$bottom="Measured place cell percentage (%)"
  scattersp <- do.call(grid.arrange, scatter_list)
  plot_h <- 1.9 * 4 / 3
  

  aligned <- align_plots(first_row, second_row, align="v")
  
  JSD_vio <- figure_2_pvalue_violin(paths, verbose=T, use_JSD=T)
  
  
  first_row <- grid.arrange(gest, JSD_vio[[1]] + theme(axis.text.x = element_text(angle = 0, vjust=1), legend.position = "NA"), 
                             nrow=1,widths=c(1, 1.5))
  second_row <- grid.arrange(glikelihood, align_plots(scattersp)[[1]], 
                             nrow=1,widths=c(1, 1.5))
  
  ap <- align_plots(first_row, second_row)
  
  
  pdf(file=sprintf("%s\\first_row.pdf", write_path),
      height=plot_h,
      width=plot_h * 2.5)
  plot(ap[[1]])
  dev.off() 
  
  pdf(file=sprintf("%s\\second_row.pdf", write_path),
      height=plot_h,
      width=plot_h * 2.5)
  
  plot(ap[[2]])
  dev.off() 
}

supp_figure_no_edges <- function(paths) 
{
  write_path <- sprintf("%s\\supp_figure_no_edges\\",figures_path)
  dir.create(write_path)
  
  true_data_paths <- sprintf("%s\\%s", paths, "equalized")
  
  sessions_to_use=c(5:8,13:16)
  
  unique_paths_both <- 
    unique(unlist(lapply(str_split(true_data_paths, "Matrices"), function(sr) {return(sr[[1]])})))
  
  percent_df_list <- lapply(unique_paths_both, 
                            function(p) {
                              load(sprintf("%s\\percent_df.R", p))
                              return(final_df_list)
                            })
  
  all_df <- c()
  
  for(i in 1:len(percent_df_list)) {
    measured_percent_df <- percent_df_list[[i]] 
    measurements_df <- c()
    
    for (col_idx in 1:ncol(measured_percent_df[[1]])) {
      
      
      col_sessions <- 
        do.call(rbind,
                lapply(measured_percent_df,
                       function(df_rep) {
                         return(df_rep[sessions_to_use, col_idx])
                       }))
      
      measurements_df <- cbind(measurements_df,
                               colMeans(col_sessions))
    }
    
    measurements_df <- cbind(measurements_df,
                             sessions_to_use)
    
    measurements_df <- as.data.frame(measurements_df)
    
    measurements_df <- cbind(measurements_df,
                             rep(str_split(str_split(unique_paths_both[[i]], "CA1\\\\")[[1]][2], "\\\\")[[1]][1],
                                 times=len(sessions_to_use)))
    
    all_df <- rbind(all_df,
                    measurements_df)
  }
  
  colnames(all_df) <- c(colnames(measured_percent_df[[1]]),
                        "Session", "Mice")
  
  estimation_list <- list()
  likelihood_list <- list()
  measured_percent_list <- list()
  estimated_per_percent_list <- list()
  session_metavar_list <- list()
  sessions_to_use <- c(5:8,13:16)
  
  for (working_path in true_data_paths) {
    
    # extract all paths
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
    
    # Run through all paths
    for (tpath in session_paths) {
      
      idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
      if (sum(idx == sessions_to_use) == 0){
        next
      }
      
      if (len((grep("pct_estimate_df.R", list.files(sprintf("%s\\simulations_JSD_likelihood", tpath))))) == 0) {
        print(tpath)
        print("Missing")
        next
      }
      
      
      
      final_df <- get_fit_params_figures_from_path(tpath, just_likelihood = T, use_JSD = F, edge_free = T)
      
      likelihood_idx <- which(colnames(final_df) == "likelihood")
      likelihood <- final_df[,likelihood_idx] # Put last column on different vector
      final_df <- final_df[,-likelihood_idx] # remove last column
      
      # Mistakely, dataframes were saved reversed with the JSD estimation, meaning that instea of 0.5, 0.6 ... .1 
      # They were saved as 1, 0.9 ... 0.5
      likelihood_t <- likelihood
      likelihood <- rev(likelihood)
      names(likelihood) <- names(likelihood_t)
      means_vec <- rowMeans(final_df)
      means_vec_t <- means_vec
      means_vec <- rev(means_vec)
      names(means_vec) <- names(means_vec_t) 
      
      if ("measured_pct" %in% colnames(final_df)) {
        col_idx <- which(colnames(final_df) == "measured_pct")
        measured_p <- final_df[1,col_idx]
        final_df <- final_df[,-col_idx]
      } else {
        
        #### IN older versions, after likelihood estimation, mistakenly the dataframe was saved without
        # the measured percentage, however as the likelihood, reflects the probability to receive MEASURED PERCENTAGE
        # under the distribution of simulations, via integration, we can re-integrate them and find out EXACTLY numerically
        # what was the measured percentage of that session (as it is random - equalizing the prior is a random process)
        # however, this is not an approximation, it's the exact percentage, just a numerical trick to get it
        
        percents_measured <- c()
        
        for (i in 1:nrow(final_df)) {
          
          kdf_func <- approxfun(density(final_df[i,], width=0.10, from=0, to=1))
          measured_p1 <- sapply(seq(0,0.999, by=0.001), function(v) {euclidean(integrate(kdf_func, v, 1, stop.on.error = F)$value, likelihood[i])})
          measured_p2 <- sapply(seq(0.001,1,  by=0.001), function(v) {euclidean(integrate(kdf_func, 0, v, stop.on.error = F)$value, likelihood[i])})
          
          measured_p1 <- seq(0,0.999, by=0.001)[which.min(measured_p1)]
          measured_p2 <- seq(0.001,1,  by=0.001)[which.min(measured_p2)]
          
          percents_measured <- rbind(percents_measured, c(measured_p1, measured_p2))
        }
        
        percents_tab <- table(percents_measured)
        percents_tab <- percents_tab[as.numeric(names(percents_tab)) > 0.01]
        measured_p <- as.numeric(names(which.max(percents_tab)))
        
      }
      
      measured_percent_list <- append(measured_percent_list, list(measured_p))
      estimation_list <- append(estimation_list, list(as.numeric(names(which.max(likelihood)))))
      likelihood_list <- append(likelihood_list, list(likelihood))
      estimated_per_percent_list <- append(estimated_per_percent_list, list(means_vec))
      
      
      
      mice_str <- str_split(str_split(tpath, "CA1\\\\")[[1]][2], "\\\\Matrices")[[1]]
      
      session_metavar_list <- append(session_metavar_list, list(c(idx, mice_str[1], ifelse(grepl("Left", mice_str[2]), "Left", "Right"))))
      
    }
  }
  
  likelihood_df <- do.call(rbind, likelihood_list)
  estimated_per_percent_df <- do.call(rbind, estimated_per_percent_list)
  session_df <- do.call(rbind, session_metavar_list)
  
  estimation_vec <- unlist(estimation_list)
  measured_vec <- unlist(measured_percent_list)
  
  estimated_df <- data.frame(percent_estimated=estimation_vec)
  estimated_df <- cbind(estimated_df, session_df)
  colnames(estimated_df) <- c("Estimated", "Session", "Mice", "Direction")
  
  # Create weighted estimations of % of place cells for both running directions
  both_dir_df <- ddply(estimated_df, .(Mice), 
                       function(sub_df) {
                         left <- sub_df[which(sub_df[,"Direction"] == "Left"),]
                         right <- sub_df[which(sub_df[,"Direction"] == "Right"),]
                         
                         intersecting_sessions <- as.numeric(intersect(left[,"Session"], right[,"Session"]))
                         
                         mice = sub_df[1,"Mice"]
                         
                         for (ses in intersecting_sessions) {
                           left_ses <- left[which(as.numeric(left[,"Session"]) == ses),]
                           right_ses <- right[which(as.numeric(right[,"Session"]) == ses),]
                           
                           # Get active percentage for weighted averages
                           tmp_all_df_row <- all_df[which(all_df[,"Session"] == ses & all_df[,"Mice"] == mice),]
                           
                           right_p <- tmp_all_df_row[,"Right active"]
                           left_p <- tmp_all_df_row[,"Left active"]
                           
                           
                           weighted_estimate_both_dir <- 
                             (right_ses[,"Estimated"] * right_p + left_ses[,"Estimated"] * left_p) / (left_p + right_p)
                           
                           sub_df <- rbind(sub_df,
                                           c(weighted_estimate_both_dir, ses, mice, "Both"))
                         }
                         
                         return(sub_df)
                       })
  
  estimated_per_percent_df <- cbind(estimated_per_percent_df, measured_vec)
  estimated_per_percent_df <- as.data.frame(estimated_per_percent_df)
  estimated_per_percent_df <- cbind(estimated_per_percent_df, session_df)
  
  
  
  melted_likelihood <- melt(t(apply(likelihood_df, 1, function(r) {return(r/sum(r))})))
  colnames(melted_likelihood) <- c("#", "Simulated percent", "Normalized likelihood")
  
  both_dir_df <- as.data.frame(both_dir_df)
  both_dir_df[,"Estimated"] <- as.numeric(both_dir_df[,"Estimated"])
  mean_sd_df <- ddply(both_dir_df, .(Direction), 
                      function(sub_df) {
                        return(c(mean(sub_df[,"Estimated"]),
                                 sd(sub_df[,"Estimated"])))
                      })
  
  colnames(mean_sd_df) <- c("Direction", "Estimated", "Sd")
  gest <- 
    ggplot(mean_sd_df, aes(x=Direction, y=Estimated)) + 
    geom_boxplot(data=both_dir_df, aes(group=Direction), width=0.5, fill=adjustcolor("gray80", alpha=0.8), color="gray65", size=1) + 
    geom_jitter(data=both_dir_df,
                aes(x=Direction, y=Estimated, group=Direction), position=position_jitter(0.2), color="gray20", size=0.75, alpha=0.1) + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank()) + 
    xlab("Running direction") +
    ylab("Fraction of place cells (estimated)") + 
    ylim(0,1.05)
  
  
  likelihood_mean_sd_df <- ddply(melted_likelihood, .(`Simulated percent`), 
                                 function(sub_df) {
                                   return(c(mean(sub_df[,"Normalized likelihood"]),
                                            sd(sub_df[,"Normalized likelihood"])))
                                 })
  
  colnames(likelihood_mean_sd_df) <- c("Simulated percent", "Mean", "Sd")
  
  melted_likelihood$`Simulated percent` <- factor(melted_likelihood$`Simulated percent` * 100,
                                                  levels=as.character(likelihood_mean_sd_df$`Simulated percent`  * 100))
  likelihood_mean_sd_df$`Simulated percent` <- factor(likelihood_mean_sd_df$`Simulated percent` * 100, 
                                                      levels=as.character(likelihood_mean_sd_df$`Simulated percent` * 100))
  
  glikelihood <- 
    ggplot(likelihood_mean_sd_df, aes(x=`Simulated percent`, y=Mean)) + 
    geom_bar(stat="summary", aes(group=`Simulated percent`), 
             width=0.5, fill="gray50", color="black", position = "dodge") + 
    geom_errorbar(aes(ymin=Mean, ymax=Mean + Sd, group=`Simulated percent`),
                  width=0.5) +
    # geom_jitter(data=melted_likelihood,
    #             aes(x=`Simulated percent`, y=`Normalized likelihood`, group=`Simulated percent`), 
    #             position=position_jitter(0.05), color="gray40", size=3, alpha=0.5) + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank()) + 
    xlab("Simulated percent (%)") +
    ylab("Normalized likelihood") + 
    ylim(0,1.05)
  
  
  scatter_list <- list()
  for (est_p in as.character(seq(0.5, 1, by=0.1))) {
    
    tmp_df <- data.frame(Measured=estimated_per_percent_df[,"measured_vec"],
                         Estimated=estimated_per_percent_df[,est_p],
                         Session=as.numeric(estimated_per_percent_df[,ncol(estimated_per_percent_df) - 2]))
    
    tmp_df$Session[tmp_df$Session > 8] <- tmp_df$Session[tmp_df$Session > 8] - 8
    tmp_df$Session <- factor(sprintf("%dth", tmp_df$Session), c("5th", "6th", "7th", "8th"))
    
    linear_line_df <- data.frame(x=c(0.3,1),
                                 y=c(0.3,1))
    gs <- 
      ggplot(linear_line_df) +
      geom_line(aes(x=x,y=y), linetype="longdash", color="#c0c1c3",
                size=1.5, alpha=0.5) +    
      geom_point(data=tmp_df, aes(x=Measured, y=Estimated, color=Session), size=1.5) +
      theme_light() +     
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position="NA",
            plot.margin=unit(c(0,0,0,0), "cm")) +
      scale_color_brewer(palette="RdYlBu") + 
      ylim(c(0.3,1)) + 
      xlim(c(0.3,1)) + 
      ylab("") +
      xlab("")
    print(sprintf("MSE = %f",
                  mean((tmp_df$Estimated - tmp_df$Measured) ** 2)))
    
    scatter_list <- append(scatter_list, list(gs))
  }
  
  scatter_list$nrow = 2
  #scatter_list$left="Estimated place cell percentage (%)"
  #scatter_list$bottom="Measured place cell percentage (%)"
  scattersp <- do.call(grid.arrange, scatter_list)
  plot_h <- 1.9 * 4 / 3
  
  
  aligned <- align_plots(first_row, second_row, align="v")
  
  JSD_vio <- figure_2_pvalue_violin(paths, verbose=T, use_JSD=F, edge_free = T)
  
  
  first_row <- grid.arrange(gest, JSD_vio[[1]] + theme(axis.text.x = element_text(angle = 0, vjust=1), legend.position = "NA"), 
                            nrow=1,widths=c(1, 1.5))
  second_row <- grid.arrange(glikelihood, align_plots(scattersp)[[1]], 
                             nrow=1,widths=c(1, 1.5))
  
  ap <- align_plots(first_row, second_row)
  
  
  pdf(file=sprintf("%s\\first_row.pdf", write_path),
      height=plot_h,
      width=plot_h * 2.5)
  plot(ap[[1]])
  dev.off() 
  
  pdf(file=sprintf("%s\\second_row.pdf", write_path),
      height=plot_h,
      width=plot_h * 2.5)
  
  plot(ap[[2]])
  dev.off() 
}


supp_figure_place_cells_no_equalized <- function(paths) {
  
  write_path <- sprintf("%s\\supp_figure_no_equalized\\",figures_path)
  dir.create(write_path)
  
  true_data_paths <- sprintf("%s\\%s", paths, "equalized")
  
  sessions_to_use=c(1:16)
  
  unique_paths_both <- 
    unique(unlist(lapply(str_split(true_data_paths, "Matrices"), function(sr) {return(sr[[1]])})))
  
  percent_df_list <- lapply(unique_paths_both, 
                            function(p) {
                              load(sprintf("%s\\percent_df_001.R", p))
                              return(final_df_list)
                            })
  
  all_df <- c()
  
  for(i in 1:len(percent_df_list)) {
    measured_percent_df <- percent_df_list[[i]] 
    measurements_df <- c()
    
    for (col_idx in 1:ncol(measured_percent_df[[1]])) {
      
      
      col_sessions <- 
        do.call(rbind,
                lapply(measured_percent_df,
                       function(df_rep) {
                         return(df_rep[sessions_to_use, col_idx])
                       }))
      
      measurements_df <- cbind(measurements_df,
                               colMeans(col_sessions))
    }
    
    measurements_df <- cbind(measurements_df,
                             sessions_to_use)
    
    measurements_df <- as.data.frame(measurements_df)
    
    measurements_df <- cbind(measurements_df,
                             rep(str_split(str_split(unique_paths_both[[i]], "CA1\\\\")[[1]][2], "\\\\")[[1]][1],
                                 times=len(sessions_to_use)))
    
    all_df <- rbind(all_df,
                    measurements_df)
  }
  
  colnames(all_df) <- c(colnames(measured_percent_df[[1]]),
                        "Session", "Mice")
  
  
  percent_active_df <- all_df[,c(7:9)]
  percent_place_cells_df <- all_df[,c(1:3)]
  percent_all_df <- all_df[,c(4:6)]
  
  colnames(percent_active_df) <- c("Both", "Right", "Left")
  colnames(percent_place_cells_df) <- c("Both", "Right", "Left")
  colnames(percent_all_df) <- c("Both", "Right", "Left")
  
  melted_percent_active <- melt(percent_active_df)
  melted_percent_place_cells <- melt(percent_place_cells_df)
  melted_percent_all <- melt(percent_all_df)
  
  ylabs = c("Percentage of active cells (%)",
            "Percentage of place cells (% of active cells)",
            "Percentage of place cells (% of all cells)")
  
  dfs <- list(melted_percent_active,
              melted_percent_place_cells,
              melted_percent_all)
  
  bar_plots_list <- list()
  for (df_idx in 1:3) {
    df_to_use <- dfs[[df_idx]]
    colnames(df_to_use) <- c("Direction", "Percent")
    sd_df <- ddply(df_to_use, .(Direction), function(sub_df) {return(c(mean(sub_df[,"Percent"]), sd(sub_df[,"Percent"])))})
    colnames(sd_df) <- c("Direction", "Percent", "Sd")
    gr <- 
      ggplot(sd_df, aes(x=Direction, y=Percent)) + 
      geom_jitter(data=df_to_use,
                  aes(x=Direction, y=Percent, group=Direction), position=position_jitter(0.2), color="gray20", size=1, alpha=0.1) + 
      geom_boxplot(data=df_to_use, aes(group=Direction), width=0.5, fill=adjustcolor("gray80", alpha=0.8), color="gray65", size=0.75) + 
      theme_light() +     
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.background = element_blank()) + 
      xlab("Running direction") +
      ylab(ylabs[df_idx]) + 
      ylim(0,1)
    
    bar_plots_list <- append(bar_plots_list,
                             list(gr))
  }
  
  bar_plots_list$nrow <- 1
  first_row <- do.call(grid.arrange, bar_plots_list)
}
 

supp_figure_bayesian_optimizer_minimum <- function() {

    write_path <- sprintf("%s\\supp_figure_bayesian_optimizer\\",figures_path)
    dir.create(write_path)
    
    ## approximate this function
    poly_func <- function(x)   {x ** 5 -(8 * x ** 3) + 10*x + 6}
    X_star <- as.matrix(seq(-2.5,2.5, length.out=2500))
    cov_mat_pred <- rbf_kernel(X_star, X_star, sigma=0.35)
    X_train <- as.matrix(X_star[sample(1:nrow(X_star), 25),])
    
    print(X_train)
    print ("###################")
    
    y <- poly_func(X_train)
    
    
    for (iter in 1:20){
      res <- sample_posterior(as.matrix(X_train),
                              y,
                              as.matrix(X_star),
                              mean_func,
                              sigma_n=0.35,
                              n_sample=20,
                              K_xsxs = cov_mat_pred,
                              plot_kernel = F)
      
      
      
      min_per_x <- apply(res,2,min)
      max_per_x <- apply(res,2,max)
      
      mean_per_pred <- apply(res,2,mean)
      sd_per_pred <- apply(res,2,sd)

      gamma <- (min(y) - mean_per_pred)/sd_per_pred
      prob_to_improve <- pnorm(gamma)
      top_to_improve <- order(prob_to_improve, decreasing=T)
      
      
      expected_improve <- sd_per_pred * (gamma * pnorm(gamma) + dnorm(gamma))
      top_to_improve_exp <- order(expected_improve, decreasing=T)

      
      addition_idx = 1
      left_to_add = 10
      new_X <- c()
      
      while (left_to_add > 0) {
        
        addition_idx <- addition_idx + 1
        
        if (sum(apply(X_train, 1, function(r) {all(r == X_star[top_to_improve[addition_idx],])})) > 0) {
          next
        }

          new_X <- rbind(new_X, X_star[top_to_improve_exp[addition_idx],])
        
        left_to_add <- left_to_add - 1
      }
      
      # Add new points
      X_train <- rbind(X_train, new_X)
      y <- c(y, poly_func(new_X))
      
      print(sprintf("------- Added %d new points (%d Total)", nrow(new_X), nrow(X_train)))
    }
    
    
    approx_df <- data.frame(cbind(apply(res,2, sd), colMeans(res), as.vector(X_star)))
    colnames(approx_df) <- c("SD", "Predicted", "X")
    measured_df <- data.frame(X=as.vector(X_train), Measured=y)
    
    true_df <- data.frame(X=as.vector(X_star), Y=poly_func(as.vector(X_star)))
    
    exampl_plot <- 
    ggplot(approx_df) + 
      geom_ribbon(aes(ymin=Predicted-SD, ymax=Predicted+SD, x=X), col=NA, fill="royalblue4", alpha=0.3) + 
      geom_line(aes(x=X, y=Predicted), col="royalblue4", size=1.25)  + 
      geom_line(data=true_df, aes(x=X, y=Y), alpha=0.5, size=2, linetype="dashed") +
      geom_point(data=measured_df, aes(x=X, y=Measured), col="red") +
      theme_light() +     
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.background = element_blank()) +
      ylab("Y") +
      xlab("X")
  
      load(sprintf("%s\\equalized\\session_7\\KS_simulations_estimate\\param_fit_1.000\\pred_df.R", test_path4))
      predicted <- t(predicted_df)

      values <- list(double_peak_pct=seq(0, 0.5,length.out=7),
                          average_width=seq(0.01, 0.17, by=0.01),
                          noise=seq(0.01, 0.08, by=0.01),
                          sd_width=seq(0.0001, 0.015, length.out=10))
      
       
      predicted <- t(predicted_df[,5:24])
      X_star <- predicted_df[,1:4]
      
      ind_mat <-  combn(1:len(values), 2) 
      
      min_plots <- list()
      min_plots_legend <- list()
      mean_plots <- list()
      
      for (idx in 1:ncol(ind_mat)) {
        ind <- ind_mat[,idx]
        
        v1 <- unlist(values[ind[1]])
        v2 <- unlist(values[ind[2]])
        
        mean_cont_mat <- matrix(rep(0, times=len(v1) * len(v2)), 
                                nrow=len(v1))
        
        min_cont_mat <- matrix(rep(0, times=len(v1) * len(v2)), 
                               nrow=len(v1))
        for (i in 1:len(v1))  {
          for (j in 1:len(v2)) {
            
            mean_cont_mat[i,j] <-  
              mean(predicted[,which(X_star[,ind[1]] == v1[i] & X_star[,ind[2]] == v2[j])])
    
            
            min_cont_mat[i,j] <- 
              min(predicted[,which(X_star[,ind[1]] == v1[i] & X_star[,ind[2]] == v2[j])])
          }
        }
        
        color = colorRampPalette(rev(brewer.pal(n = 7,  name = "Spectral")))
        
        rownames(mean_cont_mat) <- v1#round(v1, digits=3)
        colnames(mean_cont_mat) <- v2#round(v2, digits=3)
        rownames(min_cont_mat) <- v1#round(v1, digits=3)
        colnames(min_cont_mat) <- v2#round(v2, digits=3)
        

        pheatmap(mean_cont_mat, cluster_rows=F, cluster_cols=F, col=color(200))
        pheatmap(min_cont_mat, cluster_rows=F, cluster_cols=F, col=color(200))
        
        melt_min <- melt(min_cont_mat)
        colnames(melt_min) <- c("x", "y", "Cost")
        melt_mean <- melt(mean_cont_mat)
        colnames(melt_mean) <-c("x", "y", "Cost")
        
        thme <-       theme_light() +     
                      theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      legend.position = "NA",
                      plot.margin = unit(c(0,0,0,0), "cm")) 
        
        thme_leg <-  theme_light() +     
                     theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     axis.line = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm")) 
        
        
        # ylab(names(values[ind[2]])) +
        #   xlab(names(values[ind[1]]))

        
        gmin <- 
        ggplot(melt_min, aes(x=x,y=y,fill=Cost)) + 
          geom_tile(color="gray60", size=0.25) + 
          scale_fill_distiller(palette="Spectral") + 
          thme +
          ylab(names(values[ind[2]])) + xlab(names(values[ind[1]]))
        
        ggmin_leg <- 
          ggplot(melt_min, aes(x=x,y=y,fill=Cost)) + 
          geom_tile(color="gray60", size=0.25) + 
          scale_fill_distiller(palette="Spectral") + 
          thme_leg +
          ylab(names(values[ind[2]])) + xlab(names(values[ind[1]]))
        
        gmean <- 
          ggplot(melt_mean, aes(x=x,y=y,fill=Cost)) + 
          geom_tile(color="gray60", size=0.25) + 
          scale_fill_distiller(palette="Spectral") + 
          thme +
          ylab(names(values[ind[2]])) + xlab(names(values[ind[1]]))
          
        
        
        min_plots <- append(min_plots, list(gmin))
        mean_plots <- append(mean_plots, list(gmean))
        min_plots_legend <- append(min_plots_legend, list(gmin_leg))
      }
      
      min_plots$nrow =2
      mean_plots$nrow = 2
      min_plots_legend$nrow = 2
      mean_p <- do.call(grid.arrange, mean_plots)
      min_p <- do.call(grid.arrange, min_plots)
      min_leg_p <- do.call(grid.arrange, min_plots_legend)
      
      
      min_f <- grid.arrange(exampl_plot, min_p, nrow=1, widths=c(1,1.5))
      mean_f <- grid.arrange(exampl_plot, mean_p, nrow=1, widths=c(1,1.5))
      min_leg_f <- grid.arrange(exampl_plot, min_leg_p, nrow=1, widths=c(1,1.5))
      
      pdf(sprintf("%s\\bayesian_optimizer_mean.pdf", write_path),
          height =3,
          width = 7.5)
      plot(mean_f)
      dev.off()
      
      pdf(sprintf("%s\\bayesian_optimizer_min.pdf", write_path),
          height = 3,
          width = 7.5)
      plot(min_f)
      dev.off()
      
      pdf(sprintf("%s\\bayesian_optimizer_min_with_legend.pdf", write_path),
          height = 2.5,
          width = 6.5)
      plot(min_leg_f)
      dev.off()
}
      
      
num_shuffles_df <- function(tp, idx) {
  
  a <- get_spike_train_and_stim_trace_from_path(tp, idx, old = F, equalize_frames = T)
  spike_train <- a[[1]]
  stim_trace <- a[[2]]
  
  fr <- rowMeans(spike_train) / dt
  processed_real <- preprocess_spike_train(spike_train, stim_trace)
  
  true_cells_spike_train <- processed_real$working_cells_spike_train
  true_firing_rate <- processed_real$working_firing_rate
  true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
  
  tmp <- compute_tuning(true_cells_spike_train, stim_trace)
  stim_prob <- tmp[[1]]
  true_tuning_curve <- tmp[[2]]
  rm(tmp)
  
  true_SI <- compute_SI(stim_prob, true_tuning_curve, true_firing_rate)
  
  
  shuffles_df <- c()
  
  for (ij in 2:10) {
    
    res <- c()
    for (n_shuffles in c(500, 1000, 2000, 5000)) {
      place_cell_pval <- 
        sapply(c(1:nrow(true_cells_spike_train)),
               function(i, var2)
               {compute_place_signif(i,
                                     true_cells_spike_train,
                                     stim_trace,
                                     true_firing_rate,
                                     true_SI,
                                     shuffle_type=var2,
                                     num_shuffles = n_shuffles)},
               var2="c")
      # 
      
      place_cell_percentage <- sum(place_cell_pval < 0.05) / len(place_cell_pval)
      print(place_cell_percentage)
      res <- c(res, place_cell_percentage)
      
    }
    
    shuffles_df <- rbind(shuffles_df, res)
    
  }
  
  return(shuffles_df)
}