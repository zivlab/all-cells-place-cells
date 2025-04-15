


get_tunning_corr_df <- function(true_data_paths, 
                                ext="",
                                verbose=F,
                                output_name="",
                                file_name="properties.Rda") {
  
  sessions_to_use <- 1:16
  properties_df <- data.frame()
  correlations_df <- data.frame()
  consecutive_corr_df <- data.frame()
  median_properties_df <- data.frame()
  median_correlations_df <- data.frame()
  median_consecutive_corr_df <- data.frame()  
  
  spearman_correlations_df <- data.frame()
  spearman_consecutive_corr_df <- data.frame()
  spearman_median_correlations_df <- data.frame()
  spearman_median_consecutive_corr_df <- data.frame()  
  
  cosine_df <- data.frame()
  cosine_consecutive_df <- data.frame()
  cosine_median_df <- data.frame()
  cosine_median_consecutive_df <- data.frame()  
  
  list_of_data_frames <- list(mean_corr=correlations_df,
                              mean_consecutive_corr=consecutive_corr_df,
                              median_corr=median_correlations_df,
                              median_consecutive_corr=median_consecutive_corr_df,
                              mean_spearman_corr=correlations_df,
                              mean_spearman_consecutive_corr=consecutive_corr_df,
                              median_spearman_corr=median_correlations_df,
                              median_spearman_consecutive_corr=median_consecutive_corr_df,
                              mean_cosine=cosine_df,
                              mean_cosine_consecutive=cosine_consecutive_df,
                              median_cosine_corr=cosine_median_df,
                              median_cosine_consecutive_corr=cosine_median_consecutive_df)
  
  for (working_path in true_data_paths) {
    
    # extract all paths
    session_paths <- 
      list.dirs(working_path, recursive=F)[grep("session", list.dirs(working_path, recursive=F))]
    
    for (ext in c("", "\\fit_simulation")) { 
      
      # Run through all paths
      for (tpath in session_paths) {
        
        
        idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
        if (sum(idx == sessions_to_use) == 0){
          next
        }
        
        if (len((grep(file_name, list.files(sprintf("%s%s", tpath, ext))))) == 0) {
          print(tpath)
          missing_paths <- c(missing_paths, tpath)
          print("Missing")
          next
        }
        
        if (verbose) {
          print(sprintf("Loading %s", sprintf("%s%s\\%s", tpath, ext, file_name)))
        }
        load(sprintf("%s%s\\%s", tpath, ext, file_name), verbose=T)
        
        
        if (grepl("CA3", working_path)) {
          mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+4) 
        } else {
          mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+3) 
        }
        
        subfield <- ifelse(grepl("CA3", working_path), "CA3", "CA1")
        sim <- ifelse(grepl("fit_simulation", ext), "Simulation", "Real")
        
        for (group in names(result$tuning_properties_by_group)) {
          properties <- result$tuning_properties_by_group[[group]]
          group_property <- apply(properties[,-ncol(properties)], 
                                  2,
                                  mean)
          
          properties_df <- rbind(properties_df,
                                 c(group_property, 
                                   idx, 
                                   ifelse(grepl("Left", tpath), "Left", "Right"), 
                                   mice_str, 
                                   subfield,
                                   group,
                                   sim))
          
          median_group_property <- apply(properties[,-ncol(properties)], 
                                         2,
                                         median)
          
          median_properties_df <- rbind(median_properties_df,
                                        c(median_group_property, 
                                          idx, 
                                          ifelse(grepl("Left", tpath), "Left", "Right"), 
                                          mice_str, 
                                          subfield,
                                          group,
                                          sim))        
        }
        
        # Go through all subsampling schemes
        for (subsamples in names(result$across_subsamples)) {
          
          num_of_subsamples <- as.numeric(substr(subsamples, 
                                                 unlist(gregexpr('[0-9]{1}', subsamples)),
                                                 unlist(gregexpr('[0-9]{1}', subsamples)) + 1))
          
          # Go through all groups
          for (group in names(result$across_subsamples[[subsamples]])) {
            
            if (group == "melted_df") {
              next
            }
            
            # Go through all measured tuning properties
            for (property in names(result$across_subsamples[[subsamples]][[group]]$property_mats)) {
              
              mt <- result$across_subsamples[[subsamples]][[group]]$property_mats[[property]]
              if(!"matrix" %in% class(mt)) {next}
              
              cor_mat <- cor(mt)
              spearman_cor_mat <- cor(mt, method="spearman")
              cosine_sim_mat <- cosine(mt)
              
              
              # REMOVE main diagonal
              for (i in 1:nrow(cor_mat)) {
                cor_mat[i,i] <- NA
                cosine_sim_mat[i,i] <- NA
                spearman_cor_mat[i,i] <- NA
              }
              
              consecutive_cor <-c()
              consecutive_cosine <- c()
              consecutive_spearman <- c()
              
              for(i in 1:(nrow(cor_mat) - 1)) {
                consecutive_cor <- c(consecutive_cor,
                                     cor_mat[i, i+1])
                consecutive_cosine <- c(consecutive_cosine, cosine_sim_mat[i, i+1])
                consecutive_spearman <- c(consecutive_spearman, spearman_cor_mat[i, i+1])                
                
              }
              
              
              for (df_name in names(list_of_data_frames)) { 
                
                if (grepl("cosine", df_name)) {
                  all_comparisions <- cosine_sim_mat
                  consecutive_comparisions <- consecutive_cosine
                } else if (grepl("spearman", df_name)) {
                  all_comparisions <- spearman_cor_mat
                  consecutive_comparisions <- consecutive_spearman
                } else {
                  all_comparisions <- cor_mat
                  consecutive_comparisions <- consecutive_cor
                }
                
                
                if (grepl("consecutive", df_name)) {
                  comparisions_to_use = consecutive_comparisions
                } else {
                  comparisions_to_use = all_comparisions
                }
                
                if (grepl("median", df_name)) {
                  ctf_func <- median
                } else {
                  ctf_func <- mean
                }
                
                similarity_value <- ctf_func(comparisions_to_use, 
                                             na.rm=T)
                
                list_of_data_frames[[df_name]] <- 
                  rbind(list_of_data_frames[[df_name]],
                        c(similarity_value, 
                          idx, 
                          ifelse(grepl("Left", tpath), "Left", "Right"), 
                          mice_str, 
                          subfield,
                          group,
                          sim,
                          property,
                          num_of_subsamples))
                
              }
            }
            
            tuning_corr <- ctf_func(result$across_subsamples[[subsamples]][[group]]$tuning_corr, na.rm=T)
            correlations_df <- rbind(correlations_df,
                                     c(tuning_corr, 
                                       idx, 
                                       ifelse(grepl("Left", tpath), "Left", "Right"), 
                                       mice_str, 
                                       subfield,
                                       group,
                                       sim,
                                       "Tuning corr",
                                       num_of_subsamples))
            
            tuning_corr_median <- median(result$across_subsamples[[subsamples]][[group]]$tuning_corr, na.rm=T)
            
            median_correlations_df <- rbind(median_correlations_df,
                                            c(tuning_corr_median, 
                                              idx, 
                                              ifelse(grepl("Left", tpath), "Left", "Right"), 
                                              mice_str, 
                                              subfield,
                                              group,
                                              sim,
                                              "Tuning corr",
                                              num_of_subsamples))
            
          }
        }
      }
    }
  }
  
  colnames(properties_df) <- c("Peaks",
                               "Time bins",
                               "Firing rate",
                               "SI",
                               "SI_sec",
                               "Spatial bins",
                               "Session",
                               "Direction",
                               "Mice",
                               "Subfield",
                               "Group",
                               "isSim")
  
  for (p in c("Peaks",
              "Time bins",
              "Firing rate",
              "SI",
              "SI_sec",
              "Spatial bins")) {
    properties_df[,p] <- as.numeric(properties_df[,p])
  }
  
  
  
  for (df_name in names(list_of_data_frames)) {
    colnames(list_of_data_frames[[df_name]]) <- 
      c("Corr",
        "Session",
        "Direction",
        "Mice",
        "Subfield",
        "Group",
        "isSim",
        "Property",
        "Subsamples")
    
    list_of_data_frames[[df_name]]$Corr <- as.numeric(list_of_data_frames[[df_name]]$Corr)
    
  }
  
  
  save(file=sprintf("%s\\%s\\%s\\%s", base_path, "samples", "tuning_dataframes", output_file_name), list_of_data_frames)
  
  
  return(tuning_df)
}



for (df_name in names(list_of_data_frames)) {
  df <- list_of_data_frames[[df_name]]
  
  for (subsample in unique(df$Subsamples)) {
    
    subsample_cor_df <- df[df$Subsamples == subsample,]
    subsample_cor_df$Session <- as.numeric(subsample_cor_df$Session) %% 8
    subsample_cor_df$Session[subsample_cor_df$Session == 0] <- 8
    
    session_mean_sd_df <- 
      ddply(subsample_cor_df, 
            .(Subfield), 
            function(subfield_df) {
              by_sim_df <- 
                ddply(subfield_df, 
                      .(isSim),
                      function(sim_df) {
                        by_prop_df <- 
                          ddply(sim_df, 
                                .(Property),
                                function (prop_df) {
                                  by_group_df <- 
                                    ddply(prop_df,
                                          .(Group),
                                          function (group_df) { 
                                            return(ddply(group_df, 
                                                         .(Session), 
                                                         function(sub_df) {
                                                           corr_vec <- 
                                                             as.numeric(sub_df[,"Corr"])
                                                           res <-
                                                             c(mean=mean(corr_vec, 
                                                                         na.rm=T), 
                                                               sd=sem(corr_vec))
                                                           return(res)
                                                         }));
                                          });
                                  return(by_group_df)
                                });
                        return(by_prop_df)});
              return(by_sim_df)});
    
    
    
    group_plots <- list()
    for (group in unique(session_mean_sd_df$Group)) {
      
      property_plots <- list()
      for (property in  unique(session_mean_sd_df$Property)) {
        
        if (property == "Tuning corr") {
          thm <- theme_classic()
        } else {
          thm <- base_plot_theme
        }
        
        working_df <- session_mean_sd_df[session_mean_sd_df$Property == property &
                                           session_mean_sd_df$Group == group, ]
        
        working_df <- working_df[,c("Subfield", "isSim", "Session", "mean", "sd")]
        gtuning <- 
          ggplot(working_df, aes(x=Session, y=mean, group=interaction(isSim, Subfield))) + 
          geom_line(aes(color=interaction(isSim, Subfield)), size=1, linetype="dashed") + 
          geom_errorbar(aes(ymin=mean - sd, ymax=mean + sd, 
                            color=interaction(isSim, Subfield)), 
                        width=0.25, size=1) +
          #ylab(sprintf("Mean correlation (%s)", property)) +
          ylab("Cor (r)") +
          ylim(c(0,0.9)) + 
          ggtitle(property) +
          thm
        
        property_plots <- append(property_plots,
                                 list(gtuning))
      }
      
      property_plots$nrow <- 1
      property_plots$left <- group
      property_plots$widths <- c(rep(1, times=len(unique(session_mean_sd_df$Property)) - 1), 2)
      grp_row <- do.call(grid.arrange, property_plots)
      group_plots <- append(group_plots,
                            list(grp_row))
      
    }
    
    group_plots$nrow <- len(unique(session_mean_sd_df$Group))
    gf <- do.call(grid.arrange, group_plots)
    
    pdf(sprintf("~/../Desktop/learning_analysis_new/%s_subsamples_cor_%s.pdf", df_name, subsample),
        height=len(unique(session_mean_sd_df$Group)) * 3,
        width=len(unique(session_mean_sd_df$Property)) * 3 + 4)
    
    plot(gf)
    
    dev.off()
  }
}

dir.create("~/../Desktop/learning_analysis_new/boxplots/")
tests <- c(wilcox=wilcox.test)

for (test_name in names(tests)) {
  stest <- tests[[test_name]]
  
  for (df_name in names(list_of_data_frames)) {
    corr_df <- list_of_data_frames[[df_name]]
    
    for (subsample in unique(corr_df$Subsamples)) {
      subsample_cor_df <- corr_df[corr_df$Subsamples == subsample,]
      
      plts <- list()
      plts_box <- list()
      
      properties_to_use <- c("peaks",
                             "time_bins",
                             "SI",
                             "spatial_bins")
      #properties_to_use <- unique(subsample_cor_df$Property)[c(1,2,4,6)]
      
      for (subfield in unique(subsample_cor_df$Subfield)) {
        for (prop in properties_to_use)  { 
          
          fdf <- data.frame()
          fdf2 <- data.frame()
          
          pcs_real <- c()
          non_pcs_real <- c()
          pcs_sim <- c()
          non_pcs_sim <- c()
          
          real <- c()
          sim <- c()
          
          pcs <- c()
          non_pcs <- c()
          
          for (mice in unique(subsample_cor_df$Mice)) {
            for (session in unique(subsample_cor_df$Session)){
              for (dir in unique(subsample_cor_df$Direction)) {
                for (group in c("cyclic_pcs", "cyclic_non_pcs")) {
                  
                  
                  working_df <- subsample_cor_df[subsample_cor_df$Subfield == subfield &
                                                   subsample_cor_df$Mice == mice &   
                                                   subsample_cor_df$Session == session &
                                                   subsample_cor_df$Direction == dir &
                                                   subsample_cor_df$Property == prop &
                                                   subsample_cor_df$Group %in%  group,]
                  
                  if (nrow(working_df) < 2) {
                    next
                  }
                  
                  df <- 
                    data.frame(Real_corr=as.numeric(working_df$Corr[working_df$isSim=="Real"]),
                               MOS_corr=as.numeric(working_df$Corr[working_df$isSim=="Simulation"]),
                               Direction=dir,
                               Session=session,
                               Mice=mice,
                               Propery=prop,
                               Group=group)
                  
                  df2 <- 
                    data.frame(Corr=c(as.numeric(working_df$Corr[working_df$isSim=="Real"]),
                                      as.numeric(working_df$Corr[working_df$isSim=="Simulation"])),
                               Direction=c(dir, dir),
                               Session=c(session, session),
                               Mice=c(mice, mice),
                               Propery=c(prop, prop),
                               Group=ifelse(group == "cyclic_pcs", "PCs", "Non-PCs"),
                               isSim=c("Real", "Simulation"))
                  
                  
                  fdf <- rbind(fdf,df)
                  fdf2 <- rbind(fdf2, df2)
                  
                  if (group == "cyclic_pcs") {
                    pcs_real <- c(pcs_real, df2[df2$isSim == "Real","Corr"])
                    pcs_sim <- c(pcs_sim, df2[df2$isSim != "Real","Corr"])
                    
                    pcs <- c(pcs, c(df2[,"Corr"]))
                  } else {
                    non_pcs_real <- c(non_pcs_real, df2[df2$isSim == "Real","Corr"])
                    non_pcs_sim <- c(non_pcs_sim, df2[df2$isSim != "Real","Corr"]) 
                    
                    non_pcs <- c(non_pcs, c(df2[,"Corr"]))
                  }
                  
                  real <- c(real, df2[df2$isSim == "Real","Corr"])
                  sim <- c(sim, df2[df2$isSim != "Real","Corr"])
                  
                  
                }
              }
            }
          }
          
          all_corr_vec <- subsample_cor_df[subsample_cor_df$Property == prop,"Corr"]
          #all_corr_vec <- corr_df$Corr[corr_df$Corr]
          
          # lims <- c(min(all_corr_vec[!is.na(all_corr_vec)]),
          #           max(all_corr_vec[!is.na(all_corr_vec)]))
          
          lims <- c(min(fdf[,c(1,2)]),
                    1)
          
          ldf <- data.frame(x=c(min(lims),max(lims)),
                            y=c(min(lims),max(lims)))
          
          g <- 
            ggplot(fdf, aes(x=Real_corr, y=MOS_corr)) + 
            geom_line(data=ldf,
                      aes(x=x, y=y),
                      linetype="longdash",
                      col="gray40")+
            geom_point(aes(col=Group)) + 
            ylim(lims) + 
            xlim(lims) + 
            xlab("Correlation (real)") + 
            ylab("Correlation (MOS)") + 
            ggtitle(sprintf("%s - %s", subfield, prop)) +
            base_plot_theme 
          
          
          
          
          plts <- append(plts,
                         list(g))
          
          melted_df <- melt(fdf2, value="Corr")
          
          
          
          pr_npr <- stest(pcs_real, non_pcs_real, paired=T)
          ps_nps <- stest(pcs_sim, non_pcs_sim, paired=T)
          ps_pr <- stest(pcs_sim, pcs_real, paired=T)
          nps_npr <- stest(non_pcs_sim, non_pcs_real, paired=T)
          
          main_title <-
            sprintf("prnpr(%.3f), psnps(%.3f), pspr(%.3f), npsnpr(%.3f)",
                    pr_npr$p.value,
                    ps_nps$p.value,
                    ps_pr$p.value,
                    nps_npr$p.value)
          
          gb <- ggplot(melted_df, aes(x=Group, y=value, fill=isSim)) +
            geom_boxplot() +
            theme_classic() +
            ggtitle(main_title) + 
            theme(plot.title = element_text(size=10, face="bold")) +
            xlab("") + 
            ylab(sprintf("Correlation (%s)", prop))
          
          plts_box <- append(plts_box,
                             list(gb))
        }
      }
      
      plts$nrow <- 2
      gf <- do.call(grid.arrange, plts)
      
      pdf(file=sprintf("~/../Desktop/learning_analysis_new/boxplots/%s_corlineplots_s%s.pdf", df_name, subsample),
          height=4,
          width=8)
      plot(gf)
      dev.off()
      
      plts_box$nrow <- 2
      gf_box <- do.call(grid.arrange, plts_box)
      
      
      pdf(file=sprintf("~/../Desktop/learning_analysis_new/boxplots/%s_boxplots_s%s_%s.pdf", df_name, subsample, test_name),
          height=6,
          width=18)
      plot(gf_box)
      dev.off()
      
    }
    
  }
}

for (df_name in names(list_of_data_frames)) {
  corr_df <- list_of_data_frames[[df_name]]
  
  for (subsample in unique(corr_df$Subsamples)) {
    subsample_cor_df <- corr_df[corr_df$Subsamples == subsample,]
    
    metric_plot_lists <- list()
    
    properties_to_use <- c("peaks",
                           "time_bins",
                           "SI",
                           "spatial_bins")
    metric_functions <- list(KL=KL_dist,
                             JSD=JSD_dist,
                             dprime=dprime_dist,
                             wasserstein=wasserstein_dist)
    
    for (metric_name in names(metric_functions)) {
      metric_plot_lists[[metric_name]] <- list()
    }
    
    #properties_to_use <- unique(subsample_cor_df$Property)[c(1,2,4,6)]
    
    
    subfield_list <- list()
    
    for (subfield in unique(subsample_cor_df$Subfield)) {
      
      
      
      for (prop in properties_to_use)  { 
        
        
        comparisions_df <- data.frame()
        subsample_cor_df$Session <- as.numeric(subsample_cor_df$Session)
        
        for (session in 1:16){
          pcs_real <- c()
          non_pcs_real <- c()
          pcs_sim <- c()
          non_pcs_sim <- c()
          
          for (mice in unique(subsample_cor_df$Mice)) {
            for (env in c(0)) {
              for (dir in unique(subsample_cor_df$Direction)) {
                for (group in c("cyclic_pcs", "cyclic_non_pcs")) {
                  
                  
                  working_df <- subsample_cor_df[subsample_cor_df$Subfield == subfield &
                                                   subsample_cor_df$Mice == mice &   
                                                   subsample_cor_df$Session == (session + env) &
                                                   subsample_cor_df$Direction == dir &
                                                   subsample_cor_df$Property == prop &
                                                   subsample_cor_df$Group %in%  group,]
                  
                  if (nrow(working_df) < 2) {
                    next
                  }
                  
                  
                  df2 <- 
                    data.frame(Corr=c(as.numeric(working_df$Corr[working_df$isSim=="Real"]),
                                      as.numeric(working_df$Corr[working_df$isSim=="Simulation"])),
                               Direction=c(dir, dir),
                               Session=c(session, session),
                               Mice=c(mice, mice),
                               Propery=c(prop, prop),
                               Group=ifelse(group == "cyclic_pcs", "PCs", "Non-PCs"),
                               isSim=c("Real", "Simulation"))
                  
                  
                  if (group == "cyclic_pcs") {
                    pcs_real <- c(pcs_real, df2[df2$isSim == "Real","Corr"])
                    pcs_sim <- c(pcs_sim, df2[df2$isSim != "Real","Corr"])
                    
                  } else {
                    non_pcs_real <- c(non_pcs_real, df2[df2$isSim == "Real","Corr"])
                    non_pcs_sim <- c(non_pcs_sim, df2[df2$isSim != "Real","Corr"]) 
                    
                  }
                  
                }
              }
            }  
          }
          
          
          for (metric_name in names(metric_functions)) {
            
            
            func <- metric_functions[[metric_name]]
            
            
            comparisions_vector <- 
              c(`PC_REAL_PC_SIM`=func(pcs_real, pcs_sim),
                `NPC_REAL_NPC_SIM`=func(non_pcs_real, non_pcs_sim),
                `PC_SIM_NPC_SIM`=func(pcs_sim, non_pcs_sim),
                `PC_REAL_NPC_REAL`=func(pcs_real, non_pcs_real))
            
            for (comp in names(comparisions_vector))  {
              
              df <- data.frame(dist=comparisions_vector[comp],
                               metric=metric_name,
                               session=session,
                               group=comp)
              
              comparisions_df <- rbind(comparisions_df,
                                       df)
            }
            
            
            
          }
          
          
        }
        
        
        
        for (metric_name in unique(comparisions_df$metric)) {
          work_df <- comparisions_df[comparisions_df$metric == metric_name,]
          
          grp_names <- unique(work_df$group)
          pval_mt <- matrix(rep(1, times=len(grp_names) ** 2),
                            nrow=len(grp_names))
          rownames(pval_mt) <- grp_names
          colnames(pval_mt) <- grp_names
          
          work_df$group <- factor(work_df$group, levels=c("NPC_REAL_NPC_SIM",
                                                          "PC_REAL_NPC_REAL",
                                                          "PC_REAL_PC_SIM",
                                                          "PC_SIM_NPC_SIM"))
          
          for (grp in grp_names) {
            for (grp2 in grp_names) {
              if (grp != grp2) {
                
                g_df <- work_df[work_df$group == grp,]
                g2_df <- work_df[work_df$group == grp2,]
                
                wx <- 
                  wilcox.test(g_df[order(g_df$session),"dist"],
                              g2_df[order(g2_df$session),"dist"],
                              paired=T)
                
                pval_mt[grp, grp2] <- wx$p.value
              }
            }
          }
          
          max_dist <- max(work_df$dist)
          
          
          
          
          g <- 
            ggplot(work_df, aes(x=group, y=dist)) +
            geom_jitter(aes(color=group)) + 
            geom_boxplot(alpha=0.5, aes(fill=group)) + 
            base_plot_theme + 
            theme(axis.text.x  = element_text(angle = 45, vjust = 1, hjust=1)) + 
            ylab(sprintf("Distance (%s)", metric_name)) + 
            xlab("") +
            ggtitle(sprintf("%s - %s", subfield, prop)) +
            ylim(c(0, max_dist * 2.15)) +
            geom_line(data=data.frame(x=c(1,1,2,2),
                                      y=c(max_dist * 1.2,
                                          max_dist * 1.25, 
                                          max_dist * 1.25,
                                          max_dist * 1.2)),
                      aes(x=x, y=y), color="black") +
            geom_line(data=data.frame(x=c(3,3,4,4),
                                      y=c(max_dist * 1.2,
                                          max_dist * 1.25, 
                                          max_dist * 1.25,
                                          max_dist * 1.2)),
                      aes(x=x, y=y), color="black") +
            geom_line(data=data.frame(x=c(1,1,3,3),
                                      y=c(max_dist * 1.7,
                                          max_dist * 1.75, 
                                          max_dist * 1.75,
                                          max_dist * 1.7)),
                      aes(x=x, y=y), color="black") +
            geom_line(data=data.frame(x=c(2,2,4,4),
                                      y=c(max_dist * 1.95,
                                          max_dist * 2, 
                                          max_dist * 2,
                                          max_dist * 1.95)),
                      aes(x=x, y=y), color="black") +
            geom_text(x=1.5,
                      y=max_dist * 1.32,
                      label=pval_mt["NPC_REAL_NPC_SIM", "PC_REAL_NPC_REAL"]) +
            geom_text(x=3.5,
                      y=max_dist * 1.32,
                      label=pval_mt["PC_REAL_PC_SIM", "PC_SIM_NPC_SIM"]) +
            geom_text(x=2,
                      y=max_dist * 1.82,
                      label=pval_mt["NPC_REAL_NPC_SIM", "PC_SIM_NPC_SIM"]) +
            geom_text(x=3,
                      y=max_dist * 2.07,
                      label=pval_mt["PC_REAL_NPC_REAL", "PC_SIM_NPC_SIM"])
          
          
          metric_plot_lists[[metric_name]][[prop]] <- g
        }
      }
      
      subfield_list <- append(subfield_list,
                              list(metric_plot_lists))
    }
    
    for (metric_name in names(metric_functions)) {
      
      final_plt_list <- list() 
      for (sbf in c(1:len(subfield_list))) {
        final_plt_list <- append(final_plt_list,
                                 subfield_list[[sbf]][[metric_name]])
      }
      
      final_plt_list$nrow <- 2
      gf <- do.call(arrangeGrob, final_plt_list)
      dir.create(sprintf("~/../Desktop/learning_analysis_new/comparisions_16/%s", metric_name))
      pdf(file=sprintf("~/../Desktop/learning_analysis_new/comparisions_16/%s/comp_%s%s.pdf", metric_name, df_name, subsample),
          height=9,
          width=12)
      plot(gf)
      dev.off()
    }
  }
}


for (sbf in c("CA1", "CA3")) {
  plt_all <- list()
  
  # 
  # for (group in unique(properties_df$Group)) {
  #   
  #   properties <- c("Peaks", "Time bins",
  #                   "SI", "Spatial bins")
  #   group_prop_df <- properties_df[properties_df$Group == group,]
  #   
  #   group_prop_df$Session <- as.numeric(group_prop_df$Session) %% 8 
  #   group_prop_df$Session[group_prop_df$Session == 0] <- 8
  #   
  #   fdf <- 
  #   ddply(group_prop_df, .(Subfield),
  #         function(subfield_df) {
  #           ses_df <- 
  #           ddply(subfield_df, .(Session),
  #                 function(by_ses_df) {
  #                   
  #                   sim_df <- ddply(by_ses_df, .(isSim),
  #                                   function(by_sim_df) {
  #                                     
  #                                     mn <- colMeans(by_sim_df[,properties])
  #                                     psd <- apply(by_sim_df[,properties], 2, sd)
  #                                     mnf <- c(mn, psd, by_sim_df[1,c("Session", "Subfield", "Group", "isSim")])
  #                                     
  #                                     
  #                                     mnf <- unlist(mnf)
  #                                     mnf_names <- c(sprintf("%s %s", properties, c("mean")),
  #                                                    sprintf("%s %s", properties, c("sd")),
  #                                                    "Session", "Subfield", "Group", "isSim")
  # 
  #                                     names(mnf) <- mnf_names
  #                                     
  #                                   
  #                                     return(mnf)
  #                                     
  #                                   })
  #                   
  #                   return(sim_df)
  #                 })
  #           
  #           return(ses_df)
  #         })
  #   
  #   
  #   plt_list <- list()
  #   
  #   fdf2 <- fdf[fdf$Subfield == sbf,] 
  #   
  #   for (prop in properties) {
  #     
  #     work_df <- fdf2[,c(sprintf("%s %s", prop, c("mean", "sd")), c("Session", "Subfield", "Group", "isSim"))]
  #     
  #     colnames(work_df) <- c("Mean", "Sd", "Session", "Subfield", "Group", "isSim")
  #     work_df$Session <- as.numeric(work_df$Session)
  #     work_df$Mean <- as.numeric(work_df$Mean)
  #     work_df$Sd <- as.numeric(work_df$Sd)
  #     
  #     g <- 
  #       ggplot(work_df, aes(x=Session, y=Mean, group=interaction(isSim, Subfield))) + 
  #       geom_line(aes(color=interaction(isSim, Subfield)), size=1, linetype="dashed") + 
  #        geom_errorbar(aes(ymin=Mean - Sd, ymax=Mean + Sd, 
  #                          color=interaction(isSim, Subfield)), 
  #                      width=0.25, size=1) +
  #       #ylab(sprintf("Mean correlation (%s)", property)) +
  #       ylim(c(0.8 * min(properties_df[which(!is.na(as.numeric(properties_df[,prop]))),prop]),
  #              1.2 * max(properties_df[which(!is.na(as.numeric(properties_df[,prop]))),prop]))) +
  #       ylab(prop) +
  #       ggtitle(prop) +
  #       base_plot_theme
  #     
  #     plt_list <- append(plt_list,
  #                        list(g))
  #   }
  #   
  #   plt_list$nrow <- 1
  #   plt_list$left <- group
  #   gf <- do.call(grid.arrange, plt_list)
  #   
  #   plt_all <- append(plt_all,
  #                     list(gf))
  # }
  # 
  # 
  # plt_all$nrow <- 7
  # pltff <- do.call(grid.arrange, plt_all)
  # 
  # 
  # pdf(sprintf("~/../Desktop/learning_analysis_figures/line_plot_%s.pdf", sbf),
  #     width= 4 * 2.5,
  #     height = 7 * 2.5) 
  # plot(pltff)
  # dev.off()
  
  
  
  plt_box_all <- list()
  plt_box_by_session <- list()
  
  for (group in unique(properties_df$Group)) {
    
    properties <- c("Peaks", "Time bins",
                    "SI", "Spatial bins")
    group_prop_df <- properties_df[properties_df$Group == group,]
    
    
    plt_list <- list()
    plt_list_by_sess <- list()
    
    wilk_list <- list()
    
    fdf2 <- group_prop_df[group_prop_df$Subfield == sbf,] 
    
    pval_mat <- c()
    
    for (prop in properties) {
      
      work_df <- fdf2[,c(prop, c("Session", "Subfield", "Group", "isSim", "Mice", "Direction"))]
      
      colnames(work_df) <- c("Property","Session", "Subfield", "Group", "isSim", "Mice", "Direction")
      work_df$Session <- as.numeric(work_df$Session)
      work_df$Property <- as.numeric(work_df$Property)
      
      
      real_v <- c()
      sim_v <- c()
      
      session_pvalues <- list()
      
      for (s in unique(work_df$Session)) {
        tmp_rv <- c()
        tmp_sv <- c()
        
        for (m in unique(work_df$Mice)) {
          for (d in unique(work_df$Direction)) {
            rv <- 
              work_df[work_df$Mice == m &
                        work_df$Session == s &
                        work_df$Direction == d &
                        work_df$isSim == "Real",
                      "Property"]
            
            sv <- 
              work_df[work_df$Mice == m &
                        work_df$Session == s &
                        work_df$Direction == d &
                        work_df$isSim != "Real",
                      "Property"]
            
            
            real_v <- c(real_v, rv)
            sim_v <- c(sim_v, sv)
            
            tmp_rv <- c(tmp_rv, rv)
            tmp_sv <- c(tmp_sv, sv)
            
            
          }
        }
        
        session_pvalues <- append(session_pvalues,
                                  list(list(tt=t.test(tmp_sv, tmp_rv, paired=T),
                                            wx=wilcox.test(tmp_sv, tmp_rv, paired=T))))
      }
      
      pvals <- 
        do.call(rbind,lapply(session_pvalues, function(s) {c(s$tt$p.value, s$wx$p.value)}))
      pval_mat <- cbind(pval_mat, 
                        pvals)
      
      work_df$Session <- factor(work_df$Session, levels=1:16)  
      
      
      
      g_by_sess <- 
        ggplot(work_df, aes(x=Session, y=Property, fill=isSim)) + 
        geom_boxplot(aes(x=Session, y=Property))+
        ylim(c(0.8 * min(properties_df[which(!is.na(as.numeric(properties_df[,prop]))),prop]),
               1.2 * max(properties_df[which(!is.na(as.numeric(properties_df[,prop]))),prop]))) +
        ylab(prop) +
        ggtitle(prop) +
        base_plot_theme
      
      wx <- 
        wilcox.test(real_v,
                    sim_v,
                    paired=T)
      
      tt <- t.test(real_v, sim_v, paired = T)
      
      g <- 
        ggplot(work_df, aes(x=isSim, y=Property, fill=isSim)) + 
        geom_boxplot(aes(x=isSim, y=Property))+
        ylim(c(0.8 * min(properties_df[which(!is.na(as.numeric(properties_df[,prop]))),prop]),
               1.2 * max(properties_df[which(!is.na(as.numeric(properties_df[,prop]))),prop]))) +
        ylab(prop) +
        ggtitle(sprintf("w%.4f, t%.4f",
                        wx$p.value,
                        tt$p.value)) +
        base_plot_theme
      
      plt_list <- append(plt_list,
                         list(g))
      
      plt_list_by_sess <- append(plt_list_by_sess,
                                 list(g_by_sess))
    }
    
    colnames(pval_mat) <- paste(rep(properties, each=2), c("T-test", "Wilcox"))
    rownames(pval_mat) <- sprintf("Session %d", 1:nrow(pval_mat))
    
    brks <- c(0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.45, 0.65, 0.85, 1)
    ph_pv <- 
      pheatmap(pval_mat, cluster_rows=F,
               cluster_cols=F, 
               breaks=brks, 
               col=rev(brewer.pal(len(brks) - 1, "RdYlBu")))
    
    plt_list_by_sess <- append(plt_list_by_sess,
                               list(ph_pv[[4]]))
    plt_list$nrow <- 1
    plt_list$left <- group
    gf <- do.call(grid.arrange, plt_list)
    
    plt_list_by_sess$nrow <- 1
    plt_list_by_sess$left <- group
    gf_by_sess <- do.call(grid.arrange, plt_list_by_sess)
    
    plt_box_all <- append(plt_box_all,
                          list(gf))
    
    plt_box_by_session <- append(plt_box_by_session,
                                 list(gf_by_sess))
    
  }
  
  
  plt_box_all$nrow <- 7
  pltff <- do.call(grid.arrange, plt_box_all)
  
  plt_box_by_session$nrow <- 7
  pltff_by_sess <- do.call(grid.arrange, plt_box_by_session)
  
  
  pdf(sprintf("~/../Desktop/learning_analysis_figures/boxplot_both_dir_new_%s.pdf", sbf),
      width= 4 * 2.5,
      height = 7 * 2.5) 
  plot(pltff)
  dev.off()
  
  pdf(sprintf("~/../Desktop/learning_analysis_figures/boxplot_both_dir_by_session_%s.pdf", sbf),
      width= 5 * 3,
      height = 7 * 3) 
  plot(pltff_by_sess)
  dev.off()
  
}





for (subsample in unique(correlations_df$Subsamples)) {
  subsample_cor_df <- correlations_df[correlations_df$Subsamples == subsample,]
  subsample_cor_df$Session <- as.numeric(subsample_cor_df$Session) %% 8
  subsample_cor_df$Session[subsample_cor_df$Session == 0] <- 8
  session_mean_sd_df <- 
    ddply(subsample_cor_df, 
          .(Subfield), 
          function(subfield_df) {
            by_sim_df <- 
              ddply(subfield_df, 
                    .(isSim),
                    function(sim_df) {
                      by_prop_df <- 
                        ddply(sim_df, 
                              .(Property),
                              function (prop_df) {
                                by_group_df <- 
                                  ddply(prop_df,
                                        .(Group),
                                        function (group_df) { 
                                          return(ddply(group_df, 
                                                       .(Session), 
                                                       function(sub_df) {
                                                         corr_vec <- 
                                                           as.numeric(sub_df[,"Corr"])
                                                         res <-
                                                           c(mean=mean(corr_vec, 
                                                                       na.rm=T), 
                                                             sd=sem(corr_vec))
                                                         return(res)
                                                       }));
                                        });
                                return(by_group_df)
                              });
                      return(by_prop_df)});
            return(by_sim_df)});
  
  
  
  
  group_plots <- list()
  for (group in unique(session_mean_sd_df$Group)) {
    
    property_plots <- list()
    for (property in  unique(session_mean_sd_df$Property)) {
      
      if (property == "Tuning corr") {
        thm <- theme_classic()
      } else {
        thm <- base_plot_theme
      }
      
      working_df <- session_mean_sd_df[session_mean_sd_df$Property == property &
                                         session_mean_sd_df$Group == group, ]
      
      working_df <- working_df[,c("Subfield", "isSim", "Session", "mean", "sd")]
      gtuning <- 
        ggplot(working_df, aes(x=Session, y=mean, group=interaction(isSim, Subfield))) + 
        geom_line(aes(color=interaction(isSim, Subfield)), size=1, linetype="dashed") + 
        geom_errorbar(aes(ymin=mean - sd, ymax=mean + sd, 
                          color=interaction(isSim, Subfield)), 
                      width=0.25, size=1) +
        #ylab(sprintf("Mean correlation (%s)", property)) +
        ylab("Cor (r)") +
        ylim(c(0,0.9)) + 
        ggtitle(property) +
        thm
      
      property_plots <- append(property_plots,
                               list(gtuning))
    }
    
    property_plots$nrow <- 1
    property_plots$left <- group
    property_plots$widths <- c(rep(1, times=len(unique(session_mean_sd_df$Property)) - 1), 2)
    grp_row <- do.call(grid.arrange, property_plots)
    group_plots <- append(group_plots,
                          list(grp_row))
    
  }
  
  group_plots$nrow <- len(unique(session_mean_sd_df$Group))
  gf <- do.call(grid.arrange, group_plots)
  
  pdf(sprintf("~/../Desktop/subsamples_cor_%s.pdf", subsample),
      height=len(unique(session_mean_sd_df$Group)) * 3,
      width=len(unique(session_mean_sd_df$Property)) * 3 + 4)
  
  plot(gf)
  
  dev.off()
}


