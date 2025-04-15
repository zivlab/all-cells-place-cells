
figure_4_boxplots  <- function (directions_to_use = c("Both", "Left", "Right")) {
  
  write_path <- sprintf("%s\\figure_4\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_4\\boxplots\\",figures_path)
  
  sizes=c(big=3.5,
          medium=3,
          small=2.5)
  
  dir.create(write_path)
  dir.create(sprintf("%s\\estimations", write_path))
  dir.create(sprintf("%s\\estimations\\big\\", write_path))
  dir.create(sprintf("%s\\estimations\\medium\\", write_path))
  dir.create(sprintf("%s\\estimations\\small\\", write_path))
  dir.create(sprintf("%s\\naive_measurements", write_path))
  dir.create(sprintf("%s\\naive_measurements\\big\\", write_path))
  dir.create(sprintf("%s\\naive_measurements\\medium\\", write_path))
  dir.create(sprintf("%s\\naive_measurements\\small\\", write_path))
  
  true_data_paths <- sprintf("%s\\%s", all_data_paths[1:18], "equalized")
  
  sessions_to_use=c(5:8,13:16)
  
  color_plot <- function(g, x_axes) {
    gf <-  g + 
      scale_fill_manual(values=
                          c(adjustcolor(c("royalblue2", "gray80"), alpha=0.8),
                            adjustcolor(c("royalblue2", "gray80"), alpha=0.8),
                            adjustcolor(c("royalblue2", "gray80"), alpha=0.8)),
                        breaks=c("Both.CA3", "Both.CA1",
                                 "Left.CA3", "Left.CA1",
                                 "Right.CA3", "Right.CA1")) +
      scale_color_manual(values=
                           c("royalblue3","gray65",
                                         "royalblue3","gray65",
                                         "royalblue3","gray65"),
                                         breaks=c("Both.CA3", "Both.CA1",
                                                  "Left.CA3", "Left.CA1",
                                                  "Right.CA3", "Right.CA1")) + 
      base_plot_theme + 
      theme(text=element_text(size=15)) +
      xlab(x_axes) +
      ylab("Fraction") + 
      ylim(0,1.05)
    
    return(gf)
  }
  
  all_df <- get_all_df(true_data_paths, sessions_to_use)
  
  estimation_configurations <- list(KS_cyclic=list(path="KS_simulations_likelihood_c", edge_free=F),
                                    KS_random=list(path="KS_simulations_likelihood_r", edge_free=F),
                                    JSD_cyclic=list(path="JSD_simulations_likelihood_c", edge_free=F),
                                    JSD_random=list(path="JSD_simulations_likelihood_r", edge_free=F),
                                    JSD_cyclic_edge_free=list(path="JSD_simulations_likelihood_c", edge_free=T),
                                    JSD_random_edge_free=list(path="JSD_simulations_likelihood_r", edge_free=T))

  estimated_statistics_df <- data.frame()
    
  for (conf_name in names(estimation_configurations)) {
    
    conf = estimation_configurations[[conf_name]]
    estimated_res <- get_estimated_df(all_df, true_data_paths, sessions_to_use,
                                      estimated_df_folder = conf$path,
                                      edge_free=conf$edge_free)
    
    both_dir_res <- get_both_dir_df(estimated_df = estimated_res$estimated_df, all_df=all_df)
    
    both_dir_df <- both_dir_res$both_dir_df
    
    #### statistics ##### 
    
    for (direc in unique(both_dir_df$Direction)) {
      ca1_df <- both_dir_df[both_dir_df$Direction == direc & 
                            both_dir_df$Subfield == "CA1",]
      
      ca3_df <- both_dir_df[both_dir_df$Direction == direc & 
                              both_dir_df$Subfield == "CA3",]
      
      wilc <- wilcox.test(ca1_df$Estimated, ca3_df$Estimated, alternative="greater")
      estimated_statistics_df <- rbind(estimated_statistics_df,
                                       data.frame(W=wilc$statistic,
                                                  Pv=wilc$p.value,
                                                  alternative=wilc$alternative,
                                                  conf_name=conf_name,
                                                  signif_code=signif.num(wilc$p.value),
                                                  direction=direc))
    }
    
    mean_sd_df <- both_dir_res$mean_sd_df
    session_mean_sd_df <- both_dir_res$session_mean_sd_df
    
    mean_sd_df <- mean_sd_df[mean_sd_df$Direction %in% directions_to_use, ]
    both_dir_df <- both_dir_df[both_dir_df$Direction %in% directions_to_use, ]
    
    if (all(directions_to_use == "Both")) {
      g_est <- ggplot(mean_sd_df, aes(x=Subfield, y=Estimated)) 
      x_axes = "Subfield"
    } else{
      g_est <- ggplot(mean_sd_df, aes(x=Direction, y=Estimated))
      x_axes = "Running direction"
    }
    g_est <- 
      g_est +
      geom_hline(yintercept=1, linetype="dashed", size=0.8) +
      geom_boxplot(data=both_dir_df, 
                   aes(group=interaction(Direction,Subfield),
                       color=interaction(Direction,Subfield),
                       fill=interaction(Direction,Subfield)), 
                   width=0.85,
                   position = position_dodge(1),
                   size=1) + 
      geom_jitter(data=both_dir_df,
                  aes(y=Estimated, 
                      group=interaction(Direction,Subfield), 
                      fill=interaction(Direction,Subfield)),
                  position=position_jitterdodge(1), 
                  color="gray20", 
                  size=0.75, 
                  alpha=0.1) 
    
    g_est <- color_plot(g_est, x_axes) + ggtitle("Place cells (estimated)") +
      theme(plot.title=element_text(size=11))
    
    
    for (size_name in names(sizes)) {
      size = sizes[[size_name]]
      pdf(file=sprintf("%s\\estimations\\%s\\%s.pdf",
                       write_path,
                       size_name,
                       conf_name),
          height=size,
          width=size)
      
      plot(g_est)
      dev.off()
    }
   }

  write.csv(estimated_statistics_df,
            file=sprintf("%s\\estimated_measurements_statistics.csv",write_path))
    
  conf = estimation_configurations[[1]]
  
  naive_configurations <- list(cyclic_eq_prior_weighted=list(path="percent_df_cyclic_SI_eq_prior.R", weighted=T),
                               cyclic_eq_prior_either_dir=list(path="percent_df_cyclic_SI_eq_prior.R", weighted=F),
                               random_eq_prior_weighted=list(path="percent_df_rand_SI_eq_prior.R", weighted=T),
                               random_eq_prior_either_dir=list(path="percent_df_rand_SI_eq_prior.R", weighted=F),
                               cyclic_non_eq_prior_weighted=list(path="percent_df_cyclic_SI_reg.R", weighted=T),
                               cyclic_non_eq_prior_either_dir=list(path="percent_df_cyclic_SI_reg.R", weighted=F),
                               random_non_eq_prior_weighted=list(path="percent_df_rand_SI_reg.R", weighted=T),
                               random_non_eq_prior_either_dir=list(path="percent_df_rand_SI_reg.R", weighted=F))
  
  naive_statistics_df <- data.frame()
  
  for (naive_conf_name in names(naive_configurations)) { 
    naive_conf = naive_configurations[[naive_conf_name]]
    
    all_df <- get_all_df(true_data_paths, 
                         sessions_to_use,
                         file_name = naive_conf$path)
    
    estimated_res <- get_estimated_df(all_df, true_data_paths, sessions_to_use)
    
    both_dir_res <- get_both_dir_df(estimated_df = estimated_res$estimated_df, all_df=all_df)
    
    both_dir_df <- both_dir_res$both_dir_df
    mean_sd_df <- both_dir_res$mean_sd_df
    session_mean_sd_df <- both_dir_res$session_mean_sd_df
    
    mean_sd_df <- mean_sd_df[mean_sd_df$Direction %in% directions_to_use, ]
    both_dir_df <- both_dir_df[both_dir_df$Direction %in% directions_to_use, ]
    
    percent_active_df <- all_df[,c(7:9)]
    percent_place_cells_df <- all_df[,c(1:3)]
    percent_all_df <- all_df[,c(4:6)]
    
    metadata_df <- all_df[,c("Session", "Mice", "Subfield")]
    percent_active_df <- cbind(percent_active_df, metadata_df)
    percent_all_df <- cbind(percent_all_df, metadata_df)
    
    percent_place_cells_df <- cbind(percent_place_cells_df, metadata_df)  
    colnames(percent_active_df) <- c("Both", "Right", "Left", "Session", "Mice", "Subfield")
    colnames(percent_place_cells_df) <- c("Both", "Right", "Left", "Session", "Mice", "Subfield")
    colnames(percent_all_df) <- c("Both", "Right", "Left", "Session", "Mice", "Subfield")
    
    melted_percent_active <- melt(percent_active_df, id.vars = c("Session", "Mice", "Subfield"))
    melted_percent_place_cells <- melt(percent_place_cells_df, id.vars = c("Session", "Mice", "Subfield"))
    melted_percent_all <- melt(percent_all_df, id.vars = c("Session", "Mice", "Subfield"))
    
    ylabs = c("Percentage of active cells (%)",
              "Percentage of place cells (% of active cells)",
              "Percentage of place cells (% of all cells)")
    
    dfs <- list(melted_percent_active,
                melted_percent_place_cells,
                melted_percent_all)
    
    weighted_avg_both_df_colnames <- c("Active",
                                       "Of_Active",
                                       "Of_All")
    
    plot_names <- c("ActiveCells",
                   "PlaceCellsOfActive",
                   "PlaceCellsOfAll")
    
    
    plot_titles <- c("Active cells",
                     "Place cells (of active)",
                     "Place cells (of all)")
    
    
    for (df_idx in 1:3) {
      df_to_use <- dfs[[df_idx]]
      colnames(df_to_use) <- c("Session", "Mice", "Subfield", "Direction", "Percent")
      
      # Set second env session count refrence to 1...8
      df_to_use$Session <- df_to_use$Session %% 8
      df_to_use$Session[df_to_use$Session == 0] <- 8
      
      for (direc in unique(df_to_use$Direction)) {
        ca1_df <- df_to_use[df_to_use$Direction == direc & 
                            df_to_use$Subfield == "CA1",]
        
        ca3_df <- df_to_use[df_to_use$Direction == direc & 
                            df_to_use$Subfield == "CA3",]
        
        wilc <- wilcox.test(ca1_df$Percent, ca3_df$Percent, alternative="greater")
        naive_statistics_df <- rbind(naive_statistics_df,
                                     data.frame(W=wilc$statistic,
                                                Pv=wilc$p.value,
                                                alternative=wilc$alternative,
                                                conf_name=naive_conf_name,
                                                signif_code=signif.num(wilc$p.value),
                                                direction=direc,
                                                category=plot_titles[df_idx]))
      }
      
      if (naive_conf$weighted) {
        # Remove both 
        df_to_use <- 
          df_to_use[df_to_use$Direction != "Both",]
        
        weighted_both_df <- 
          both_dir_res$both_dir_df_all[,c("Session", "Mice", "Subfield", weighted_avg_both_df_colnames[df_idx])]
        
        # Reorganize column names
        weighted_both_df$Direction <- rep("Both", nrow(weighted_both_df))
        weighted_both_df <- weighted_both_df[,c("Session", "Mice", "Subfield", "Direction", weighted_avg_both_df_colnames[df_idx])]
        
        colnames(weighted_both_df) <- c("Session", "Mice", "Subfield", "Direction", "Percent")
        df_to_use <- rbind(df_to_use, weighted_both_df)
      }
      
      
      df_to_use <- df_to_use[df_to_use$Direction %in% directions_to_use, ]
      
      sd_df <- ddply(df_to_use, .(Subfield),
                     function(subfield_df) {
                       return(ddply(subfield_df, .(Direction), 
                                    function(sub_df) {
                                      return(c(mean(sub_df[,"Percent"]), sd(sub_df[,"Percent"])))}))
                     })
      colnames(sd_df) <- c("Subfield", "Direction", "Percent", "Sd")
      
      
      if (all(directions_to_use == "Both")) {
        g <- ggplot(sd_df, aes(x=Subfield, y=Percent))
        x_axes = "Subfield"
      } else{
        g <- ggplot(sd_df, aes(x=Direction, y=Percent)) 
        x_axes = "Running direction"
      }
      
      g <- 
        g +
        geom_hline(yintercept=1, linetype="dashed", size=0.8) +
        geom_boxplot(data=df_to_use, 
                     aes(group=interaction(Direction,Subfield),
                         color=interaction(Direction,Subfield),
                         fill=interaction(Direction,Subfield)), 
                     width=0.85,
                     position = position_dodge(1),
                     size=1) + 
        geom_jitter(data=df_to_use,
                    aes(y=Percent, 
                        group=interaction(Direction,Subfield), 
                        fill=interaction(Direction,Subfield)),
                    position=position_jitterdodge(1), 
                    color="gray20", 
                    size=0.75, 
                    alpha=0.1) 
      
      g <- color_plot(g, x_axes) + ggtitle(plot_titles[df_idx]) +
        theme(plot.title=element_text(size=11))
        
       
      
      for (size_name in names(sizes)) {
        size = sizes[[size_name]]
        
        dir.create(sprintf("%s\\naive_measurements\\%s\\%s\\",
                           write_path,
                           size_name,
                           plot_names[[df_idx]]))
        
        pdf(file=sprintf("%s\\naive_measurements\\%s\\%s\\%s.pdf",
                         write_path,
                         size_name,
                         plot_names[[df_idx]],
                         naive_conf_name),
            height=size,
            width=size)
        
        plot(g)
        dev.off()
      }
    }
  }
  
  write.csv(naive_statistics_df,
            file=sprintf("%s\\naive_measurements_statistics.csv",write_path))
}

figure_4_likelihood_plots <- function() {

  write_path <- sprintf("%s\\figure_4\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_4\\likelihood_plots\\",figures_path)
  dir.create(write_path)
  dir.create(sprintf("%s\\figure_4\\likelihood_plots\\big\\",figures_path))
  dir.create(sprintf("%s\\figure_4\\likelihood_plots\\small\\",figures_path))
  dir.create(sprintf("%s\\figure_4\\likelihood_plots\\medium\\",figures_path))
  
  sizes=c(big=3,
          medium=2.5,
          small=2)
  
  dir.create(write_path)
  
  true_data_paths <- sprintf("%s\\%s", all_data_paths[1:18], "equalized")
  
  sessions_to_use=c(5:8,13:16)
  
  all_df <- get_all_df(true_data_paths, sessions_to_use)
  
  estimation_configurations <- list(KS_cyclic=list(path="KS_simulations_likelihood_c", edge_free=F),
                                    KS_random=list(path="KS_simulations_likelihood_r", edge_free=F),
                                    JSD_cyclic=list(path="JSD_simulations_likelihood_c", edge_free=F),
                                    JSD_random=list(path="JSD_simulations_likelihood_r", edge_free=F),
                                    JSD_cyclic_edge_free=list(path="JSD_simulations_likelihood_c", edge_free=T),
                                    JSD_random_edge_free=list(path="JSD_simulations_likelihood_r", edge_free=T))
  
  
  for (conf_name in names(estimation_configurations)) {
    
    conf = estimation_configurations[[conf_name]]
    estimated_res <- get_estimated_df(all_df, true_data_paths, sessions_to_use,
                                      estimated_df_folder = conf$path,
                                      edge_free=conf$edge_free)
    
    likelihood_df_all <- estimated_res$likelihood_df
    
    metadata_df <- do.call(rbind, estimated_res$session_metavar_list)
    
    ind_list <- list(CA1=which(metadata_df[,4] == "CA1"),
                     CA3=which(metadata_df[,4] == "CA3"))
    
    for (subfield_name in names(ind_list)) {
      ind <- ind_list[[subfield_name]]
      
      likdf_tmp <- likelihood_df_all[ind,]
      likelihood_df <- (t(apply(likdf_tmp, 1, function(r) {return(r/sum(r))})))
      likelihood_df <- as.data.frame(likelihood_df)
      likelihood_sessions <- as.numeric(metadata_df[ind,1])
      likelihood_sessions[likelihood_sessions > 8] <-  likelihood_sessions[likelihood_sessions > 8] - 8
      likelihood_df <- cbind(likelihood_df, likelihood_sessions)
      colnames(likelihood_df) <- c(as.character(seq(0.5,1,by=0.1)), "sessions")
      melted_likelihood <- melt(likelihood_df, id.vars="sessions")
      colnames(melted_likelihood) <- c("Session", "Simulated percent", "Normalized likelihood")
      
      likelihood_mean_sd_df <- ddply(melted_likelihood, .(`Simulated percent`), 
                                     function(sub_df) {
                                       return(c(mean(sub_df[,"Normalized likelihood"], na.rm=T),
                                                sd(sub_df[,"Normalized likelihood"], na.rm=T)))
                                     })
      
      colnames(likelihood_mean_sd_df) <- c("Simulated percent", "Mean", "Sd")
      # 
      # 
      # glikelihood <- 
      #   ggplot(likelihood_mean_sd_df, aes(x=`Simulated percent`, y=Mean)) + 
      #   geom_bar(stat="summary", aes(group=`Simulated percent`), 
      #            width=0.4, fill="gray50", color="black", position = "dodge") + 
      #   geom_errorbar(aes(ymin=Mean, ymax=Mean + Sd, group=`Simulated percent`),
      #                 size=1, width=0.3) +
      #   theme_light() +     
      #   base_plot_theme +
      #   xlab("Ground Truth Fraction") +
      #   ylab("Normalized likelihood") + 
      #   ylim(0,1.05) + 
      #   ggtitle(subfield_name) + 
      #   theme(plot.title=element_text(size=11),
      #         text=element_text(size=15))
        
      lower_limit <- min(likelihood_mean_sd_df$Mean - likelihood_mean_sd_df$Sd) * 1.1
      if(is.na(lower_limit) || lower_limit  >= -0.1) {
        lower_limit <- -0.1
      }
      
      glikelihood <-  ggplot() +
        geom_point(data=likelihood_mean_sd_df, aes(x=`Simulated percent`, y=Mean),
                   color="#132046", stroke=0)  +
        geom_line(data=likelihood_mean_sd_df, aes(x=1:nrow(likelihood_mean_sd_df), y=Mean),
                  size=.3,  color="#132046")  +
        geom_errorbar(data=likelihood_mean_sd_df, aes(x=`Simulated percent`,
                                                      ymin=Mean - Sd, 
                                                      ymax=Mean + Sd, 
                                                      group=`Simulated percent`),
                      color="#132046",
                      width=.65) +
                          base_plot_theme +
                          xlab("Ground Truth Fraction") +
                          ylab("Normalized likelihood") +
                          ylim(c(lower_limit,1.05)) + 
                          ggtitle(subfield_name) +
                          theme(plot.title=element_text(size=11),
                                text=element_text(size=15))
      likdf_normed <- likelihood_df[,1:6]
      glikelihood_wl <- glikelihood
      
      
      for (idx in 1:nrow(likdf_normed)) {
        tmp_df <- data.frame(Likelihood=unlist(likdf_normed[idx,]),
                             Sim=1:ncol(likdf_normed))
        glikelihood_wl <- 
        glikelihood_wl + 
          geom_line(data=tmp_df, aes(x=Sim, y=Likelihood), col=adjustcolor("#132046", alpha=.1))
      }
      for (size_name in names(sizes)) {
        size = sizes[[size_name]]
    
        dir.create(sprintf("%s\\%s\\%s\\",
                           write_path,
                           size_name,
                           subfield_name))
        
        pdf(file=sprintf("%s\\%s\\%s\\%s.pdf",
                         write_path,
                         size_name,
                         subfield_name,
                         conf_name),
            height=size,
            width=size)
        
        if (size_name == "small") {
          thm <- theme(axis.text.x = element_text(angle = 45, vjust=0.55))
          plot(glikelihood + thm)
        } else {
          plot(glikelihood)
        }
        
        
        dev.off()
        
        pdf(file=sprintf("%s\\%s\\%s\\with_lines_%s.pdf",
                         write_path,
                         size_name,
                         subfield_name,
                         conf_name),
            height=size,
            width=size)
        
        if (size_name == "small") {
          thm <- theme(axis.text.x = element_text(angle = 45, vjust=0.55))
          plot(glikelihood_wl + thm)
        } else {
          plot(glikelihood_wl)
        }
        
        
        dev.off()
      }
    }
  }
}

figure_4_scatter_plots <- function() {
  write_path <- sprintf("%s\\figure_4\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_4\\scatter_plots_and_mse\\",figures_path)
  dir.create(write_path)
  dir.create(sprintf("%s\\figure_4\\scatter_plots_and_mse\\big\\",figures_path))
  dir.create(sprintf("%s\\figure_4\\scatter_plots_and_mse\\small\\",figures_path))
  dir.create(sprintf("%s\\figure_4\\scatter_plots_and_mse\\medium\\",figures_path))
  dir.create(sprintf("%s\\figure_4\\scatter_plots_and_mse\\medium_1\\",figures_path))
  
  sizes=c(big=3,
          medium=2.5,
          medium_1=2.2,
          small=2)

  scatter_sizes=c(big=2,
          medium=1.5,
          medium_1=1.35,
          small=1.25)  
  
  
  true_data_paths <- sprintf("%s\\%s", all_data_paths[1:18], "equalized")
  
  sessions_to_use=c(5:8,13:16)
  
  all_df <- get_all_df(true_data_paths, sessions_to_use)
  
  residual_statistics_df <- data.frame()
  residual_values_df <- data.frame()
  estimation_configurations <- list(KS_cyclic=list(path="KS_simulations_likelihood_c", edge_free=F),
                                    KS_random=list(path="KS_simulations_likelihood_r", edge_free=F),
                                    JSD_cyclic=list(path="JSD_simulations_likelihood_c", edge_free=F),
                                    JSD_random=list(path="JSD_simulations_likelihood_r", edge_free=F),
                                    JSD_cyclic_edge_free=list(path="JSD_simulations_likelihood_c", edge_free=T),
                                    JSD_random_edge_free=list(path="JSD_simulations_likelihood_r", edge_free=T))
  
  
  for (conf_name in names(estimation_configurations)) {
    
    conf = estimation_configurations[[conf_name]]
    estimated_res <- get_estimated_df(all_df, true_data_paths, sessions_to_use,
                                      estimated_df_folder = conf$path,
                                      edge_free=conf$edge_free)
    
    likelihood_df_all <- estimated_res$likelihood_df
    
    metadata_df <- do.call(rbind, estimated_res$session_metavar_list)
    
    ind_list <- list(CA1=which(metadata_df[,4] == "CA1"),
                     CA3=which(metadata_df[,4] == "CA3"))
    
    estimated_per_percent_all <- estimated_res$estimated_per_percent_df
    estimated_per_percent_all <- cbind(estimated_per_percent_all, estimated_res$measured_vec)
    estimated_per_percent_all <- cbind(estimated_per_percent_all, estimated_res$session_df[,1])
    estimated_per_percent_all <- as.data.frame(estimated_per_percent_all)
    colnames(estimated_per_percent_all) <- c(as.character(seq(0.5, 1, by=0.1)), "measured", "session")
    
    
    
    for (subfield_name in names(ind_list)) {
      ind <- ind_list[[subfield_name]]
      estimated_per_percent_df <- estimated_per_percent_all[ind,]
      
      squared_residuals_df <- c()
      squared_residuals_df_for_test  <- c()
      scatter_list <- list()
      
      for (est_p in as.character(seq(0.5, 1, by=0.1))) {
        
        tmp_df <- data.frame(Measured=as.numeric(estimated_per_percent_df[,"measured"]),
                             Estimated=as.numeric(estimated_per_percent_df[,est_p]),
                             Session=as.numeric(estimated_per_percent_df[,"session"]),
                             Mice=metadata_df[ind,2])
        # 
         tmp_df$Session[tmp_df$Session > 8] <- tmp_df$Session[tmp_df$Session > 8] - 8
         tmp_df$Session <- factor(sprintf("%dth", tmp_df$Session), c("5th", "6th", "7th", "8th"))
        
        linear_line_df <- data.frame(x=c(0.3,1),
                                     y=c(0.3,1))
        gs <- 
          ggplot(linear_line_df) +
          geom_line(aes(x=x,y=y), linetype="dashed", color="#c0c1c3",
                    size=1, alpha=0.5) +    
          geom_point(data=tmp_df, aes(x=Measured, y=Estimated, color=Session), size=0.8) +
          theme_light() +     
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.border = element_blank(),
                panel.background = element_blank(),
                legend.position="NA",
                text=element_text(size=15),
                plot.margin=unit(c(0,0,0,0), "cm")) +
          scale_color_manual(breaks=c("5th","6th","7th","8th"),
                             values=c(rev(brewer.pal(4, "RdYlBu")))) + 
          ylim(c(0.3,1)) + 
          xlim(c(0.3,1)) + 
          ylab("Frac") +
          xlab("Frac")
        
        squared_residuals <- (tmp_df$Estimated - tmp_df$Measured) ** 2
        #squared_residuals <- ddply(tmp_df, .(Session), function(mice_df) {mean((mice_df$Estimated - mice_df$Measured) ** 2)})[,2]
        squared_residuals_df <- rbind(squared_residuals_df,
                                      data.frame(residual=squared_residuals,
                                                 gt_frac=rep(est_p, len(squared_residuals))))
        squared_residuals_df_for_test <-
                                  cbind(squared_residuals_df_for_test,
                                        squared_residuals)
        scatter_list <- append(scatter_list, list(gs))
      }
      
      scatter_list$nrow = 2
      scattersp <- do.call(arrangeGrob, scatter_list)
      
      
    
      
      for (size_name in names(scatter_sizes)) {
        size = scatter_sizes[[size_name]]
        dir.create(sprintf("%s\\%s\\%s",write_path, size_name, subfield_name))
        dir.create(sprintf("%s\\%s\\%s\\scatter_plots",
                           write_path,
                           size_name,
                           subfield_name))
        
        pdf(file=sprintf("%s\\%s\\%s\\scatter_plots\\%s.pdf",
                              write_path,
                              size_name,
                              subfield_name,
                              conf_name),
            height=size * 2,
            width=size * 3.2)
        
        plot(scattersp)
        dev.off()
      }
      
      sd_df <- ddply(squared_residuals_df, .(gt_frac), function(frac_df) {c(mean(frac_df[,"residual"]),
                                                                            mean(frac_df[,"residual"]) -
                                                                            sd(frac_df[,"residual"]),
                                                                            mean(frac_df[,"residual"]) +
                                                                            sd(frac_df[,"residual"]),
                                                                            sd(frac_df[,"residual"]),
                                                                            sem(frac_df[,"residual"]))})
      colnames(sd_df) <- c("gt_frac", "mse", "lower", "upper", "sd", "sem")
      
      sd_df$conf = conf_name
      sd_df$subfield = subfield_name
      sd_df$N <- nrow(squared_residuals_df) / nrow(sd_df)
      residual_values_df <- rbind(residual_values_df,
                                  sd_df)
      g_residual_box <- 
      ggplot(squared_residuals_df,
                   aes(x=gt_frac, y=residual)) +
        # geom_linerange(data=sd_df, aes(x=gt_frac, ymin=lower, ymax=upper, y=mse),
        #              size=1,
        #              color="gray65") + 
        # geom_point(data=sd_df, aes(x=gt_frac, y=mse)) +
        geom_boxplot(color="gray65", fill=adjustcolor("gray80", alpha=0.8), size=1, width=0.5) +
        geom_jitter(color="gray20", 
                    size=0.75, 
                    alpha=0.1) +
        ylab("MSE") +
        xlab("Ground Truth Fraction") +
        base_plot_theme + 
        theme(plot.title=element_text(size=11),
              text=element_text(size=15)) +
        ylim(c(0, max(squared_residuals_df_for_test) * 1.3)) + 
        ggtitle(subfield_name)
      
      
      
      pos_y <- 
      c(max(squared_residuals_df_for_test)  * 0.4,
        max(squared_residuals_df_for_test)  * 0.6,
        max(squared_residuals_df_for_test)  * 0.8,
        max(squared_residuals_df_for_test)  * 1,
        max(squared_residuals_df_for_test)  * 1.2)
      
      pos_y <- rev(pos_y)
      
      #min_mse <- which.min(colMeans(squared_residuals_df_for_test))
      min_mse <- which.max(colMeans(likelihood_df_all[ind,]))
      
      comp_idx = 5
      
      g_signif <- g_residual_box
      for (comp in  c(6:1)) {
        
        if (comp == min_mse) {
          next
        }
        
        pv <- wilcox.test(squared_residuals_df_for_test[,min_mse], 
                          squared_residuals_df_for_test[,comp],
                          alternative="less",
                          correct=F,
                          paired=T)
        x_coords <- c(min_mse,min_mse,comp,comp)
        
        if (comp > min_mse) {
          x_coords <- rev(x_coords)
          names(x_coords) <- c()
        }
        
        cords_df <- data.frame(x=x_coords,
                               y=c(pos_y[comp_idx] * 1.05, pos_y[comp_idx], 
                                   pos_y[comp_idx], pos_y[comp_idx] * 1.05))
        
        corrected_pv <- pv$p.value * 5
        if (corrected_pv >= 1) {
            corrected_pv = 1
        }
        
        sig_label <- signif.num(corrected_pv)
        
        residual_statistics_df <- rbind(residual_statistics_df,
                                        data.frame(W=pv$statistic,
                                                   pvalue=pv$p.value,
                                                   alternative=pv$alternative,
                                                   subfield=subfield_name,
                                                   comparision=sprintf("%d - %d", min_mse, comp),
                                                   conf=conf_name,
                                                   corrected_pv = corrected_pv,
                                                   code=sig_label,
                                                   method=pv$method))
        
        
        if(sig_label == "n.s.") {
          jump_label = 1.15
        } else {
          jump_label = 1.05
        }
        
        g_signif <- 
        g_signif + geom_line(data=cords_df,
                             aes(x=x,y=y)) +
                   geom_text(x=mean(cords_df$x),
                             y=max(cords_df$y) * jump_label,
                             label=sig_label)
        comp_idx <- comp_idx - 1
        
      }
      
      for (size_name in names(sizes)) {
        size = sizes[[size_name]]
        dir.create(sprintf("%s\\%s\\%s",write_path, size_name, subfield_name))
        dir.create(sprintf("%s\\%s\\%s\\mse_plots",
                           write_path,
                           size_name,
                           subfield_name))
        
        pdf(file=sprintf("%s\\%s\\%s\\mse_plots\\%s.pdf",
                         write_path,
                         size_name,
                         subfield_name,
                         conf_name),
            height=size,
            width=size)
        
        plot(g_residual_box)
        dev.off()
        
        pdf(file=sprintf("%s\\%s\\%s\\mse_plots\\signif_codes_%s.pdf",
                         write_path,
                         size_name,
                         subfield_name,
                         conf_name),
            height=size,
            width=size)
        
        plot(g_signif)
        dev.off()
      }
    }
  }
  
  
  write.csv(residual_statistics_df,
            file=sprintf("%s\\residuals_statistics.csv",write_path), quote=T)
  
  write.csv(residual_values_df,
            file=sprintf("%s\\residuals_values.csv",write_path), quote=T)
}

figure_4_activity_plots <- function() {
  write_path <- sprintf("%s\\figure_4\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_4\\activity_plots\\",figures_path)
  dir.create(write_path)
  dir.create(sprintf("%s\\figure_4\\activity_plots\\big\\",figures_path))
  dir.create(sprintf("%s\\figure_4\\activity_plots\\small\\",figures_path))
  dir.create(sprintf("%s\\figure_4\\activity_plots\\medium\\",figures_path))
  dir.create(sprintf("%s\\figure_4\\activity_plots\\medium_1\\",figures_path))
  
  sizes=c(big=3,
          medium=2.5,
          medium_1=2.2,
          small=2)
  
  scatter_sizes=c(big=2,
                  medium=1.5,
                  small=1.25)  
  
  
  true_data_paths <- sprintf("%s\\%s", all_data_paths[1:18], "equalized")
  
  sessions_to_use=c(5:8,13:16)
  
  act_df <- 
  create_pct_by_activity_df(true_data_paths,
                            indices = sessions_to_use)
  
  shuffle_dfs <- c(cyclic="cyclic_df", 
                   random="random_df")
  
  for (shuffle_type_name in names(shuffle_dfs)) {
    shuffle_type <- shuffle_dfs[[shuffle_type_name]]
    
    trace_colors = c(CA1="#CE3736", 
                     CA3="#3F5DAB")
      
      for (subfield in c("CA1", "CA3")) {
        
        indices <- which(act_df$subfields == subfield)
        
        df <- act_df[[shuffle_type]][,indices]
        
        fraction_sd <- apply(df, 1, sd, na.rm=T)
        fraction_mean <- rowMeans(df, na.rm=T) 
        bins <- as.numeric(names(fraction_sd))
        df <- df[order(bins),]
        fraction_sd <- fraction_sd[order(bins)]
        fraction_mean <- fraction_mean[order(bins)]
        bins <- bins[order(bins)]
        
        mean_sd_df <- data.frame(y=fraction_mean,
                                 sd=fraction_sd,
                                 x=bins)
        g <- ggplot(mean_sd_df, aes(x=x, y=y)) + base_plot_theme
        
        
        gcor2 <- 
          ggplot(act_df[1:11,],aes(x=as.numeric(bins),y=pct)) + 
          geom_errorbar(aes(ymin=pct-sd, ymax=pct+sd), size=1, width=2) +
          geom_line(size=1.2) +
          geom_point(size=3) +
          theme_light() +
         # ylim(c(0.25, 0.7))+
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.position="NA",
                panel.border = element_blank(),
                panel.background = element_blank()) +
          geom_hline(yintercept=1, linetype="dashed") +
          labs(x="Ca2+ events", 
               y = "Fraction of place cells (%)")
        
        
        for (i in 1:ncol(df)) {
          tmp_df <- data.frame(x=bins,
                               y=df[,i])
          
          g <- g +
            geom_line(data=tmp_df, aes(x=x,y=y),
                      color=adjustcolor(trace_colors[[subfield]],
                                        alpha=0.05),
                      size=0.7)
        }
        
        gf <- 
          g + geom_hline(yintercept=1, linetype="dashed", size=0.8) + 
          geom_line(size=1,
                    linetype="dashed",
                    col=trace_colors[[subfield]]) + 
          geom_errorbar(aes(x=x, ymin=y-sd, ymax=y+sd),
                        col=trace_colors[[subfield]],
                        size=1,
                        width=5) +
          geom_point(col=trace_colors[[subfield]],
                     size=2.5) +
          ylim(c(0.2,1.1)) + 
          xlim(c(0,110)) + 
          ggtitle(subfield) + 
          ylab("Naive fraction") + 
          xlab("Firing events (#)") + 
          theme(plot.title=element_text(size=11),
                text=element_text(size=15))
        
        for (size_name in names(sizes)) {
          size <- sizes[[size_name]]
          pdf(file=sprintf("%s\\%s\\%s_%s.pdf",
                           write_path,
                           size_name,
                           shuffle_type_name,
                           subfield),
              height=size,
              width=size)
          plot(gf)
          dev.off()
        }
      }
  }
}

figure_4_diluted_pcs <- function() {
  write_path <- sprintf("%s\\figure_4\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_4\\diluted_pcs\\",figures_path)
  dir.create(write_path)
  dir.create(sprintf("%s\\figure_4\\activity_plots\\big\\",figures_path))
  dir.create(sprintf("%s\\figure_4\\activity_plots\\small\\",figures_path))
  dir.create(sprintf("%s\\figure_4\\activity_plots\\medium\\",figures_path))
  dir.create(sprintf("%s\\figure_4\\activity_plots\\medium_1\\",figures_path))
  
  sizes=c(big=3,
          medium=2.5,
          medium_1=2.2,
          small=2)
  
  
  
  true_data_paths <- sprintf("%s\\%s", all_data_paths[1:8], "equalized")
  file_name = "diluted_pcs.Rda"
  sessions_to_use <- 1:16
  results_df <- data.frame()
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
      
      if (len((grep(file_name, list.files(sprintf("%s%s", tpath, ""))))) == 0) {
        print(tpath)
        missing_paths <- c(missing_paths, tpath)
        print("Missing")
        next
      }
      
      if (verbose) {
        print(sprintf("Loading %s", sprintf("%s%s\\%s", tpath, "", file_name)))
      }
      
      load(sprintf("%s%s\\%s", tpath, "", file_name), verbose=F)
      
      
      if (grepl("CA3", working_path)) {
        mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', tpath))
        mice_str <- substr(tpath, mice_str_index, mice_str_index+4) 
      } else {
        mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', tpath))
        mice_str <- substr(tpath, mice_str_index, mice_str_index+3) 
      }
      
      subfield <- ifelse(grepl("CA3", working_path), "CA3", "CA1")
    
      
      results_df <- 
        rbind(results_df,
              data.frame(Diluted=mean(final_result$diluted_frac),
                         Top=final_result$top_20_frac,
                         Lower=final_result$lower_20_frac,
                         Mice=mice_str,
                         Session=idx))
    }
  }
  
  melted <- melt(results_df[,-4], values = c("Diluted", "Top", "Lower"), id.vars ="Session")
  gbox <- 
    ggplot(melted, aes(x=variable, y=value))+ 
    geom_boxplot() +
    base_plot_theme +
    xlab("") +
    ylab("Fraction of place cells (%)") +
    theme(plot.title=element_text(size=11),
          text=element_text(color="black",size=15),
          axis.text=element_text(color="black"),
          axis.ticks = element_line(color="black"))
  
  # for (ridx in 1:nrow(results_df)) {
  #   row_df <- 
  #   data.frame(x=1:3, 
  #              y=unlist(as.vector(results_df[ridx,1:3])))
  #   
  #   gbox <- gbox + geom_line(data=row_df, aes(x=x,y=y))
  #   
  # }
  for (size_name in names(sizes)) {
    size <- sizes[[size_name]]
    dir.create(sprintf("%s\\%s\\",write_path, size_name))
    pdf(file=sprintf("%s\\%s\\CA1_dilution.pdf",
                     write_path,
                     size_name),
        height=size,
        width=size)
    plot(gbox)
    dev.off()
  }
  
}
