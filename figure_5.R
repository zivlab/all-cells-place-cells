
figure_5_lineplots  <- function () {
  
  write_path <- sprintf("%s\\figure_5\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_5\\lineplots\\",figures_path)
  directions_to_use = c("Both")
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
  
  sessions_to_use=c(1:16)
  
  all_df <- get_all_df(true_data_paths, sessions_to_use)
  
  estimations_statistics_df <- data.frame()
  
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
    
    both_dir_res <- get_both_dir_df(estimated_df = estimated_res$estimated_df, all_df=all_df, dev_func = sem)
    
    both_dir_df <- both_dir_res$both_dir_df
    
    #### statistics 
    for (direc in unique(both_dir_df$Direction)) {
      anova_df <- both_dir_df[both_dir_df$Direction == direc,]
      colnames(anova_df) <- c("Frac", "Session", "Mice", "#", "Subfield")
      ca3_anova_df <- anova_df[anova_df$Subfield == "CA3",]
      ca1_anova_df <- anova_df[anova_df$Subfield == "CA1",]
      
      model.aov <- aov(data=anova_df,
                       formula=Frac ~ 
                         Subfield * Session + 
                         Error(Mice/(Subfield*Session)))
      
      pvals <- summary(model.aov)
      pvals_df <- as.data.frame(pvals$`Error: Mice:Session`[[1]])
      pvals_df$Group = "Two way"
      pvals_df$Measurement <- "Estimated"
      pvals_df$Direction = direc
      pvals_df$Conf = conf_name
      pvals_df$Factor <- rownames(pvals_df)
      rownames(pvals_df) <- c()
      estimations_statistics_df <- rbind(estimations_statistics_df, pvals_df)
      
      print("2way")
      print(signif.num(pvals$`Error: Mice`[[1]]$`Pr(>F)`[1:2]))
      print(signif.num(pvals$`Error: Mice:Session`[[1]]$`Pr(>F)`[1:2]))
      
      print("CA1")
      model.aov <- aov(data=ca1_anova_df, formula=Frac~Session + Error(factor(Mice)))
      pvals <- summary(model.aov)
      print(signif.num(pvals$`Error: Within`[[1]]$`Pr(>F)`[1]))
      
      pvals_df <- as.data.frame(pvals$`Error: Within`[[1]])
      pvals_df$Group <- "CA1"
      pvals_df$Measurement <- "Estimated"
      pvals_df$Conf = conf_name
      pvals_df$Direction = direc
      pvals_df$Factor <- rownames(pvals_df)
      rownames(pvals_df) <- c()
      estimations_statistics_df <- rbind(estimations_statistics_df, pvals_df)
      
      print("CA3")
      model.aov <- aov(data=ca3_anova_df, formula=Frac~Session + Error(factor(Mice)))
      pvals <- summary(model.aov)
      print(signif.num(pvals$`Error: Within`[[1]]$`Pr(>F)`[1]))
      
      pvals_df <- as.data.frame(pvals$`Error: Within`[[1]])
      pvals_df$Group <- "CA3"
      pvals_df$Measurement <- "Estimated"
      pvals_df$Conf = conf_name
      pvals_df$Direction = direc
      pvals_df$Factor <- rownames(pvals_df)
      rownames(pvals_df) <- c()
      estimations_statistics_df <- rbind(estimations_statistics_df, pvals_df)
    }
    ####
    
    mean_sd_df <- both_dir_res$mean_sd_df
    session_mean_sd_df <- both_dir_res$session_mean_sd_df
    
    mean_sd_df <- mean_sd_df[mean_sd_df$Direction %in% directions_to_use, ]
    both_dir_df <- both_dir_df[both_dir_df$Direction %in% directions_to_use, ]
    session_mean_sd_df <- session_mean_sd_df[session_mean_sd_df$Direction %in% directions_to_use, ]
    
    g_est <- 
      ggplot(session_mean_sd_df, aes(x=Session, y=Estimated)) + 
      # geom_point(data=both_dir_df, aes(x=Session, y=Estimated, group=Subfield), 
      #            fill="gray75", color="gray75",
      #            size=.5, alpha=0.2, position=position_dodge(0.4)) +
      geom_line(aes(group=Subfield, color=Subfield), size=1, position=position_dodge(0.4)) +
      
      geom_errorbar(aes(ymin=Estimated - Sd, ymax=Estimated + Sd, group=Subfield, color=Subfield), 
                    size=0.75, width=1, position=position_dodge(0.4)) +
      geom_point(aes(group=Subfield, color=Subfield), size=2, position=position_dodge(0.4)) + 
      theme_light() +     
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position="NA",
            text=element_text(size=15),
            plot.title=element_text(size=11)) + 
      scale_color_manual(values=c(CA1="#CE3736", 
                                  CA3="#3F5DAB")) +
      xlab("Session") +
      ylab("Estimated fraction") + 
      ylim(0,1.1) + 
      geom_hline(yintercept=1, linetype="dashed") +
      ggtitle("Place cells (estimated)")
    
    
    
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
  
  
  estimations_statistics_df$signif_code <- signif.num(estimations_statistics_df$`Pr(>F)`)
  
  write.csv(estimations_statistics_df,
            file=sprintf("%s\\estimated_measurements_statistics.csv",write_path))
  
  naive_statistics_df <- data.frame()
  conf = estimation_configurations[[1]]
  
  naive_configurations <- list(cyclic_eq_prior_weighted=list(path="percent_df_cyclic_SI_eq_prior.R", weighted=T),
                               cyclic_eq_prior_either_dir=list(path="percent_df_cyclic_SI_eq_prior.R", weighted=F),
                               random_eq_prior_weighted=list(path="percent_df_rand_SI_eq_prior.R", weighted=T),
                               random_eq_prior_either_dir=list(path="percent_df_rand_SI_eq_prior.R", weighted=F),
                               cyclic_non_eq_prior_weighted=list(path="percent_df_cyclic_SI_reg.R", weighted=T),
                               cyclic_non_eq_prior_either_dir=list(path="percent_df_cyclic_SI_reg.R", weighted=F),
                               random_non_eq_prior_weighted=list(path="percent_df_rand_SI_reg.R", weighted=T),
                               random_non_eq_prior_either_dir=list(path="percent_df_rand_SI_reg.R", weighted=F))
  
  for (naive_conf_name in names(naive_configurations)) { 
    naive_conf = naive_configurations[[naive_conf_name]]
    
    all_df <- get_all_df(true_data_paths, 
                         sessions_to_use,
                         file_name = naive_conf$path)
    
    estimated_res <- get_estimated_df(all_df, true_data_paths, sessions_to_use)
    
    both_dir_res <- get_both_dir_df(estimated_df = estimated_res$estimated_df, all_df=all_df, dev_func=sem)
    
    both_dir_df <- both_dir_res$both_dir_df
    mean_sd_df <- both_dir_res$mean_sd_df
    session_mean_sd_df <- both_dir_res$session_mean_sd_df
    
    mean_sd_df <- mean_sd_df[mean_sd_df$Direction %in% directions_to_use, ]
    both_dir_df <- both_dir_df[both_dir_df$Direction %in% directions_to_use, ]
    
    both_dir_df_all <- both_dir_res$both_dir_df_all
    session_mean_sd_df_all <- both_dir_res$session_mean_sd_df_al
    
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
    
    ylabs <- c("Fraction",
               "Naive fraction",
               "Naive fraction")
    
    
    for (i  in 1:len(dfs)) {
      df_idx <- i
      pn <- weighted_avg_both_df_colnames[i]
      
      work_df <- dfs[[i]]
      
      ### Statistics 
      for(direc in unique(work_df$variable)) {
        
        
        anova_df <- work_df[work_df$variable == direc,]
        anova_df$SubID <- 1:nrow(anova_df)
        anova_df$Session <- as.character(as.numeric(anova_df$Session) %% 8)
        anova_df$Session[anova_df$Session == "0"] <- "8"
        colnames(anova_df) <- c("Session","Mice", "Subfield", "#", "Frac", "SubID")
        
        
        ca3_anova_df <- anova_df[anova_df$Subfield == "CA3",]
        ca1_anova_df <- anova_df[anova_df$Subfield == "CA1",]
        
        model.aov <- aov(data=anova_df,
                         formula=Frac ~ 
                           Subfield * Session + 
                           Error(Mice/(Subfield*Session)))
        
        pvals <- summary(model.aov)
        pvals_df <- as.data.frame(pvals$`Error: Mice:Session`[[1]])
        pvals_df$Group = "Two way"
        pvals_df$Measurement <- plot_titles[i]
        pvals_df$Conf = naive_conf_name
        pvals_df$Direction = direc
        pvals_df$Factor <- rownames(pvals_df)
        rownames(pvals_df) <- c()
        naive_statistics_df <- rbind(naive_statistics_df, pvals_df)
        
        print("2way")
        print(signif.num(pvals$`Error: Mice`[[1]]$`Pr(>F)`[1:2]))
        print(signif.num(pvals$`Error: Mice:Session`[[1]]$`Pr(>F)`[1:2]))
        
        print("CA1")
        model.aov <- aov(data=ca1_anova_df, formula=Frac~Session + Error(factor(Mice)))
        pvals <- summary(model.aov)
        print(signif.num(pvals$`Error: Within`[[1]]$`Pr(>F)`[1]))
        
        pvals_df <- as.data.frame(pvals$`Error: Within`[[1]])
        pvals_df$Group <- "CA1"
        pvals_df$Measurement <- plot_titles[i]
        pvals_df$Conf = naive_conf_name
        pvals_df$Direction = direc
        pvals_df$Factor <- rownames(pvals_df)
        rownames(pvals_df) <- c()
        naive_statistics_df <- rbind(naive_statistics_df, pvals_df)
        
        print("CA3")
        model.aov <- aov(data=ca3_anova_df, formula=Frac~Session + Error(factor(Mice)))
        pvals <- summary(model.aov)
        print(signif.num(pvals$`Error: Within`[[1]]$`Pr(>F)`[1]))
        
        pvals_df <- as.data.frame(pvals$`Error: Within`[[1]])
        pvals_df$Group <- "CA3"
        pvals_df$Measurement <- plot_titles[i]
        pvals_df$Conf = naive_conf_name
        pvals_df$Direction = direc
        pvals_df$Factor <- rownames(pvals_df)
        rownames(pvals_df) <- c()
        naive_statistics_df <- rbind(naive_statistics_df, pvals_df)
      }
      
      work_df <- work_df[work_df$variable %in% directions_to_use,]
      work_df$Session <- as.character(as.numeric(work_df$Session) %% 8)
      work_df$Session[work_df$Session == "0"] <- "8"
      
      if (naive_conf$weighted) {
        session_mean_work_df <- session_mean_sd_df_all[,c("Subfield", "Session", "Direction", sprintf("%s_%s", pn, c("mean", "sd")))]
        colnames(session_mean_work_df) <- c("Subfield", "Session", "Direction", "Mean", "Sd")
        
      } else {
        
        
        session_mean_work_df <- 
          ddply(work_df, .(Subfield), function(subfield_df) {
            both_df  <- ddply(subfield_df, .(Session), function(session_df){
              
              session_df
              return(c(mean(session_df[,"value"]), 
                       sd(session_df[,"value"])))
            })
          })
        colnames(session_mean_work_df) <- c("Subfield", "Session", "Mean", "Sd")
      }
      
      
      
      g <- 
        ggplot(session_mean_work_df, aes(x=Session, y=Mean)) + 
        geom_line(aes(group=Subfield, color=Subfield), size=1, position=position_dodge(0.4)) +
        geom_errorbar(aes(ymin=Mean - Sd, ymax=Mean + Sd, group=Subfield, color=Subfield), 
                      size=.75, width=1, position=position_dodge(0.4)) +
        geom_point(aes(group=Subfield, color=Subfield), size=2, position=position_dodge(0.4)) + 
        theme_light() +     
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="NA",
              text=element_text(size=15),
              plot.title=element_text(size=11))   + 
        scale_color_manual(values=c(CA1="#CE3736", 
                                    CA3="#3F5DAB")) +
        xlab("Session") +
        ylab(ylabs[i]) + 
        ylim(0,1.1) + 
        ggtitle(plot_titles[i])
      
      
      
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
  
  naive_statistics_df$signif_code <- signif.num(naive_statistics_df$`Pr(>F)`)
  
  write.csv(naive_statistics_df,
            file=sprintf("%s\\naive_measurements_statistics.csv",write_path))
}

figure_5_likelihood_plots <- function() {
  
  write_path <- sprintf("%s\\figure_5\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_5\\likelihood_plots_new\\",figures_path)
  dir.create(write_path)
  dir.create(sprintf("%s\\figure_5\\likelihood_plots_new\\big\\",figures_path))
  dir.create(sprintf("%s\\figure_5\\likelihood_plots_new\\small\\",figures_path))
  dir.create(sprintf("%s\\figure_5\\likelihood_plots_new\\medium\\",figures_path))
  dir.create(sprintf("%s\\figure_5\\likelihood_plots_new\\medium_1\\",figures_path))
  dir.create(sprintf("%s\\figure_5\\likelihood_plots_new\\medium_2\\",figures_path))
  
  sizes=c(big=2.25,
          medium=2,
          medium_1=1.75,
          medium_2=1.5,
          small=1.25)
  
  dir.create(write_path)
  
  true_data_paths <- sprintf("%s\\%s", all_data_paths[1:18], "equalized")
  
  sessions_to_use=c(1:16)
  
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
      #likelihood_sessions <- ceiling(likelihood_sessions / 2)
      likelihood_df <- cbind(likelihood_df, likelihood_sessions)
      colnames(likelihood_df) <- c(as.character(seq(0.5,1,by=0.1)), "sessions")
      melted_likelihood <- melt(likelihood_df, id.vars="sessions")
      colnames(melted_likelihood) <- c("Session", "Simulated percent", "Normalized likelihood")
      
      
      
      gt_max_likelihood_by_session_df <- 
          ddply(melted_likelihood, 
                 .(Session), 
                 function(df) {
                   idf <- ddply(df, 
                                .(`Simulated percent`), 
                                function(idf) { mean(idf[,"Normalized likelihood"])}); 
                   return(which.max(idf[,"V1"]))
                   })
      
      
      colnames(gt_max_likelihood_by_session_df) <- c("Session", "GT")
      
      gt_max_likelihood_by_session_df$GT <- seq(0.5,1,length.out=6)[gt_max_likelihood_by_session_df$GT]
      
      
      gmax_likelihood_by_session <- 
        ggplot(gt_max_likelihood_by_session_df,
               aes(x=Session, y=GT)) +
        geom_point(stat="summary") + 
        geom_line(stat="summary") + 
          base_plot_theme + 
          theme(legend.position="top",
                axis.text=element_text(color="black"))
      
      gsessions_all <- 
        ggplot(melted_likelihood, 
               aes(x=`Simulated percent`, 
                   y=`Normalized likelihood`, 
                   color=Session, 
                   group=Session)) + 
        geom_point(stat="summary") + 
        geom_line(stat="summary") + 
        scale_color_distiller(palette="RdYlBu") +
        base_plot_theme +
        theme(legend.position="top",
              axis.text=element_text(color="black"))
      session_combinations <- list(all=list(sessions=c(1:8),width=4),
                                   first_last=list(sessions=c(1,2),width=2),
                                   all_all=list(sessions=1:8, widths=4))
      
      for (combination_name in names(session_combinations)) {
        
        sessions_to_use_for_plot <- session_combinations[[combination_name]]$sessions
        ses_plots <- list()
        
        
        melted_likelihood$Session
        
        
        for (ses in sessions_to_use_for_plot) {
          melted_ses_likelihood <- melted_likelihood[melted_likelihood$Session == ses,]
          likelihood_mean_sd_df <- ddply(melted_ses_likelihood, .(`Simulated percent`), 
                                         function(sub_df) {
                                           return(c(mean(sub_df[,"Normalized likelihood"]),
                                                    sd(sub_df[,"Normalized likelihood"])))
                                         })
          
          colnames(likelihood_mean_sd_df) <- c("Simulated percent", "Mean", "Sd")
          
          lower_limit <- min(likelihood_mean_sd_df$Mean - likelihood_mean_sd_df$Sd) * 1.1
          if(is.na(lower_limit) || lower_limit  >= -0.1) {
            lower_limit <- -0.1
          }
          
          likelihood_mean_sd_df$`Simulated percent` <- str_replace(as.character(likelihood_mean_sd_df$`Simulated percent`), "0", "")
          
          
          glikelihood <-  ggplot() +
            geom_point(data=likelihood_mean_sd_df, aes(x=`Simulated percent`, y=Mean),
                       color="black", stroke=0)  +

            geom_line(data=likelihood_mean_sd_df, aes(x=1:nrow(likelihood_mean_sd_df), y=Mean),
                      size=.3, linetype="dashed")  +

            geom_ribbon(data=likelihood_mean_sd_df, aes(x=`Simulated percent`,
                                                          ymin=Mean - Sd,
                                                          ymax=Mean + Sd,
                                                          group=`Simulated percent`),
                          width=.5) +
            theme_light() +
            base_plot_theme +
            xlab("") +
            ylab("") +
            ylim(c(lower_limit,1.05)) +
            theme(plot.title=element_text(size=11),
                  text=element_text(size=13.5),
                  axis.text.x = element_text(size=10.5, angle = 40, vjust=.5),
                  plot.margin=unit(c(0,0,0,0), "cm"))
          
          
          
          # Swapped for overlay - note for self
          glikelihood <- 
          ggplot(likelihood_mean_sd_df, aes(x=`Simulated percent`,y=Mean)) + 
            geom_point() + 
            geom_line(aes(group=1)) + 
            geom_ribbon(aes(group=1, x=`Simulated percent`, ymin=Mean-Sd, ymax=Mean+Sd),
                        color=NA, alpha=.2) +
            theme_light() +     
            base_plot_theme +
            xlab("") +
            ylab("") + 
            ylim(c(lower_limit,1.05)) + 
            theme(plot.title=element_text(size=11),
                  text=element_text(size=13.5),
                  axis.text.x = element_text(size=10.5, angle = 40, vjust=.5),
                  plot.margin=unit(c(0,0,0,0), "cm")) 
          
          
          
          #ggtitle(sprintf("Session = %d", ses))
          
          ses_plots <- append(ses_plots,
                              list(glikelihood))
        }
        
        
        
        ses_plots$nrow <- 4
        ses_plots$align <- "v"
        gf <- do.call(plot_grid, ses_plots)
        
        
        
        
        for (size_name in names(sizes)) {
          size = sizes[[size_name]]
          
          
          dir.create(sprintf("%s\\%s\\%s\\",
                             write_path,
                             size_name,
                             combination_name))
          
          dir.create(sprintf("%s\\%s\\%s\\%s\\",
                             write_path,
                             size_name,
                             combination_name,
                             subfield_name))
          
          pdf(file=sprintf("%s\\%s\\%s\\%s\\ribbon_%s.pdf",
                           write_path,
                           size_name,
                           combination_name,
                           subfield_name,
                           conf_name),
              width=size * 1,
              height=size * session_combinations[[combination_name]]$width)
          
          plot(gf)
          dev.off()
          
          pdf(file=sprintf("%s\\%s\\%s\\%s\\all_overlaid_%s.pdf",
                           write_path,
                           size_name,
                           combination_name,
                           subfield_name,
                           conf_name),
              width=size * 1,
              height=size * 1)
          
          plot(gsessions_all)
          dev.off()
          
          
        
          
          pdf(file=sprintf("%s\\%s\\%s\\%s\\max_likelihood_by_session%s.pdf",
                           write_path,
                           size_name,
                           combination_name,
                           subfield_name,
                           conf_name),
              width=size * 1,
              height=size * 1)
          
          plot(gmax_likelihood_by_session)
          
          dev.off()
          
          
        }
      }
    }
  }
}
