supp_figure_different_cost_functions_boxplots  <- function () {
  
  directions_to_use = c("Left", "Right")
  write_path <- sprintf("%s\\supp_figure_different_cost_functions\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\supp_figure_different_cost_functions\\boxplots\\",figures_path)
  
  sizes=c(big=3,
          medium=2.5,
          medium_1=2.2,
          small=2)
  
  dir.create(write_path)
  dir.create(sprintf("%s", write_path))
  dir.create(sprintf("%s\\big\\", write_path))
  dir.create(sprintf("%s\\medium\\", write_path))
  dir.create(sprintf("%s\\medium_1\\", write_path))
  dir.create(sprintf("%s\\small\\", write_path))
  
  true_data_paths <- sprintf("%s\\%s", all_data_paths[1:18], "equalized")
  
  sessions_to_use=c(5:8,13:16)
  
  
  
  estimation_configurations <- list(KS_cyclic=list(path="KS_simulations_likelihood_c", edge_free=F),
                                    KS_random=list(path="KS_simulations_likelihood_r", edge_free=F),
                                    JSD_cyclic=list(path="JSD_simulations_likelihood_c", edge_free=F),
                                    JSD_random=list(path="JSD_simulations_likelihood_r", edge_free=F),
                                    JSD_cyclic_edge_free=list(path="JSD_simulations_likelihood_c", edge_free=T),
                                    JSD_random_edge_free=list(path="JSD_simulations_likelihood_r", edge_free=T))
  
  values_df <- data.frame()
  statistics_df <- data.frame()
  directions_to_use <- c("Both", "Left", "Right")
  
  subfield_indices <- list("CA1"=1:8,"CA3"=9:18)
  for (subfield_name in names(subfield_indices)) {
    indices = subfield_indices[[subfield_name]]
    all_df <- get_all_df(true_data_paths[indices], sessions_to_use)
    
    for (conf_name in names(estimation_configurations)) {
      
      conf = estimation_configurations[[conf_name]]
      estimated_res <- get_estimated_df(all_df, true_data_paths[indices], sessions_to_use,
                                        estimated_df_folder = conf$path,
                                        edge_free=conf$edge_free)
      
      both_dir_res <- get_both_dir_df(estimated_df = estimated_res$estimated_df, all_df=all_df, dev_func = sem)
      
      both_dir_df <- both_dir_res$both_dir_df
      mean_sd_df <- both_dir_res$mean_sd_df
      session_mean_sd_df <- both_dir_res$session_mean_sd_df
      
      mean_sd_df <- mean_sd_df[mean_sd_df$Direction %in% directions_to_use, ]
      both_dir_df <- both_dir_df[both_dir_df$Direction %in% directions_to_use, ]
      mean_sd_df$Conf = conf_name
      values_df <- rbind(values_df, mean_sd_df)
      
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
                     width=0.5,
                     position = position_dodge(0.75),
                     size=1)
      
      g_est <- color_plot(g_est, x_axes) +
        theme(plot.title=element_text(size=11))
      
      
      #### Statistics 
      directions <- unique(both_dir_df$Direction)
      direction_x_indices <- c("Both"=1, "Left"=2, "Right"=3)
      combination_matrix <- combn(len(directions), 2)
      gfinal <- g_est 
      
      for (correct_continuouity in c(T, F)) {
        for (comparision in 1:ncol(combination_matrix)) {
          
          direc1 <- directions[combination_matrix[1,comparision]]
          direc2 <- directions[combination_matrix[2,comparision]]
          
          print(sprintf("Comparing %s - %s",
                        direc1,
                        direc2))
          
          df1 <- both_dir_df[both_dir_df$Direction == direc1,]
          df1 <- df1[order(paste(df1$TwoEnvSession, df1$Mice, sep="_")),]
          
          df2 <- both_dir_df[both_dir_df$Direction == direc2,]
          df2 <- df2[order(paste(df2$TwoEnvSession, df2$Mice, sep="_")),]
          
          print(all(df1$TwoEnvSession  ==  df2$TwoEnvSession))
          print(all(df1$Mice  ==  df2$Mice))
          
          wilk <- wilcox.test(df1$Estimated, df2$Estimated, paired=T, correct = correct_continuouity)
          wilk_correction <- wilk$p.value * ncol(combination_matrix)
          wilk_correction <- ifelse(wilk_correction > 1, 1, wilk_correction)
          
          d1_x <- direction_x_indices[direc1]
          d2_x <- direction_x_indices[direc2]
          
          
          pval_df <- data.frame(Subfield=subfield_name,
                                N=nrow(df1),
                                Comp=sprintf("%s-%s",direc1,direc2),
                                statistic=wilk$statistic,
                                pval=wilk$p.value,
                                alternative=wilk$alternative,
                                method=wilk$method,
                                corrected_pval=wilk_correction,
                                signif=signif.num(wilk$p.value),
                                correcSignif=signif.num(wilk_correction),
                                conf_name=conf_name,
                                continuouity_correction=correct_continuouity)
          
          statistics_df <- rbind(statistics_df,
                                 pval_df)
          
          if (correct_continuouity) {
            gfinal <-
              gfinal + 
              geom_text(data=data.frame(x = (d1_x + d2_x) / 2,
                                        y = 1.05,
                                        label=signif.num(wilk_correction)),
                        aes(x=x,y=y,label=label))
          }
        }
      }
      
      for (size_name in names(sizes)) {
        size = sizes[[size_name]]
        pdf(file=sprintf("%s\\%s\\%s_%s_signif_codes.pdf",
                         write_path,
                         size_name,
                         subfield_name,
                         conf_name),
            height=size,
            width=size)
        
        plot(gfinal)
        dev.off()
        
        pdf(file=sprintf("%s\\%s\\%s_%s.pdf",
                         write_path,
                         size_name,
                         subfield_name,
                         conf_name),
            height=size,
            width=size)
        
        plot(g_est)
        dev.off()
      }
    }
  }
  
  
  write.csv(file=sprintf("%s//statistics_df.csv",write_path), statistics_df)
  write.csv(file=sprintf("%s//values_df.csv",write_path), values_df)
}

supp_figure_different_cost_functions_likelihood_plots <- function() {
  
  write_path <- sprintf("%s\\supp_figure_different_cost_functions\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\supp_figure_different_cost_functions\\likelihood_plots\\",figures_path)
  dir.create(write_path)
  dir.create(sprintf("%s\\supp_figure_different_cost_functions\\likelihood_plots\\big\\",figures_path))
  dir.create(sprintf("%s\\supp_figure_different_cost_functions\\likelihood_plots\\small\\",figures_path))
  dir.create(sprintf("%s\\supp_figure_different_cost_functions\\likelihood_plots\\medium\\",figures_path))
  dir.create(sprintf("%s\\supp_figure_different_cost_functions\\likelihood_plots\\medium_1\\",figures_path))
  dir.create(sprintf("%s\\supp_figure_different_cost_functions\\likelihood_plots\\medium_2\\",figures_path))
  
  
  sizes=c(big=2.5,
          medium=2,
          medium_1=1.75,
          medium_2=1.6,
          small=1.25)
  
  dir.create(write_path)
  
  true_data_paths <- sprintf("%s\\%s", all_data_paths[1:18], "equalized")
  
  sessions_to_use=c(5:8,13:16)
  
  all_df <- get_all_df(true_data_paths, sessions_to_use)
  
  estimation_configurations <- list(JSD_cyclic=list(path="JSD_simulations_likelihood_c", edge_free=F),
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
                  size=.3, linetype="dashed", color="#132046")  +
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
      }
    }
  }
}
