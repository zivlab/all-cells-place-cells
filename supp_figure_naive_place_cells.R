
supp_figure_naive_place_cells_boxplots  <- function () {

  write_path <- sprintf("%s\\supp_figure_naive_pcs\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\supp_figure_naive_pcs\\boxplots\\",figures_path)
  
  directions_to_use = c("Both", "Left", "Right")
  sizes=c(big=3.5,
          medium=3,
          medium_1=2.7,
          small=2.5)
  
  dir.create(sprintf("%s\\", write_path))
  dir.create(sprintf("%s\\big\\", write_path))
  dir.create(sprintf("%s\\medium\\", write_path))
  dir.create(sprintf("%s\\medium_1\\", write_path))
  dir.create(sprintf("%s\\small\\", write_path))
  
  
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
      theme(text=element_text(size=15, color="black"),
            axis.text = element_text(color="black")) +
      xlab(x_axes) +
      ylab("Fraction") + 
      ylim(0,1.05)
    
    return(gf)
  }
  
  
  
  naive_configurations <- list(cyclic_MI_eq_prior_either_dir=list(path="percent_df_cyclic_MI_eq_prior.R", 
                                                                  weighted=F, 
                                                                  multiple_pvals=F,
                                                                  title="MI"),
                               cyclic_non_eq_prior_either_dir=list(path="percent_df_cyclic_SI_reg.R", 
                                                                   weighted=F, 
                                                                   multiple_pvals=F,
                                                                   title="Non-unifrom prior"))
                            
  values_df <- c()  
  statistics_df <- c()
  
  for (naive_conf_name in names(naive_configurations)) {
    naive_conf = naive_configurations[[naive_conf_name]]
    
    all_df <- get_all_df(true_data_paths, 
                         sessions_to_use,
                         file_name = naive_conf$path,
                         multiple_pvals=naive_conf$multiple_pvals,
                         verbose=T)
    
    
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
    
    
    for (subfield_name in c("CA1", "CA3")){ 
      for (df_idx in 1:3) {
        
        df_to_use <- dfs[[df_idx]]
        colnames(df_to_use) <- c("Session", "Mice", "Subfield", "Direction", "Percent")
        
        df_to_use <- df_to_use[df_to_use$Subfield %in% subfield_name,]
        
        # Set second env session count refrence to 1...8
        df_to_use$Session <- df_to_use$Session %% 8
        df_to_use$Session[df_to_use$Session == 0] <- 8
        
        
        stats_df <- c()
        for (direc in unique(df_to_use$Direction)) {
          direc_df <- df_to_use[df_to_use$Direction == direc,]
          stats_df <- cbind(stats_df,
                            direc_df[order(paste(direc_df$Session, direc_df$Mice)),"Percent"])
        }
        
        colnames(stats_df) <- unique(df_to_use$Direction)
        combination_matrix <- combn(len(colnames(stats_df)), 2)
        
        for (comparision in 1:ncol(combination_matrix)) {
          idx1 <- combination_matrix[1, comparision]
          idx2 <- combination_matrix[2, comparision]
          wilc <- wilcox.test(stats_df[,idx1], stats_df[,idx2], paired=T, correct=F)
          direc_1 <- colnames(stats_df)[idx1]
          direc_2 <- colnames(stats_df)[idx2]
          pvals_df <- data.frame(Comp=sprintf("%s-%s",
                                              direc_1,
                                              direc_2),
                                 method=wilc$method,
                                 W=wilc$statistic,
                                 Pval=wilc$p.value,
                                 Alt=wilc$alternative,
                                 PvalCorr=wilc$p.value * ncol(combination_matrix),
                                 Signif=signif.num(wilc$p.value),
                                 CorrSign=signif.num(ifelse(wilc$p.value * ncol(combination_matrix) > 1, 1,wilc$p.value * ncol(combination_matrix))),
                                 Subfield=subfield_name,
                                 Df_name=plot_titles[[df_idx]],
                                 ConfName=naive_conf_name)
          statistics_df <- rbind(statistics_df,pvals_df)
        }
        
        
        #df_to_use <- df_to_use[df_to_use$Direction %in% directions_to_use, ]
        
        sd_df <- ddply(df_to_use, .(Subfield),
                       function(subfield_df) {
                         return(ddply(subfield_df, .(Direction), 
                                      function(sub_df) {
                                        return(c(mean(sub_df[,"Percent"]), sem(sub_df[,"Percent"])))}))
                       })
        colnames(sd_df) <- c("Subfield", "Direction", "Percent", "Sd")
        
        sd_df$df_name <- plot_titles[df_idx]
        sd_df$conf <- naive_conf_name
        values_df <- rbind(values_df, sd_df)
        
        
        g <- ggplot(sd_df, aes(x=Direction, y=Percent)) 
        x_axes = "Running direction"
        
        
        g <- 
          g +
          geom_hline(yintercept=1, linetype="dashed", size=0.8) +
          geom_boxplot(data=df_to_use, 
                       aes(group=interaction(Direction,Subfield),
                           color=interaction(Direction,Subfield),
                           fill=interaction(Direction,Subfield)), 
                       width=0.5,
                       position = position_dodge(0.75),
                       size=1)
        # geom_jitter(data=df_to_use,
        #             aes(y=Percent, 
        #                 group=interaction(Direction,Subfield), 
        #                 fill=interaction(Direction,Subfield)),
        #             position=position_jitterdodge(1), 
        #             color="gray20", 
        #             size=0.75, 
        #             alpha=0.1) 
        
        g <- color_plot(g, x_axes) + 
          ggtitle(sprintf("%s - %s",
                          plot_titles[df_idx],
                          naive_conf$title)) +
          theme(plot.title=element_text(size=11))
        
        
        
        for (size_name in names(sizes)) {
          size = sizes[[size_name]]
          
          dir.create(sprintf("%s\\%s\\%s\\",
                             write_path,
                             size_name,
                             plot_names[[df_idx]]))
          
          pdf(file=sprintf("%s\\%s\\%s\\%s_%s.pdf",
                           write_path,
                           size_name,
                           plot_names[[df_idx]],
                           subfield_name,
                           naive_conf_name),
              height=size,
              width=size * 1.5)
          
          plot(g)
          dev.off()
        }
      }
    }

  }
  
  write.csv(file=sprintf("%s\\Values.csv",write_path), values_df)
  write.csv(file=sprintf("%s\\Statistics.csv",write_path), statistics_df)
}

supp_figure_naive_place_cells_pvalue_range  <- function () {
  
  write_path <- sprintf("%s\\supp_figure_naive_pcs\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\supp_figure_naive_pcs\\pvalue_ranges\\",figures_path)
  
  directions_to_use = c("Both", "Left", "Right")
  sizes=c(big=2.5,
          big_1=2.25,
          medium=2,
          medium_1=1.75,
          medium_2=1.5,
          small=1.25)
  
  dir.create(sprintf("%s\\", write_path))
  dir.create(sprintf("%s\\big\\", write_path))
  dir.create(sprintf("%s\\big_1\\", write_path))
  dir.create(sprintf("%s\\medium\\", write_path))
  dir.create(sprintf("%s\\medium_1\\", write_path))
  dir.create(sprintf("%s\\medium_2\\", write_path))
  dir.create(sprintf("%s\\small\\", write_path))
  
  
  true_data_paths <- sprintf("%s\\%s", all_data_paths[1:18], "equalized")
  
  sessions_to_use=c(5:8,13:16)
  
  naive_configurations <- list(eq_prior=list(path="multiple_pvalS_percent_df_cyclic_SI_eq_prior.R", 
                                             weighted=F, 
                                             multiple_pvals=T),
                               regular_prior=list(path="multiple_pvals_percent_df_cyclic_SI_reg.R", 
                                                  weighted=F, 
                                                  multiple_pvals=T))
  
  
  for (naive_conf_name in names(naive_configurations)) {
    naive_conf = naive_configurations[[naive_conf_name]]
    
    all_df <- get_all_df(true_data_paths, 
                         sessions_to_use,
                         file_name = naive_conf$path,
                         multiple_pvals=naive_conf$multiple_pvals,
                         verbose=T)
    
    all_df_new <- c()
    for (pval in names(all_df)) {
      tmp <- all_df[[pval]]
      all_df_new <- rbind(all_df_new,
                          cbind(tmp, rep(as.numeric(pval), times=nrow(tmp))))
    }
    
    all_df <- all_df_new 
    colnames(all_df)[13] <- "Pvalue"
    
    percent_active_df <- all_df[,c(7:9)]
    percent_place_cells_df <- all_df[,c(1:3)]
    percent_all_df <- all_df[,c(4:6)]
    
    metadata_df <- all_df[,c("Session", "Mice", "Subfield", "Pvalue")]
    percent_active_df <- cbind(percent_active_df, metadata_df)
    percent_all_df <- cbind(percent_all_df, metadata_df)
    
    percent_place_cells_df <- cbind(percent_place_cells_df, metadata_df)  
    colnames(percent_active_df) <- c("Both", "Right", "Left", "Session", "Mice", "Subfield", "Pvalue")
    colnames(percent_place_cells_df) <- c("Both", "Right", "Left", "Session", "Mice", "Subfield", "Pvalue")
    colnames(percent_all_df) <- c("Both", "Right", "Left", "Session", "Mice", "Subfield", "Pvalue")
    
    melted_percent_active <- melt(percent_active_df, id.vars = c("Session", "Mice", "Subfield", "Pvalue"))
    melted_percent_place_cells <- melt(percent_place_cells_df, id.vars = c("Session", "Mice", "Subfield", "Pvalue"))
    melted_percent_all <- melt(percent_all_df, id.vars = c("Session", "Mice", "Subfield", "Pvalue"))
    
    ylabs = c("Percentage of active cells (%)",
              "Percentage of place cells (% of active cells)",
              "Percentage of place cells (% of all cells)")
    
    dfs <- list(melted_percent_active,
                melted_percent_place_cells,
                melted_percent_all)
    
    
    plot_names <- c("ActiveCells",
                    "PlaceCellsOfActive",
                    "PlaceCellsOfAll")
    
    
    plot_titles <- c("Active cells",
                     "Place cells (of active)",
                     "Place cells (of all)")
    
    
    
    for (sbf in unique(melted_percent_active$Subfield)) {
      for (df_idx in 2:3) {
        
        df_to_use <- dfs[[df_idx]]
        df_to_use <- df_to_use[df_to_use$Subfield == sbf,]
        colnames(df_to_use) <- c("Session", "Mice", "Subfield", "Pvalue", "Direction", "Percent")
        
        # Set second env session count refrence to 1...8
        df_to_use$Session <- df_to_use$Session %% 8
        df_to_use$Session[df_to_use$Session == 0] <- 8
        
        
        
        df_to_use <- df_to_use[df_to_use$Direction %in% directions_to_use, ]
        
        sd_df <- ddply(df_to_use, .(Pvalue),
                       function(pval_df) {
                         return(ddply(pval_df, .(Direction), 
                                      function(direction_df) {
                                        return(ddply(direction_df, .(Subfield),
                                                     function(sub_df) {
                                                        return(c(mean(sub_df[,"Percent"]), sd(sub_df[,"Percent"])))
                                                     }))
                                        }))
                       })
        colnames(sd_df) <- c("Pvalue", "Direction", "Subfield", "Percent", "Sd")
        
        
        g <- 

        ggplot(sd_df, aes(x=Pvalue, 
                          y=Percent, 
                          group=interaction(Subfield, Direction), 
                          color=interaction(Subfield, Direction))) + 
          geom_line(position=position_dodge(.005)) +
          geom_point(position=position_dodge(.005)) +
          geom_errorbar(aes(ymin=Percent-Sd, ymax=Percent+Sd),
                        position=position_dodge(.005),
                        width=0.02) +
          ylim(0,1.1) +
          geom_hline(yintercept=1, linetype="dashed", size=0.8) +
          base_plot_theme + 
          theme(text=element_text(size=15, color="black"),
                axis.text = element_text(color="black"),
                plot.title=element_text(size=11)) +
          ylab("Fraction") +
          xlab(TeX("$P_{value} Threshold$")) +
          ggtitle(sprintf("%s - %s", 
                          plot_titles[[df_idx]],
                          sbf)) +
          scale_color_manual(values=c("#F05A28", "#652D90", "gray20"),
                             breaks=sprintf("%s.%s", sbf,
                                            c("Right", "Left", "Both")))
          
        
        for (size_name in names(sizes)) {
          size = sizes[[size_name]]
          
          dir.create(sprintf("%s\\%s\\%s\\",
                             write_path,
                             size_name,
                             plot_names[[df_idx]]))
          
          pdf(file=sprintf("%s\\%s\\%s\\%s_%s_wide.pdf",
                           write_path,
                           size_name,
                           plot_names[[df_idx]],
                           sbf,
                           naive_conf_name),
              height=size,
              width=size * 1.5)
          
          plot(g)
          dev.off()
          
          pdf(file=sprintf("%s\\%s\\%s\\%s_%s.pdf",
                           write_path,
                           size_name,
                           plot_names[[df_idx]],
                           sbf,
                           naive_conf_name),
              height=size,
              width=size)
          
          plot(g)
          dev.off()
        }
      } 
    }
  }
}


