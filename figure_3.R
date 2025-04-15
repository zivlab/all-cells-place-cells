figure_3_histograms <- function() {
  write_path <- sprintf("%s\\figure_3\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_3\\histograms\\",figures_path)
  dir.create(write_path)
  
  nreps=15
  
  tests <- c(KS=ks.test)#,
             #wilcox=wilcox.test)
  
  visualization_splines <- c(raw=-1,
                             #`spline_1e_4`=1e-4,
                             `spline_1e_5`=1e-5,
                             `spline_2e_5`=2e-5,
                             `spline_5e_5`=5e-5)#,
                             #`spline_8e_5`=8e-5)
  
  
  cost_functions <- list(KS=list(likelihood="KS_simulations_likelihood_c",
                                 estimation="KS_simulations_estimate"))#,
                         # JSD=list(likelihood="JSD_simulations_likelihood_c",
                         #          estimation="JSD_simulations_estimate"))
  
  
  for (np in all_data_paths[2:8]) { 
    for (ses in c(5:8,13:16)) {
      
      session_path <- sprintf("%s\\equalized\\session_%d", np, ses)
      b <- get_spike_train_and_stim_trace_from_path(np, ses, equalize_prior=T)
      
      
      for (cost_func_name in names(cost_functions)) {
        
        params_ex <- get_fit_params(session_path,
                                    likelihood_path = cost_functions[[cost_func_name]]$likelihood,
                                    estimation_path = cost_functions[[cost_func_name]]$estimation)
        params <- params_ex$params
        
        print(params_ex)
        
        for (statistical_test in names(tests)) {
          for (vis_spline in names(visualization_splines)) {
            
            test_to_use <- tests[[statistical_test]]
            vis_spline_to_use <- visualization_splines[vis_spline]
            
            
            dir.create(sprintf("%s\\%s_%s\\",
                               write_path,
                               vis_spline,
                               statistical_test))
                               
                               
            dir.create(sprintf("%s\\%s_%s\\%s",
                               write_path,
                               vis_spline,
                               statistical_test,
                               cost_func_name))
            
            dir.create(sprintf("%s\\%s_%s\\%s\\big_histograms",
                               write_path,
                               vis_spline,
                               statistical_test,
                               cost_func_name))
            dir.create(sprintf("%s\\%s_%s\\%s\\small_histograms",
                               write_path, 
                               vis_spline,
                               statistical_test,
                               cost_func_name))
            
            for (i in 1:nreps) {
              
              dir_path_big_hists <- sprintf("%s\\%s_%s\\%s\\big_histograms\\%d",
                                            write_path,
                                            vis_spline,
                                            statistical_test,
                                            cost_func_name,
                                            i)
              
              dir_path_small_hists <- sprintf("%s\\%s_%s\\%s\\small_histograms\\%d",
                                              write_path, 
                                              vis_spline,
                                              statistical_test,
                                              cost_func_name,
                                              i)
              dir.create(dir_path_big_hists)
              dir.create(dir_path_small_hists)
              
              graphs <- plot_gen_place_cell_histograms_new(average_width = params["average_width"],
                                                           sd_width = params["sd_width"],
                                                           noise_level = params["noise"],
                                                           double_peak_percent = params["double_peak_pct"],
                                                           true_spike_train = b[[1]],
                                                           stim_trace = b[[2]],
                                                           place_cell_percentage = params_ex$pct,
                                                           smooth_lambda = -1,
                                                           breaks_peaks=9,
                                                           breaks_SI=21,
                                                           pval_func = test_to_use,
                                                           visualization_spline = vis_spline_to_use,
                                                           density_p = T,
                                                           ttle_size = 7,
                                                           title=T)
              
              thm1 <- theme(plot.margin=unit(c(0,0,0,0), "cm"))
              thm_bigger_text <- theme(plot.margin=unit(c(0,0,0,0), "cm"),
                                       text=element_text(size=14))
              graphs_t <- lapply(graphs, function(pl) {pl + thm1})
              small_graphs_t  <- lapply(graphs, function(pl) {pl + thm1 + ggtitle("")})
              graphs_t$nrow <- 2
              small_graphs_t$nrow  <- 1
              
              
              graphs_t_bigger_text <- lapply(graphs, function(pl) {pl + thm_bigger_text})
              small_graphs_t_bigger_text  <- lapply(graphs, function(pl) {pl + thm_bigger_text + ggtitle("")})
              graphs_t_bigger_text$nrow <- 2
              small_graphs_t_bigger_text$nrow  <- 1
              
              hists_p <- do.call(arrangeGrob, graphs_t)
              small_hists_p <- do.call(arrangeGrob, small_graphs_t)
              hists_p_bigger_text <- do.call(arrangeGrob, graphs_t_bigger_text)
              small_hists_p_bigger_text <- do.call(arrangeGrob, small_graphs_t_bigger_text)
              
              if (grepl("CA3", np)) {
                mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', np))
                mice_str <- substr(np, mice_str_index, mice_str_index+4) 
              } else {
                mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', np))
                mice_str <- substr(np, mice_str_index, mice_str_index+3) 
              }
              
              dir <- ifelse(grepl("Right", np), "R", "L")
              
              pdf(sprintf("%s\\simulation_hists_%s_%d_%s.pdf",
                          dir_path_big_hists,
                          mice_str, 
                          ses, 
                          dir), 
                  height=4, width=4); 
              plot(hists_p); 
              dev.off() 
              
              pdf(sprintf("%s\\simulation_hists_%s_%d_%s.pdf",
                          dir_path_small_hists,
                          mice_str, 
                          ses, 
                          dir), 
                  height=1.5, width=6); 
              plot(small_hists_p); 
              dev.off() 

              pdf(sprintf("%s\\simulation_hists_%s_%d_%s_bigger_text.pdf",
                          dir_path_big_hists,
                          mice_str, 
                          ses, 
                          dir), 
                  height=4, width=4); 
              plot(hists_p_bigger_text); 
              dev.off() 
              
              pdf(sprintf("%s\\simulation_hists_%s_%d_%s_bigger_text.pdf",
                          dir_path_small_hists,
                          mice_str, 
                          ses, 
                          dir), 
                  height=1.5, width=6); 
              plot(small_hists_p_bigger_text); 
              dev.off()               
            }
          }
        }
      }
    }
  }
}

figure_3_likelihood_plots <- function() {
  write_path <- sprintf("%s\\figure_3\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_3\\likelihood_plots\\",figures_path)
  dir.create(write_path)
  
  cost_functions <- list(KS=list(likelihood="KS_simulations_likelihood_c",
                                 estimation="KS_simulations_estimate"),
                         JSD=list(likelihood="JSD_simulations_likelihood_c",
                                  estimation="JSD_simulations_estimate"))
  
  for (cost_func_name in names(cost_functions)) {
    write_path_f <- sprintf("%s\\%s", write_path, cost_func_name)
    dir.create(write_path_f)
    
    for (np in all_data_paths) { 
      for (ses in 1:16) {
        
        session_path <- sprintf("%s\\equalized\\session_%d", np, ses) 
        likelihood_df <- get_fit_params(session_path,
                                        just_likelihood = T,
                                        likelihood_path = cost_functions[[cost_func_name]]$likelihood,
                                        estimation_path = cost_functions[[cost_func_name]]$estimation)
        
        measured_percent <- likelihood_df[1,"measured_pct"]
        
        likelihood_vec <- likelihood_df[,"likelihood"]
        likelihood_df <- likelihood_df[,-which(colnames(likelihood_df) == "likelihood")] # Remove last
        likelihood_df <- likelihood_df[,-which(colnames(likelihood_df) == "measured_pct")] # Remove last
        
        
        # This is just to re-extract measured place cell percentage,
        # In older versions of the analysis this was saved with the df as well
        rownames(likelihood_df) <- rev(rownames(likelihood_df))
        names(likelihood_vec) <- rownames(likelihood_df)
        mdf <- melt(likelihood_df)
        mdf$Var1 <- factor(mdf$Var1)
        mdf$SimP <- factor(sprintf("%.2f",as.numeric(as.vector(mdf$Var1))),
                           levels=sprintf("%.2f",seq(0.4,1,by=0.1)))
        
        mdf$MeasuredP <- mdf$value
        
        bar_df <- data.frame(Likelihood=likelihood_vec, 
                             Percentage=factor(sprintf("%.2f",as.numeric(names(likelihood_vec))),
                                               levels=sprintf("%.2f",seq(0.4,1,by=0.1))))
        
        bar_df$Likelihood <- bar_df$Likelihood / sum(bar_df$Likelihood)
        
        limits_of_plots <- c(min(mdf$value) - 0.05, max(mdf$value) + 0.1)
        gbox <- 
          ggplot(mdf) +
          geom_violin(aes(y=MeasuredP, group=SimP, fill=SimP, x=SimP), size=1, alpha=0.5) + 
          geom_jitter(aes(y=MeasuredP, group=SimP, x=SimP), size=0.05, position=position_jitter(0.15), alpha=0.25)  +
          stat_summary(fun=mean, geom="point", size=1.5, color="red", aes(x=SimP, y=MeasuredP)) + 
          scale_fill_brewer(palette="RdYlBu") +
          theme_light() +
          xlab("Ground Truth fraction") +
          ylab("Measured fraction") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.text.x = element_text(angle = 45, vjust=0.5),
                legend.position="NA",
                panel.border = element_blank(),
                text = element_text(size = 14),
                panel.background = element_blank()) +
          geom_hline(yintercept=(measured_percent), linetype="dashed", col="gray60", size=2) +
          geom_text(y=(measured_percent + 0.01), x=2.5, label="Real data",col="gray60",size=3) 
        
        gbar <- ggplot(bar_df, aes(x=Percentage, y=Likelihood)) +
          geom_bar(stat="identity", fill="gray50", color="black", width=0.5) + 
          theme_light() +
          xlab("Ground Truth fraction") +
          ylab("Normalized Likelihood") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle = 45, vjust=0.5),
                axis.line = element_line(colour = "black"),
                legend.position="NA",
                panel.border = element_blank(),
                text = element_text(size = 14),   
                panel.background = element_blank())
        
        margin <- theme(plot.margin=unit(c(0,0,0,0,0), "cm"))
        margin_big_text <- theme(plot.margin=unit(c(0,0,0,0,0), "cm"),
                                 text=element_text(size=15))
        
        rest <- grid.arrange(gbox + margin, gbar + margin, nrow=1)
        rest_15 <- grid.arrange(gbox + margin_big_text, 
                                gbar + margin_big_text, nrow=1)
        
        if (grepl("CA3", np)) {
          mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', np))
          mice_str <- substr(np, mice_str_index, mice_str_index+4) 
        } else {
          mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', np))
          mice_str <- substr(np, mice_str_index, mice_str_index+3) 
        }
        
        dir <- ifelse(grepl("Right", np), "R", "L")
        
        pdf(sprintf("%s\\box_and_bar_%s_%d_%s.pdf",write_path_f, mice_str, ses, dir), 
            height=3, width=5); 
        plot(rest); 
        dev.off() 
        
        pdf(sprintf("%s\\box_and_bar_%s_%d_%s_bigger_text.pdf",write_path_f, mice_str, ses, dir), 
            height=3, width=5); 
        plot(rest_15); 
        dev.off()         
      }
    }
  }
}

figure_3_simulation_vs_real_realisations <- function() {
  write_path <- sprintf("%s\\figure_3\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_3\\sim_vs_real_realisations\\",figures_path)
  dir.create(write_path)
  
  samples_path <- sprintf("%s\\samples", base_path)
  sim_vs_real_sampels <- list.files(samples_path, full.names = T)
  for (sample_path in sim_vs_real_sampels) {
    load(sample_path, verbose=T)
    
    simulation_vs_real_df <- final_result_2
    names(simulation_vs_real_df)[[7]] <- "real"
    simulation_vs_real_df$real$Time <- 
      simulation_vs_real_df$real$Time / max(simulation_vs_real_df$real$Time) * 20
    
    g_real_cyc <-  ggplot(simulation_vs_real_df$real, aes(x=Time, y=Cyclic)) + 
      geom_line(size=1) + 
      theme_light() +
      ylim(c(0, 0.95)) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position="NA",
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.margin=unit(c(0,0,0,0), "cm")) +
      labs(x="Time (m)", y="Fraction")
    #labs(x="Sample duration (minutes)", y = "Fraction of place cells (%)") +
    #xlim(0, max(final_result[[7]][,"Time"]))
    
    
    lev <- rev(as.numeric(names(simulation_vs_real_df)[-len(names(simulation_vs_real_df))]))
    
    sim_pct_colors = brewer.pal(n = len(lev),  name = "RdYlBu")
    names(sim_pct_colors) <- lev
    plots_list <- list()
    
    for (sim_pct in lev) {
      print(sim_pct)
      simulation_df <- simulation_vs_real_df[[as.character(sim_pct)]]
      g_sim_cyc <- g_real_cyc
      
      for (i in 1:(ncol(simulation_df) / 2)) {
        tmp_df_cyc <- data.frame(Fraction=simulation_df[, (i * 2)],
                                 Time=simulation_vs_real_df$real$Time)
        
        g_sim_cyc <- 
          g_sim_cyc + geom_line(data=tmp_df_cyc,
                                aes(x=Time, y=Fraction),
                                color=sim_pct_colors[as.character(sim_pct)],
                                size=1,
                                alpha=0.5) 
        
      
      }
      
      g_f <- g_sim_cyc + theme(text=element_text(size=15))
      
      plots_list <- append(plots_list, list(g_f))
    }
    
    plots_list <- lapply(plots_list, function(pl) {pl + theme(plot.margin=unit(c(0,0,0,0), "cm"))})
    plots_list$nrow <- 2
    #plots_list$left <-"Fraction of place cells (%)"
    #plots_list$bottom <- "Sample duration (minutes)"
    sim_real_p <- do.call(grid.arrange, plots_list)
    
    if (grepl("CA3", sample_path)) {
      mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', sample_path))
      mice_str <- substr(sample_path, mice_str_index, mice_str_index+4) 
    } else {
      mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', sample_path))
      mice_str <- substr(sample_path, mice_str_index, mice_str_index+3) 
    }
    
    session_str_index <- unlist(gregexpr("session_", sample_path))
    
    session <- substr(sample_path, 
                      session_str_index + str_length("session_"), 
                      session_str_index + str_length("session_") + 1)
    
    session <- str_replace(session, "_", "")
    
    dir <- ifelse(grepl("Right", sample_path), "R", "L")
    cost_func <- ifelse(grepl("KS", sample_path), "KS", "JSD")
    
    dir.create(sprintf("%s//%s", write_path, cost_func))
    dir.create(sprintf("%s//%s//big//", write_path, cost_func))
    dir.create(sprintf("%s//%s//small//", write_path, cost_func))
    pdf(file=sprintf("%s/%s/big/sim_vs_real_%s_%s_%s.pdf",
                     write_path,
                     cost_func,
                     mice_str,
                     session,
                     dir),
        height=3,
        width=1.5 * 3)
    
    plot(sim_real_p)
    dev.off()
    
    pdf(file=sprintf("%s/%s/small/sim_vs_real_%s_%s_%s.pdf",
                     write_path,
                     cost_func,
                     mice_str,
                     session,
                     dir),
        height=2.5,
        width=1.35 * 3)
    
    plot(sim_real_p)
    dev.off()    
  }
}

figure_3_simulation_vs_all_sessions <- function() {
  write_path <- sprintf("%s\\figure_3\\",figures_path)
  dir.create(write_path)
  write_path <- sprintf("%s\\figure_3\\sim_vs_real_all_sessions\\",figures_path)
  dir.create(write_path)
  dir.create(sprintf("%s//pvalue", write_path))
  dir.create(sprintf("%s//duration", write_path))
  
  subfield_path_list <- list(CA1=1:8,
                             CA3=9:18)
  for (subfield_name in names(subfield_path_list)) { 
    
  true_data_paths <- sprintf("%s\\%s", all_data_paths[subfield_path_list[[subfield_name]]], "equalized")
  ses_ind <- c(1:16)
  #by_activity <- plot_place_cell_pct_by_activity_plot(true_data_paths, cyclic_only=T, fit_sim = "fit_simulation_new_v_fit\\")
  
  
  if (subfield_name == "CA1") {
    configurations <- list(list(ext="", prefix=""),
                           list(ext="_new", prefix=""),
                           list(ext="", prefix="JSD_"),
                           list(ext="", prefix="JSD_2_"),
                           list(ext="", prefix="JSD_3_"),
                           list(ext="", prefix="JSD_4_")) 
  } else {
    configurations <- list(list(ext="", prefix=""),
                           list(ext="", prefix="JSD_"),
                           list(ext="", prefix="JSD_2_"),
                           list(ext="", prefix="JSD_3_"),
                           list(ext="", prefix="JSD_4_"))    
  }
  
  

  for (conf in configurations){
  by_duration_simulation <- 
        plot_place_cell_pct_by_sample_duration_figure(ses_ind = ses_ind, 
                                                      true_data_paths, 
                                                      bin_size = 2.5, 
                                                      ext = conf$ext,
                                                      prefix = conf$prefix,
                                                      cyclic_only=T, 
                                                      fit_sim = "fit_simulation\\", verbose=F)
  
  by_duration_real <- plot_place_cell_pct_by_sample_duration_figure(true_data_paths, bin_size = 2.5, 
                                                                    ext = "", 
                                                                    cyclic_only=T, 
                                                                    fit_sim = "", 
                                                                    verbose=F, 
                                                                    ses_ind = ses_ind)
  
  overlaid_duration <- by_duration_simulation[[1]] + 
                         geom_line(data=by_duration_real[[4]], 
                                   linetype="dashed", 
                                   aes(x=time, y=fraction), 
                                   size=1, 
                                   color="black") + 
                         geom_ribbon(data=by_duration_real[[4]], 
                                     aes(ymin=fraction-sd, ymax=fraction+sd), 
                                     fill=adjustcolor("black", alpha=0.3), 
                                     color=NA) +
                         theme(text=element_text(size=15))
  
  pdf(file=sprintf("%s//duration//%s_%s%sduration_real_vs_sim.pdf",
                   write_path,
                   subfield_name,
                   conf$ext,
                   conf$prefix),
      height=2.5,
      width=2.5)
  
  plot(overlaid_duration)
  dev.off()
  }
  
  by_pval_simulation <- plot_pval_by_sample_duration_figure(true_data_paths, 
                                                            ext = "", 
                                                            fit_sim = "fit_simulation\\", 
                                                            verbose=F, 
                                                            prefix="JSD_1_",
                                                            type1=F, 
                                                            ses_ind = ses_ind)

  by_pval_real <- plot_pval_by_sample_duration_figure(true_data_paths,
                                                      ext = "", 
                                                      fit_sim = "", 
                                                      prefix="JSD_1_",
                                                      verbose=F, type1 = F, 
                                                      ses_ind = ses_ind)
  
  
  #overlaid_pval <- 
          by_pval_simulation[[1]] + 
                    geom_line(data=by_pval_real[[2]], 
                              linetype="dashed", 
                              aes(x=time, y=pvd), 
                              size=1, color="black") + 
                    geom_ribbon(data=by_pval_real[[2]], 
                                aes(ymin=pvd-pvdsd, ymax=pvd+pvdsd), 
                                fill=adjustcolor("black", alpha=0.3), 
                                color=NA) +
                    theme(text=element_text(size=15))
  
  plot(overlaid_pval)
  pdf(file=sprintf("%s//pvalue//%s_pval_real_vs_sim.pdf",
                   write_path,
                   subfield_name),
      height=2.5,
      width=2.5)
  
  plot(overlaid_pval)
  dev.off()
  }
}




validation_figures <- function() {
  
  
  
  prfx = "09"
  path = sprintf("/Users/itayta/Downloads/yaniv_place_cells_project//simulated_dataset%s", prfx)
  
  load(sprintf("%s/%s", path, '\\KS_simulations_likelihood_c\\spike_train'), verbose=T)
  load(sprintf("%s/%s", path, '\\KS_simulations_likelihood_c\\stim_trace'), verbose=T)
  load(sprintf("%s/%s", path, '\\KS_simulations_likelihood_c\\pred_df.R'), verbose=T)
  #load(sprintf("%s/%s", path, 'JSD_simulations_estimate//param_fit_0.700\\\\pct_estimate_df.R'), verbose=T)
  
  fit_params_df <- extract_simulation_params(path, 
                                             plot=F,  
                                             absolute_min = T, 
                                             estimation_path = estimation_folder, 
                                             pct_range =  seq(0.5, 1, by=0.1),
                                             file_name="pct_estimate_df.R") 
  
  test_to_use <- ks.test
  vis_spline_to_use <- 2e-5
  
  params <- fit_params_df["1.0",]
  params_ex <- list()
  params_ex$pct <- 1
  
  b <- list(spike_train, stim_trace)
  
  graphs <- plot_gen_place_cell_histograms_new(average_width = params["average_width"],
                                               sd_width = params["sd_width"],
                                               noise_level = params["noise"],
                                               double_peak_percent = params["double_peak_pct"],
                                               true_spike_train = b[[1]],
                                               stim_trace = b[[2]],
                                               place_cell_percentage = params_ex$pct,
                                               smooth_lambda = -1,
                                               breaks_peaks=9,
                                               breaks_SI=21,
                                               pval_func = test_to_use,
                                               visualization_spline = vis_spline_to_use,
                                               density_p = T,
                                               ttle_size = 7,
                                               title=T)
  
  thm1 <- theme(plot.margin=unit(c(0,0,0,0), "cm"))
  thm_bigger_text <- theme(plot.margin=unit(c(0,0,0,0), "cm"),
                           text=element_text(size=14))
  graphs_t <- lapply(graphs, function(pl) {pl + thm1})
  small_graphs_t  <- lapply(graphs, function(pl) {pl + thm1 + ggtitle("")})
  graphs_t$nrow <- 2
  small_graphs_t$nrow  <- 1
  
  
  graphs_t_bigger_text <- lapply(graphs, function(pl) {pl + thm_bigger_text})
  small_graphs_t_bigger_text  <- lapply(graphs, function(pl) {pl + thm_bigger_text + ggtitle("")})
  graphs_t_bigger_text$nrow <- 2
  small_graphs_t_bigger_text$nrow  <- 1
  
  hists_p <- do.call(arrangeGrob, graphs_t)
  small_hists_p <- do.call(arrangeGrob, small_graphs_t)
  hists_p_bigger_text <- do.call(arrangeGrob, graphs_t_bigger_text)
  small_hists_p_bigger_text <- do.call(arrangeGrob, small_graphs_t_bigger_text)
  
  plot(hists_p)
  
  
  final_df <- final_df[6:1,]
  likelihood_df <- final_df
  
  measured_percent <- likelihood_df[1,"measured_pct"]
  
  likelihood_vec <- likelihood_df[,"likelihood"]
  likelihood_df <- likelihood_df[,-which(colnames(likelihood_df) == "likelihood")] # Remove last
  likelihood_df <- likelihood_df[,-which(colnames(likelihood_df) == "measured_pct")] # Remove last
  
  
  # This is just to re-extract measured place cell percentage,
  # In older versions of the analysis this was saved with the df as well
  #rownames(likelihood_df) <- rev(rownames(likelihood_df))
  names(likelihood_vec) <- rownames(likelihood_df)
  mdf <- melt(likelihood_df)
  mdf$Var1 <- factor(mdf$Var1)
  mdf$SimP <- factor(sprintf("%.2f",as.numeric(as.vector(mdf$Var1))),
                     levels=sprintf("%.2f",seq(0.4,1,by=0.1)))
  
  mdf$MeasuredP <- mdf$value
  
  bar_df <- data.frame(Likelihood=likelihood_vec, 
                       Percentage=factor(sprintf("%.2f",as.numeric(names(likelihood_vec))),
                                         levels=sprintf("%.2f",seq(0.4,1,by=0.1))))
  
  bar_df$Likelihood <- bar_df$Likelihood / sum(bar_df$Likelihood)
  
  limits_of_plots <- c(min(mdf$value) - 0.05, max(mdf$value) + 0.1)
  gbox <- 
    ggplot(mdf) +
    geom_violin(aes(y=MeasuredP, group=SimP, fill=SimP, x=SimP), size=1, alpha=0.5) + 
    geom_jitter(aes(y=MeasuredP, group=SimP, x=SimP), size=0.05, position=position_jitter(0.15), alpha=0.25)  +
    stat_summary(fun=mean, geom="point", size=1.5, color="red", aes(x=SimP, y=MeasuredP)) + 
    scale_fill_brewer(palette="RdYlBu") +
    theme_light() +
    xlab("Ground Truth fraction") +
    ylab("Measured fraction") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, vjust=0.5),
          legend.position="NA",
          panel.border = element_blank(),
          text = element_text(size = 14),
          panel.background = element_blank()) +
    geom_hline(yintercept=(measured_percent), linetype="dashed", col="gray60", size=2) +
    geom_text(y=(measured_percent + 0.01), x=2.5, label="Real data",col="gray60",size=3) 
  
  gbar <- ggplot(bar_df, aes(x=Percentage, y=Likelihood)) +
    geom_bar(stat="identity", fill="gray50", color="black", width=0.5) + 
    theme_light() +
    xlab("Ground Truth fraction") +
    ylab("Normalized Likelihood") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, vjust=0.5),
          axis.line = element_line(colour = "black"),
          legend.position="NA",
          panel.border = element_blank(),
          text = element_text(size = 14),   
          panel.background = element_blank())
  
  margin <- theme(plot.margin=unit(c(0,0,0,0,0), "cm"))
  margin_big_text <- theme(plot.margin=unit(c(0,0,0,0,0), "cm"),
                           text=element_text(size=15))
  
  rest <- grid.arrange(gbox + margin, gbar + margin, nrow=1)
  rest_15 <- grid.arrange(gbox + margin_big_text, 
                          gbar + margin_big_text, nrow=1)
  
  fpath = "/Users/itayta/Downloads/yaniv_place_cells_project/allcells_figures/new_figures_revisiting/final_figures_for_submission/supplemental_figures/validation_figure"
  
  
  pdf(sprintf("%s/%s_hists.pdf", fpath, prfx), 
      height=4, width=4); 
    plot(small_hists_p_bigger_text); 
  dev.off() 
  
  pdf(sprintf("%s/%s_lklhood_15.pdf", fpath, prfx), 
      height=4, width=8); 
  plot(rest_15); 
  dev.off() 
  
  pdf(sprintf("%s/%s_lklhood.pdf", fpath, prfx), 
      height=4, width=8); 
  plot(rest); 
  dev.off() 
}

  
  
