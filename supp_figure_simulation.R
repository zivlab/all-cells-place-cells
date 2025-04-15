


sample_posterior_new <- function(x1, x2, x1x1, sigma_n=.015) {
  
  
  x2x2 <- rbf_kernel(x2, sigma=sigma_n)
  x2x1 <- rbf_kernel(x2, Xtag=x1, sigma=sigma_n)
  x1x2 <-  rbf_kernel(x1, Xtag=x2, sigma=sigma_n)
  #x1x1 <- rbf_kernel(x1, sigma=sigma_n)
  
  y = rep(0, nrow(x1x1))
  gt <- poly_func(x2)
  invx2x2 <- inv(x2x2)
  
  mu_conditioned = x1x2 %*% invx2x2 %*% gt
  cov_mat_conditioned = x1x1 - x1x2 %*% invx2x2 %*% x2x1
    
  var <- unlist(lapply(1:ncol(cov_mat_conditioned), function(i) {cov_mat_conditioned[i,i]}))
  
  return(cbind(ypred=mu_conditioned, sdpred=var))
}
supp_figure_bayesian_optimizer_minimum <- function() {
  
  write_path <- sprintf("%s\\supp_figure_simulations\\",figures_path)
  dir.create(write_path)
  
  write_path <- sprintf("%s\\example_polynom", write_path)
  dir.create(write_path)
  dir.create(sprintf("%s\\big\\",write_path))
  dir.create(sprintf("%s\\medium\\",write_path))
  dir.create(sprintf("%s\\small\\", write_path))
  
  sizes=c(big=3,
          medium=2.5,
          small=2)
  
  ## approximate this function
  poly_func <- function(x)   {x ** 5 -(8 * x ** 3) + 10*x + 6}
  X_star <- as.matrix(seq(-2.5,2.5, length.out=2500))
  sgman=.5
  cov_mat_pred <- rbf_kernel(X_star, X_star, sigma=sgman)
  
  sampled <- sort(sample(500:1500, 5))
  X_train <- as.matrix(X_star[sampled,])
  
  x1x1 <- cov_mat_pred[-sampled,-sampled]
  
  print(X_train)
  print ("###################")
  
  y <- poly_func(X_train)
  
  
  plist <- list()
  
  mdf <- as.data.frame(cbind(as.vector(X_star), poly_func(X_star), rep(0, 2500)))
  colnames(mdf) <- c("x", "gt", "y")
  
  tmpdf <- as.data.frame(data.frame(x=X_train, y=y))
  
  
  g <- 
  ggplot(mdf, aes(x=x, y=gt)) + 
    geom_line() + geom_line(aes(x=x,y=y)) + 
    geom_ribbon(aes(x=x, ymin=y-2.5, ymax=y+2.5), color=NA, alpha=.5) +
    geom_point(data=tmpdf, aes(x=x,y=y), col="cyan", size=2) + 
    base_plot_theme +
    ylab("Y") + xlab("X")
  
  plot(g)
  
  plist <- append(plist, list(g))
  
  
  iter <- 1
  for (iter in 1:22){
    #iter <- iter + 1
    
    x1 = as.matrix(X_star[-sampled,])
    x2 = as.matrix(X_star[sampled,])
    x1x1 <- cov_mat_pred[-sampled,-sampled]
    
    res <- sample_posterior_new(x1=x1,
                                x2=x2,
                                x1x1 = x1x1,
                                sigma_n = sgman)

    
    min_per_x <- max(res[,1])
    mean_per_pred <- res[,1]
    sd_per_pred <- sqrt(res[,2])
    
    
    all_y_mean = c(mean_per_pred, y)[order(c(x1,x2))]
    all_y_sd = c(sd_per_pred, rep(0, times=len(y)))[order(c(x1,x2))]
    

    non_zero_sd_ind <- all_y_sd != 0
    gamma <- (min(all_y_mean) - all_y_mean + 0.01)
    
    Z <- rep(0, len(all_y_sd))
    Z[non_zero_sd_ind] <- gamma[non_zero_sd_ind] / all_y_sd[non_zero_sd_ind]
    
    EI <- rep(0, len(all_y_sd))
    EI[non_zero_sd_ind] <- (gamma[non_zero_sd_ind]) * pnorm(Z[non_zero_sd_ind]) + dnorm(Z[non_zero_sd_ind]) * all_y_sd[non_zero_sd_ind]
    
    plot(EI)
    
    #expected_improve <-  gamma * pnorm(gamma / (all_y_sd + 10 ** -40)) + all_y_sd * dnorm(gamma / (all_y_sd + 10 ** -40))
    
    #expected_improve <- all_y_sd * (gamma * pnorm(gamma) + dnorm(gamma))
    top_to_improve_exp <- order(EI, decreasing=T)
    
    
    addition_idx = 1
    left_to_add = 15
    new_X <- c()
    
    while (left_to_add > 0) {


      if (top_to_improve_exp[addition_idx] %in% sampled) {
        addition_idx <- addition_idx + 1
        next
      }

      new_X <- c(new_X, top_to_improve_exp[addition_idx])
      addition_idx <- addition_idx + 1
      left_to_add <- left_to_add - 1
    }
    
    #new_X <- top_to_improve_exp[1:left_to_add]
    
    
    
    mdf <- as.data.frame(cbind(as.vector(X_star), 
                               poly_func(X_star), 
                               all_y_mean,
                               6*all_y_sd,
                               EI))
    colnames(mdf) <- c("x", "gt", "y", "sd", "ei")
    
    tmpdf <- as.data.frame(data.frame(x=c(X_star)[sampled], y=y))
    newx <- as.data.frame(data.frame(x=c(X_star)[new_X], 
                                     y=poly_func(c(X_star)[new_X]),
                                     ei=EI[new_X]))
    
    g <- 
      ggplot(mdf, aes(x=x, y=gt)) + 
      geom_line() + geom_line(aes(x=x,y=y)) + 
      geom_ribbon(aes(x=x, ymin=y-sd, ymax=y+sd), color=NA, alpha=.5) +
      geom_point(data=tmpdf, aes(x=x,y=y), col="cyan", size=2) + 
      base_plot_theme +
      ylab("Y") + xlab("X") +
      geom_point(data=newx, aes(x=x,y=y), col="red") +
      #geom_vline(xintercept=X_star[new_X,][1], col="red", linetype="dashed", linewidth=.25) +
      #geom_vline(xintercept=X_star[new_X,][2], col="red", linetype="dashed", linewidth=.25) +
      ggtitle(sprintf("Iter: %d, Sampled:%d, Added:%d", iter, len(sampled), 2))
      
      
    g2 <- 
      ggplot(mdf, aes(x=x, y=ei)) + 
        geom_line() + 
        base_plot_theme +
        geom_point(data=newx, aes(x=x,y=ei), col="red") +
        #geom_vline(xintercept=X_star[new_X,][1], col="red", linetype="dashed", linewidth=.25) +
        #geom_vline(xintercept=X_star[new_X,][2], col="red", linetype="dashed", linewidth=.25) +
        ylab("Expected to improve") + xlab("X")
      
    gf <- plot_grid(g,g2,nrow=2, align = "y")
    plot(gf)
    plist <- append(plist, list(gf))
    
    
    # Add new points
    sampled <- sort(c(sampled,new_X))
    #sampled <- unique(sampled)
    #y <- c(y, poly_func(X_star[new_X,]))
    y <- poly_func(as.matrix(X_star[sampled,]))
    
    #print(sprintf("------- Added %d new points (%d Total)", nrow(new_X), nrow(X_train)))
  }
  
  res <- sample_posterior(as.matrix(X_train),
                          y,
                          as.matrix(X_star),
                          mean_func,
                          sigma_n=0.35,
                          n_sample=20,
                          K_xsxs = cov_mat_pred,
                          plot_kernel = F)
  
  
  # approx_df <- data.frame(cbind(apply(res,2, sd), colMeans(res), as.vector(X_star)))
  # colnames(approx_df) <- c("SD", "Predicted", "X")
  # measured_df <- data.frame(X=as.vector(X_train), Measured=y)
  # 
  # true_df <- data.frame(X=as.vector(X_star), Y=poly_func(as.vector(X_star)))
  # 
  # polynom_plot <- 
  #   ggplot(approx_df) + 
  #   geom_ribbon(aes(ymin=Predicted-SD, ymax=Predicted+SD, x=X), col=NA, fill="royalblue4", alpha=0.3) + 
  #   geom_line(aes(x=X, y=Predicted), col="royalblue4", size=.75)  + 
  #   geom_line(data=true_df, aes(x=X, y=Y), alpha=0.5, size=.75, linetype="dashed") +
  #   geom_point(data=measured_df, aes(x=X, y=Measured), col="red", stroke=0) +
  #   theme_light() +     
  #   theme(panel.grid.major = element_blank(), 
  #         panel.grid.minor = element_blank(),
  #         axis.line = element_line(colour = "black"),
  #         panel.border = element_blank(),
  #         panel.background = element_blank(),
  #         text=element_text(size=14)) + 
  #   geom_vline(xintercept=approx_df[which.min(approx_df$Predicted),]$X,
  #              size=.5,
  #              col="black",
  #              linetype="dashed") + 
  #     geom_hline(yintercept=approx_df[which.min(approx_df$Predicted),]$Predicted,
  #            size=.5,
  #            col="black",
  #            linetype="dashed") + 
  #   ylab("Y") +
  #   xlab("X")
  
  for (size_name in names(sizes)) {
    size = sizes[size_name]
    
    
    for (iter in 1:len(plist)) {

      
      png(sprintf("%s\\%s\\example_iter_%d.png",
                  write_path,
                  size_name,
                  iter),
          height=size,
          width=size,
          units="in",
          res=300)
      plot(plist[[iter]])
      dev.off()
            
      pdf(sprintf("%s\\%s\\example_iter_%d.pdf",
                  write_path,
                  size_name,
                  iter),
        height=size,
        width=size)
      plot(plist[[iter]])
      dev.off()
      
    }
    # pdf(sprintf("%s\\%s\\example_polynom.pdf",
    #             write_path,
    #             size_name),
    #   height=size,
    #   width=size)
    # plot(polynom_plot)
    # dev.off()

  }

}

supp_figure_bayesian_optimizer_param_space_minimum <- function() {
  write_path <- sprintf("%s\\supp_figure_simulations\\",figures_path)
  dir.create(write_path)
  
  write_path <- sprintf("%s\\param_space\\", write_path)
  dir.create(write_path)
  dir.create(sprintf("%s\\big\\",write_path))
  dir.create(sprintf("%s\\medium\\",write_path))
  dir.create(sprintf("%s\\medium_1\\",write_path))
  dir.create(sprintf("%s\\small\\", write_path))
  
  sizes=c(big=3,
          medium=2.5,
          medium_1=2.2,
          small=2)
  
  datasets_to_use <- 
                  list(list(dataset=4, param="1.000", session="7"),
                       list(dataset=3, param="1.000", session="8"),
                       list(dataset=5, param="1.000", session="7"),
                       list(dataset=6, param="1.000", session="8"),
                       list(dataset=5, param="1.000", session="13"),
                       list(dataset=6, param="1.000", session="14"),
                       list(dataset=4, param="1.000", session="7"),
                       list(dataset=3, param="1.000", session="15"),
                       list(dataset=4, param="1.000", session="16"),
                       list(dataset=6, param="1.000", session="6"),
                       list(dataset=5, param="1.000", session="16"),
                       list(dataset=3, param="1.000", session="14"),
                       list(dataset=1, param="0.600", session="3"),
                       list(dataset=2, param="0.700", session="4"),
                       list(dataset=1, param="0.600", session="2"),
                       list(dataset=2, param="0.800", session="2"),
                       list(dataset=7, param="0.700", session="4"),
                       list(dataset=8, param="0.600", session="3"),
                       list(dataset=7, param="0.800", session="2"))
  
  for (dataset in datasets_to_use){
    load(sprintf("%s\\equalized\\session_%s\\KS_simulations_estimate\\param_fit_%s\\pred_df.R", 
                 all_data_paths[dataset$dataset],
                 dataset$session,
                 dataset$param))
    
    

    
    predicted <- t(predicted_df[,5:24])
    predicted <- (predicted - min(predicted)) / (max(predicted)-min(predicted))
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
              plot.margin = unit(c(0,0,0,0), "cm"),
              text=element_text(size=13)) 
      
      thme_leg <-  theme_light() +     
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              plot.margin = unit(c(0,0,0,0), "cm"),
              text=element_text(size=13)) 
      
      
      # ylab(names(values[ind[2]])) +
      #   xlab(names(values[ind[1]]))
      
      
      gmin <- 
        ggplot(melt_min, aes(x=x,y=y,fill=Cost)) + 
        geom_tile(color="NA", size=0) + 
        scale_fill_distiller(palette="Spectral") + 
        thme +
        ylab(names(values[ind[2]])) + xlab(names(values[ind[1]]))
      
      gmin_leg <- 
        ggplot(melt_min, aes(x=x,y=y,fill=Cost)) + 
        geom_tile(color="NA", size=0) + 
        scale_fill_distiller(palette="Spectral") + 
        thme_leg +
        ylab(names(values[ind[2]])) + xlab(names(values[ind[1]]))
      
      gmean <- 
        ggplot(melt_mean, aes(x=x,y=y,fill=Cost)) + 
        geom_tile(color="NA", size=0) + 
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
    mean_p <- do.call(plot_grid, mean_plots)
    min_p <- do.call(plot_grid, min_plots)
    min_leg_p <- do.call(plot_grid, min_plots_legend)
    
    
    for (size_name in names(sizes)) {
      size <- sizes[[size_name]]
      
      dir.create(sprintf("%s\\%s\\dataset_%s_param_%s_session_%s",
                         write_path,
                         size_name,
                         dataset$dataset,
                         dataset$param,
                         dataset$session))
      
      pdf(file=sprintf("%s\\%s\\dataset_%s_param_%s_session_%s\\params_min.pdf",
                       write_path,
                       size_name,
                       dataset$dataset,
                       dataset$param,
                       dataset$session),
          width=size * 2,
          height=size * 1.2)
      plot(min_p)
      dev.off()
      
      pdf(file=sprintf("%s\\%s\\dataset_%s_param_%s_session_%s\\params_mean.pdf",
                       write_path,
                       size_name,
                       dataset$dataset,
                       dataset$param,
                       dataset$session),
          width=size * 2,
          height=size * 1.2)
      plot(mean_p)
      dev.off()
      
      pdf(file=sprintf("%s\\%s\\dataset_%s_param_%s_session_%s\\params_min_legend.pdf",
                       write_path,
                       size_name,
                       dataset$dataset,
                       dataset$param,
                       dataset$session),
          width=size * 2.5,
          height=size * 1.2)
      plot(min_leg_p)
      dev.off()
    }
  }
}

supp_figure_best_fit <- function() {
  write_path <- sprintf("%s\\supp_figure_simulations\\",figures_path)
  dir.create(write_path)
  
  write_path <- sprintf("%s\\simulations_fit\\", write_path)
  dir.create(write_path)
  dir.create(sprintf("%s\\big\\",write_path))
  dir.create(sprintf("%s\\medium\\",write_path))
  dir.create(sprintf("%s\\medium_1\\",write_path))
  dir.create(sprintf("%s\\medium_2\\",write_path))
  dir.create(sprintf("%s\\small\\", write_path))
  
  sizes = c(big=2,
            medium=1.75,
            medium_1=1.5,
            medium_2=1.25,
            small=1.15)
  comp_df_all <- create_comparision_of_properties()
  write.csv(file=sprintf("%s\\comparision_df.Rda", write_path), comp_df_all)
  metrics <- unique(comp_df_all$metric)
  properties <- unique(comp_df_all$property)
  subfields <- unique(comp_df_all$Subfield)
  
  statistics_df <- data.frame()
  values_df <- data.frame()
  
  pooled_sessions_statistics_df <- data.frame()
  pooled_sessions_values_df <- data.frame()
  
  for (metric_to_use in metrics) {
    metric_df <- comp_df_all[comp_df$metric == metric_to_use,]
  
    statistical_plots <- list()
    violin_plots <- list()
    histograms <- list()
    density_plots <- list()
    compare_boxplots <- list()
    
    for (subfield in subfields) {
      subfield_df <- metric_df[metric_df$Subfield == subfield,]
      
      
      for (property in properties) {
        property_df <- subfield_df[subfield_df$property == property,]
        property_df$param <- unlist(lapply(str_split(property_df$param, " "), function(lst) {lst[[1]]}))
        
        property_df$value <- (property_df$value - min(property_df$value))/
                              (max(property_df$value) - min(property_df$value))
        statistical_test <- ks.test(property_df[property_df$param != "Chosen",]$value,
                                    property_df[property_df$param == "Chosen",]$value,
                                    alternative = "less")
       
        paired_wilcox_df <- data.frame()
        
        statistics_df <- rbind(statistics_df, 
                               data.frame(property=property,
                                          metric=metric_to_use,
                                          subfield=subfield,
                                          pvalue=statistical_test$p.value,
                                          statistic=statistical_test$statistic,
                                          alternative=statistical_test$alternative,
                                          method=statistical_test$method,
                                          signif_code=signif.num(statistical_test$p.value)))
        
        values_df <- rbind(values_df, 
                               data.frame(property=property,
                                          metric=metric_to_use,
                                          subfield=subfield,
                                          mean_opt=mean(property_df[property_df$param == "Chosen",]$value),
                                          mean_sampled=mean(property_df[property_df$param != "Chosen",]$value),
                                          sd_opt=sd(property_df[property_df$param == "Chosen",]$value),
                                          sd_sampled=sd(property_df[property_df$param != "Chosen",]$value),
                                          sem_opt=sem(property_df[property_df$param == "Chosen",]$value),
                                          sem_sampled=sem(property_df[property_df$param != "Chosen",]$value)))
        gstat <- 
          ggplot()  +
          ggtitle(label=statistical_test$statistic,
                  subtitle=signif.num(statistical_test$p.value))
        ghist <-
        ggplot(property_df, aes(x=value)) +
          geom_histogram(aes(group=param, fill=param), bins=100, alpha=0.8) +
          theme_classic() + 
          base_plot_theme + 
          ggtitle(sprintf("%s - %s", subfield, property)) +
          theme(text=element_text(size=13.5),
                plot.title=element_text(size=11),
                plot.margin=unit(c(0,0,0,0), "cm")) +
          xlab(metric_to_use)  +
          scale_fill_manual(values=c("gray20", "#F05A28")) +
          ylab("Count")
        
        gdens <-
          ggplot(property_df, aes(x=value)) +
          geom_density(aes(group=param, fill=param), color=NA, alpha=0.3) +
          theme_classic() + 
          base_plot_theme + 
          ggtitle(sprintf("%s - %s", subfield, property)) +
          theme(text=element_text(size=13.5),
                plot.title=element_text(size=11),
                plot.margin=unit(c(0,0,0,0), "cm")) +
          xlab(metric_to_use) +
          scale_fill_manual(values=c("gray20", "#F05A28")) +
          ylab("Density")
        
        gviolin <- 
          ggplot(property_df, aes(x=param, y=value)) +
          geom_violin(aes(group=param, fill=param, alpha=.3), size=.75) +
          theme_classic() + 
          stat_summary(fun=mean, geom="point", size=1.5, color="black", aes(x=param, y=value)) + 
          base_plot_theme + 
          ggtitle(sprintf("%s - %s", subfield, property)) +
          theme(text=element_text(size=13.5),
                plot.title=element_text(size=11),
                plot.margin=unit(c(0,0,0,0), "cm")) +
          ylab(metric_to_use) +
          scale_fill_manual(values=c("gray20", "#F05A28")) + 
          scale_x_discrete(breaks=c("Chosen", "Param"), 
                           labels=c("Optimized", "Sampled")) + 
          xlab("")
        
        
        statistical_plots <- append(statistical_plots, list(gstat))
        violin_plots <- append(violin_plots, list(gviolin))
        histograms <- append(histograms, list(ghist))
        density_plots <- append(density_plots, list(gdens))
        
        paired_wilcox_df <- data.frame()
        gpaired = ggplot() +
        base_plot_theme + 
          ggtitle(sprintf("%s - %s", subfield, property)) +
          theme(text=element_text(size=13.5),
                plot.title=element_text(size=11),
                plot.margin=unit(c(0,0,0,0), "cm")) + 
          ylab(metric_to_use) +
          xlab("")
        
        
        
        for (session in unique(property_df$Session)) {
        #for (mice in unique(property_df$Mice)) {
          #for (direc in unique(property_df$Direction)) {
            
              
              wilc_df <- property_df[property_df$Session == session,]
                                     
              
              
               
                chosen <- wilc_df[wilc_df$param == "Chosen",]
                
                param <- wilc_df[wilc_df$param == "Param",]
                
                paired_wilcox_df <- rbind(paired_wilcox_df,
                                          data.frame(Random=mean(param$value),
                                                     Chosen=mean(chosen$value),
                                                     Session=session))
                
                gpaired <- gpaired +
                  geom_line(data=data.frame(x=c("Optimized", "Sampled"),
                                            y=c(mean(chosen$value), mean(param$value))),
                            aes(x=x,y=y, group=1), alpha=.3) +
                  geom_point(data=data.frame(x=c("Optimized", "Sampled"),
                                             y=c(mean(chosen$value), mean(param$value))),
                             aes(x=x,y=y, group=1))

            }
        
        wilk_test <- wilcox.test(paired_wilcox_df$Chosen, paired_wilcox_df$Random, paired=T, alternative="less") 
        
        max_val <- max(paired_wilcox_df[,1:2])
        
        
        signif_df = data.frame(y=max_val, x=1.5, label=signif.num(wilk_test$p.value))
        
        gpaired <- 
        gpaired + 
          geom_text(data=signif_df, aes(x=x,y=y,label=label))
        
        
        
        compare_boxplots <- append(compare_boxplots, list(gpaired))
        
        pooled_sessions_statistics_df <- 
          rbind(pooled_sessions_statistics_df, 
                               data.frame(property=property,
                                          metric=metric_to_use,
                                          subfield=subfield,
                                          pvalue=wilk_test$p.value,
                                          statistic=wilk_test$statistic,
                                          alternative=wilk_test$alternative,
                                          method=wilk_test$method,
                                          signif_code=signif.num(wilk_test$p.value)))
        
        pooled_sessions_values_df <- rbind(pooled_sessions_values_df, 
                           data.frame(property=property,
                                      metric=metric_to_use,
                                      subfield=subfield,
                                      mean_opt=mean(paired_wilcox_df$Chosen),
                                      mean_sampled=mean(paired_wilcox_df$Random),
                                      sd_opt=sd(paired_wilcox_df$Chosen),
                                      sd_sampled=sd(paired_wilcox_df$Random),
                                      sem_opt=sem(paired_wilcox_df$Chosen),
                                      sem_sampled=sem(paired_wilcox_df$Random),
                                      Nsessions=len(paired_wilcox_df$Session)))
      }
    }
    
    statistical_plots$nrow <- 2
    statistical_plots$align <- "v"
    violin_plots$nrow <- 2
    violin_plots$align <- "v"
    histograms$nrow <- 2
    histograms$align <- "v"
    density_plots$nrow <- 2
    density_plots$align <- "v"
    compare_boxplots$nrow <- 2
    compare_boxplots$align <- "v"
    
    
    gstat_all <- do.call(plot_grid, statistical_plots)
    gviol_all <- do.call(plot_grid, violin_plots)
    ghist_all <- do.call(plot_grid, histograms)
    gdens_all <- do.call(plot_grid, density_plots)
    gcomp_all <- do.call(plot_grid, compare_boxplots)
    
    for (size_name in names(sizes)) {
      size = sizes[[size_name]]
      
      dir.create(sprintf("%s\\%s\\%s",
                 write_path,
                 size_name,
                 metric_to_use))
      
      pdf(file=sprintf("%s\\%s\\%s\\violin.pdf",
                       write_path,
                       size_name,
                       metric_to_use),
          height=size*2,
          width=size*4)
      plot(gviol_all)
      dev.off()
      
      pdf(file=sprintf("%s\\%s\\%s\\hists.pdf",
                       write_path,
                       size_name,
                       metric_to_use),
          height=size*2,
          width=size*4)
      plot(ghist_all)
      dev.off()      
      
      pdf(file=sprintf("%s\\%s\\%s\\dens.pdf",
                       write_path,
                       size_name,
                       metric_to_use),
          height=size*2,
          width=size*4)
      plot(gdens_all)
      dev.off()
      
      pdf(file=sprintf("%s\\%s\\%s\\stats.pdf",
                       write_path,
                       size_name,
                       metric_to_use),
          height=size*2,
          width=size*4)
      plot(gstat_all)
      dev.off()      
      
      pdf(file=sprintf("%s\\%s\\%s\\compared_lineplots.pdf",
                       write_path,
                       size_name,
                       metric_to_use),
          height=size*2,
          width=size*4)
      plot(gcomp_all)
      dev.off()      
    }
  }
  
  write.csv(pooled_sessions_values_df, file=sprintf("%s//session_values.csv", write_path))
  write.csv(pooled_sessions_statistics_df, file=sprintf("%s//session_statistics.csv", write_path))
  write.csv(values_df, file=sprintf("%s//dist_all_values.csv", write_path))
  write.csv(statistics_df, file=sprintf("%s//dist_all_statistics.csv", write_path))
}

supp_figure_simulated_rasters <- function() {
  write_path <- sprintf("%s\\supp_figure_simulations\\",figures_path)
  dir.create(write_path)
  
  write_path <- sprintf("%s\\simulated_rasters\\", write_path)
  dir.create(write_path)
  
  for (ext in c("\\fit_simulation")) {
    for (p in all_data_paths[c(3,4,5)]) {
      for (session in c(1,2,7,13)) {
        
        a <- get_spike_train_and_stim_trace_from_path(p,session)
        spike_train <- a[[1]]
        stim_trace <- a[[2]]
        
        fr <- rowMeans(spike_train) / dt
        processed_real <- preprocess_spike_train(spike_train, stim_trace)
          
        true_cells_spike_train <- processed_real$working_cells_spike_train
        true_firing_rate <- processed_real$working_firing_rate
        true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
          
        params <- get_fit_params(sprintf("%s\\equalized\\session_%d", p, session))
        print(params)
        simulated_tuning_curve <-
             generate_tuning_curves_cost(n = nrow(spike_train),
                                        percentage = params$pct,
                                        average_width = params$params["average_width"], 
                                        sd_width = params$params["sd_width"],
                                        fixed_fr=fr,
                                        noise=params$params["noise"],
                                        double_peak_pct = params$params["double_peak_pct"],
                                        n_bins=24,
                                        plot=F)
          
        pois_factor <- currate_spike_train_cost(simulated_tuning_curve, 
                                                  true_time_bins_per_cells,
                                                  stim_trace,
                                                  verbose = T,
                                                  jump=0.05)
          
        generated_spike_train <- 
            generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                                       stim_trace = stim_trace,
                                       factor=pois_factor,
                                       fs=1)
          
        pval_by_time = pct_by_sample_size_figures(stim_trace = stim_trace,
                                                    spike_train=generated_spike_train)
        spike_train <- generated_spike_train
        
        decreasing_rand <- apply(t(apply(pval_by_time$rand, 1, diff)), 1, function(r) {sum(r < 0)})
        
        session_length = 20
        
        
        run_df <- data.frame(Time=((1:len(stim_trace) / len(stim_trace)) * session_length), 
                             Position=stim_trace * 4) 
        run_df$Position[which(run_df$Position == 4)] <- 0
        actual_times <- (pval_by_time$duration / (ncol(spike_train) * dt)) * session_length
        
        
        if (grepl("CA3", p)) {
          mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', p))
          mice_str <- substr(p, mice_str_index, mice_str_index+4) 
        } else {
          mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', p))
          mice_str <- substr(p, mice_str_index, mice_str_index+3) 
        }
        
        
        dir <- ifelse(grepl("Left", p), "L", "R")
        issim <- ifelse(ext == "", "", "_simulation")
        mice_path <- sprintf("%s\\mice_%s_%s_%d%s",
                             write_path,
                             mice_str,
                             dir,
                             session,
                             issim)
        dir.create(mice_path)
        
        save(file=sprintf("%s\\pval_by_time_df.Rda", mice_path),
             pval_by_time)
        
        pcs <- get_place_cells(spike_train, stim_trace)
        group_indices <- list()
        group_indices$random_non_pcs <- pcs$ind[pcs$random_non_pcs]
        group_indices$random_pcs <- pcs$ind[pcs$random_pcs]
        group_indices$cyclic_pcs <- pcs$ind[pcs$cyclic_pcs]
        group_indices$cyclic_non_pcs <- pcs$ind[pcs$cyclic_non_pcs]
        group_indices$switchers <- pcs$ind[pcs$switchers]
        groups <- names(group_indices)
        
        for (group in groups) {
          mice_group_path <- sprintf("%s\\%s", mice_path, group)
          dir.create(mice_group_path)
          dir.create(sprintf("%s\\regular_size", mice_group_path))
          dir.create(sprintf("%s\\big_size", mice_group_path))
          dir.create(sprintf("%s\\med_size", mice_group_path))
          
          for (idx in group_indices[[group]]) {
            
            #idx <- round(runif(1,1,nrow(spike_train)))
            firing_ind <- which(spike_train[idx,] > 0)
            
            cell_df <- data.frame(Positions=stim_trace[firing_ind] * 4,
                                  Times=run_df$Time[firing_ind])
            
            pval_df_all <- data.frame(Pvalue=c(pval_by_time$rand[idx,], pval_by_time$cyclic[idx,]),
                                      Time=c(actual_times,actual_times),
                                      `Shuffle type`=c(rep("Random", times=length(actual_times)),
                                                       rep("Cyclic", times=length(actual_times))))
            
            pval_df <- data.frame(Pvalue=c(pval_by_time$cyclic[idx,]),
                                  Time=c(actual_times),
                                  `Shuffle type`=c(rep("Cyclic", times=length(actual_times))))
            
            tuning_curve <- (unlist(lapply(sort(unique(stim_trace)), function (sb) {mean(spike_train[idx,which(stim_trace == sb)]) / dt}))) 
            barp <- barplot(tuning_curve)
            smoothed_tuning <- smooth.spline(barp, tuning_curve, all.knots = T, lambda=1e-4)
            smoothed_tuning$y[smoothed_tuning$y < 0] <- 0
            pos <- sort(unique(stim_trace)) * 4; pos[pos == 4] <- 0
            smoothed_df <- data.frame(FR=smoothed_tuning$y,
                                      Positions=pos)
            
            if (len(simulated_tuning_curve) != 0) {
              sim_barp <- barplot(simulated_tuning_curve[idx,])
              sim_smoothed_tuning <- smooth.spline(sim_barp, simulated_tuning_curve[idx,], all.knots = T, lambda=1e-4)
              sim_smoothed_tuning$y[sim_smoothed_tuning$y < 0] <- 0
              sim_smoothed_df <- data.frame(FR=sim_smoothed_tuning$y / (dt * pois_factor),
                                            Positions=pos)
            }
            
            polygon_bin <- which(pval_df$Pvalue[2:30] < 0.05)[1]
            if(len(polygon_bin) > 0) {
              polygon_df <- polygon_df <- data.frame(x=c(pval_df$Time[polygon_bin],# - (pval_df$Time[polygon_bin] - pval_df$Time[polygon_bin-1]) * 0.5,
                                                         20, 20,
                                                         pval_df$Time[polygon_bin]),# -(pval_df$Time[polygon_bin] - pval_df$Time[polygon_bin-1]) * 0.5),
                                                     y=c(0,0,96,96))
              signif_text_pos = ifelse(polygon_bin - 5 > 0, polygon_bin -5, polygon_bin +5)
            } else {
              polygon_df(x=c(20,20,20,20), y=c(0,0,96,96))
              signif_text_pos = 4
            }
            
            grundf <-
              ggplot(polygon_df, aes(x=x,y=y)) + 
              # geom_polygon(fill=adjustcolor("deepskyblue1", alpha=0.2)) +
              # #ggplot(run_df) +
              # geom_line(data=run_df, aes(x=Time, y=Position),
              #           col="gray30") + 
              theme_light() +  
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    legend.position="NA",
                    #panel.border = element_blank(),
                    panel.background = element_blank()) +
              xlim(20,0) +#plot.margin = margin(t=0,r=0,b=0,l=0,unit="cm")) +
              ylim(0,96) +
              ylab("Position (cm)") +
              xlab("Time (minutes)") +
              geom_vline(xintercept = polygon_df$x[1],
                         linetype="dashed",
                         color="black") +
              geom_point(data=cell_df, aes(x=Times,y=Positions),
                         fill="#DA1C5C",
                         color="#DA1C5C",
                         alpha=.8,
                         stroke=0,
                         size=2) + 
              coord_flip()
            
            
            gtuning_df <- ggplot(smoothed_df) + 
              geom_line(aes(x=Positions, y=FR), col="red", size=2) +
              theme_light() +
              ylab("Firing rate (spikes/sec)") + 
              xlab("Position (cm)") +
              xlim(0,96) +
              ylim(0,ifelse(max(smoothed_df$FR) > 1.5, max(smoothed_df$FR), 1.5)) +
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    legend.position="NA",
                    panel.border = element_blank(),
                    panel.background = element_blank()) 
            
            if (len(simulated_tuning_curve) != 0) { 
              
              scale_factor = 0.5
              combined_tuning_df <- rbind(sim_smoothed_df, smoothed_df)
              combined_tuning_df$Type <- rep(c("Defined tuning curve", "Measured"), each=nrow(smoothed_df))
              
              gtuning_df <- ggplot(combined_tuning_df) + 
                geom_line(aes(x=Positions, y=FR, group=Type, color=Type), size=2, alpha=0.8) +
                theme_light() +
                ylab("Firing rate (spikes/sec)") + 
                xlab("Position (cm)") +
                scale_color_manual(breaks=c("Defined tuning curve", "Measured"), values=c("royalblue4", "red")) +
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      legend.position="NA",
                      panel.border = element_blank(),
                      panel.background = element_blank()) +
                ylim(c(0,ifelse(max(smoothed_df$FR) > 1.5, max(smoothed_df$FR), 1.5))) +
                xlim(0,96)
              
            }
            
            gpvaldf <- ggplot(pval_df) + 
              geom_line(aes(x=Time, y=Pvalue, color=`Shuffle.type`)) + 
              geom_point(aes(x=Time, y=Pvalue, color=`Shuffle.type`)) + 
              theme_light() +  
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    legend.position="NA",
                    panel.border = element_blank(),
                    panel.background = element_blank()) + 
              ylab(TeX(r'($P_{value}$)')) +
              ylim(0,1) + 
              xlim(0, 20) +
              geom_hline(yintercept = 0.05, linetype="dashed", col="gray60") +
              scale_color_manual(name = "Shuffle type", 
                                 breaks = c("Cyclic", "Random"),
                                 values = c("mediumorchid2", "lightsalmon2"))
            
            thm_big_text <- theme(text=element_text(size=15))
            thm_med_text <- theme(text=element_text(size=14))
            
            grundf_bigger_text <- grundf + thm_big_text
            gtuning_df_bigger_text <- gtuning_df + thm_big_text
            gpvaldf_bigger_text <- gpvaldf + thm_big_text
            
            grundf_med_text <- grundf + thm_med_text
            gtuning_df_med_text <- gtuning_df + thm_med_text
            gpvaldf_med_text <- gpvaldf + thm_med_text
            
            res <- align_plots(grundf, 
                               gtuning_df, 
                               gpvaldf, align="v")
            
            res_big_text <- align_plots(grundf_bigger_text, 
                                        gtuning_df_bigger_text, 
                                        gpvaldf_bigger_text, 
                                        align="v")
            
            res_med_text <- align_plots(grundf_med_text, 
                                        gtuning_df_med_text, 
                                        gpvaldf_med_text, 
                                        align="v")
            
            g <- grid.arrange(res[[1]], res[[2]], res[[3]], heights=c(1.8,1,1))
            g_big <- grid.arrange(res_big_text[[1]], 
                                  res_big_text[[2]], 
                                  res_big_text[[3]], 
                                  heights=c(1.8,1,1))
            
            g_med <- grid.arrange(res_med_text[[1]], 
                                  res_med_text[[2]], 
                                  res_med_text[[3]], 
                                  heights=c(1.8,1,1))
            
            
            png(sprintf("%s\\regular_size\\neur_%d.png", mice_group_path, idx), unit="in", height=5, width=1.5, res=300); 
            plot(g); dev.off()   
            
            pdf(sprintf("%s\\regular_size\\neur_%d.pdf", mice_group_path, idx), height=5, width=1.5); 
            plot(g); dev.off()   
            
            pdf(sprintf("%s\\big_size\\neur_%d.pdf", mice_group_path, idx), height=5, width=1.5); 
            plot(g_big); dev.off()   
            
            pdf(sprintf("%s\\med_size\\neur_%d.pdf", mice_group_path, idx), height=5, width=1.5); 
            plot(g_med); dev.off()
            
          }
        }
      }
    }
  }
}

create_comparision_of_properties <- function() {
  stimulus_trace_filename = "\\stim_trace.mat"
  spike_train_filename = "\\spike_train.mat"
  
  sessions_to_use <- 1:16
  equalized_paths <- sprintf("%s\\%s", all_data_paths, "equalized")
  comp_df_all <- data.frame()
  
  for (path_idx in 1:len(all_data_paths)) {
  
    working_path <- all_data_paths[path_idx]
    equalized_path <- equalized_paths[path_idx]
    
    m <- readMat(paste(working_path, stimulus_trace_filename, sep=""))
    
    if(len(grep("Left", working_path)) > 0) {
      stim_trace <- m$position.left.trials
    } else {
      stim_trace <- m$position.right.trials
    }
    
    
    m <- readMat(paste(working_path, spike_train_filename, sep=""))
    
    if(len(grep("Left", working_path)) > 0) {
      spike_train <- m$spikes.left.trials
    } else {
      spike_train <- m$spikes.right.trials
    }
    
    
    session_paths <- 
      list.dirs(equalized_path, recursive=F)[grep("session", 
                                                  list.dirs(equalized_path, 
                                                            recursive=F))]    
    
    tmp_spike_train <- spike_train
    tmp_stimulus_trace <- stim_trace
    
    # Run through all paths
    for (tpath in session_paths) {
      for (jx in 1:2) { 
    
        
        ses_idx <- sapply(str_split(tpath, "session_"), function(l) {as.numeric(l[2])})
        
        # Ignore session
        if (sum(ses_idx == sessions_to_use) == 0){
          next
        }
        
        spike_train <- tmp_spike_train
        stim_trace <- tmp_stimulus_trace
        spike_train <- spike_train[[ses_idx]][[1]]
        spike_train <- t(as.matrix(spike_train))
        stim_trace <- stim_trace[[ses_idx]][[1]]
        stim_trace <-  as.vector(stim_trace)
        
        tmp <- equalize_prior(spike_train, stim_trace)
        spike_train <- tmp$spike_train
        stim_trace <- tmp$stim_trace
        
        comp_df <- get_cost_compare_df_by_path(tpath, spike_train, stim_trace)
        
        if (grepl("CA3", tpath)) {
          mice_str_index <- unlist(gregexpr('C[0-9]{2}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+4) 
        } else {
          mice_str_index <- unlist(gregexpr('C[0-9]{1}M[0-9]', tpath))
          mice_str <- substr(tpath, mice_str_index, mice_str_index+3)
        }
        
        dir <- ifelse(grepl("Right", tpath), "R", "L")
        subfield <- ifelse(grepl("CA3", tpath), "CA3", "CA1")
        
        metadata_df <- matrix(rep(c(dir, ses_idx, mice_str, subfield, jx), each=nrow(comp_df)), nrow=nrow(comp_df))
        colnames(metadata_df) <- c("Direction", "Session", "Mice", "Subfield", "Repetition")
        comp_df <- cbind(comp_df,
                         metadata_df)
        
        comp_df_all <- rbind(comp_df_all,
                             comp_df)
      }
    }
  }
  
  return(comp_df_all)
}

get_cost_compare_df_by_path <- function(path, true_spike_train, stim_trace) {
  
  results <- data.frame()
  param_space <- list(double_peak_pct=seq(0, 0.5,length.out=7),
                      average_width=seq(0.01, 0.17, by=0.01),
                      noise=seq(0.01, 0.08, by=0.01),
                      sd_width=seq(0.0001, 0.015, length.out=10))
  
  parameter_set <- as.matrix(cross_df(param_space))
  
  metrics_functions <- list(wasserstein=wasserstein_dist,
                            JSD=JSD_dist,
                            KL=KL_dist,
                            dprime=dprime_dist)
  
  firing_rate <-rowMeans(true_spike_train) / dt
  processed_real <- preprocess_spike_train(true_spike_train, stim_trace=stim_trace, verbose=F)
  true_cells_spike_train <- processed_real$working_cells_spike_train
  true_firing_rate <- processed_real$working_firing_rate
  true_time_bins_per_cells <- processed_real$working_time_bins_per_cells
  
  tmp <- compute_tuning(true_cells_spike_train, stim_trace)
  stim_prob <- tmp[[1]]
  true_tuning_curve <- tmp[[2]]
  rm(tmp)
  
  tuning_for_peaks <- true_tuning_curve
  
  
  true_peaks <- unlist(apply(tuning_for_peaks, 1, 
                             function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
  true_active_spatial_bins_per_cell <- apply(true_tuning_curve, 1, function(n){sum(n>0)})
  #print(true_active_spatial_bins_per_cell)
  true_SI <- compute_SI(stim_prob, true_tuning_curve, true_firing_rate)
  
  
  
  params <-   
    extract_simulation_params(path, 
                              plot=F, 
                              absolute_min=T, 
                              estimation_path="KS_simulations_estimate", 
                              pct_range = seq(0.5, 1, by=0.1))
  
  
  for (param_idx in 1:nrow(params)) {
    
    params_to_use <- params[param_idx,]
    pct <- as.numeric(rownames(params)[param_idx])
    simulated_tuning_curve <-
      generate_tuning_curves_cost(n = nrow(true_spike_train),
                                  percentage = pct,
                                  average_width = params_to_use["average_width"], 
                                  sd_width = params_to_use["sd_width"],
                                  fixed_fr=firing_rate,
                                  noise=params_to_use["noise"],
                                  double_peak_pct = params_to_use["double_peak_pct"],
                                  plot=F)
    pois_factor <- currate_spike_train_cost(simulated_tuning_curve, 
                                            true_time_bins_per_cells,
                                            stim_trace,
                                            verbose = F,
                                            jump = 0.1)
    
    for (i in 1:8) {
      simulated_tuning_curve <-
        generate_tuning_curves_cost(n = nrow(true_spike_train),
                                    percentage = pct,
                                    average_width = params_to_use["average_width"], 
                                    sd_width = params_to_use["sd_width"],
                                    fixed_fr=firing_rate,
                                    noise=params_to_use["noise"],
                                    double_peak_pct = params_to_use["double_peak_pct"],
                                    plot=F)
      
      generated_spike_train <- 
        generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                                   stim_trace = stim_trace,
                                   factor=pois_factor,
                                   fs=1)
      
      processed_generated <- preprocess_spike_train(generated_spike_train, stim_trace, verbose=F)
      generated_cells_active_st <- processed_generated$working_cells_spike_train
      generated_firing_rate <- processed_generated$working_firing_rate
      generated_time_bins <- processed_generated$working_time_bins_per_cells
      
      tmp <- compute_tuning(generated_cells_active_st, stim_trace)
      gen_stim_prob <- tmp[[1]]
      gen_tuning_curve <- tmp[[2]]
      rm(tmp)
      
      generated_peaks <- unlist(apply(gen_tuning_curve, 1, 
                                      function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
      generated_active_spatial_bins_per_cell <- apply(gen_tuning_curve, 1, function(n){sum(n>0)})
      generated_SI <- compute_SI(gen_stim_prob, gen_tuning_curve, generated_firing_rate)
      
      for (metric_name in names(metrics_functions)) {
        metric <- metrics_functions[[metric_name]] 
        
        results <- rbind(results,
                         data.frame(metric=metric_name,param=sprintf("Chosen %d", i),
                                    pct=pct,property="SI",
                                    value=metric(true_SI[[1]], generated_SI[[1]])))
        
        
        results <- rbind(results,
                         data.frame(metric=metric_name,param=sprintf("Chosen %d", i),
                                    pct=pct,property="Peaks",
                                    value=metric(true_peaks, generated_peaks)))
        
        results <- rbind(results,
                         data.frame(metric=metric_name,param=sprintf("Chosen %d", i),
                                    pct=pct,property="Sbins",
                                    value=metric(true_active_spatial_bins_per_cell, 
                                                 generated_active_spatial_bins_per_cell)))
        results <- rbind(results,
                         data.frame(metric=metric_name,param=sprintf("Chosen %d", i),
                                    pct=pct,property="Time bins",
                                    value=metric(true_time_bins_per_cells, 
                                                 generated_time_bins)))   
        
        
      }
    }
    
    params_indices <- sample(1:nrow(parameter_set),20)
    for (idx in 1:len(params_indices)) { 
      params_to_use <- parameter_set[params_indices[idx],]
      simulated_tuning_curve <-
        generate_tuning_curves_cost(n = nrow(true_spike_train),
                                    percentage = pct,
                                    average_width = params_to_use["average_width"], 
                                    sd_width = params_to_use["sd_width"],
                                    fixed_fr=firing_rate,
                                    noise=params_to_use["noise"],
                                    double_peak_pct = params_to_use["double_peak_pct"],
                                    plot=F)
      
      generated_spike_train <- 
        generate_spike_trains_cost(tuning_curves = simulated_tuning_curve,
                                   stim_trace = stim_trace,
                                   factor=pois_factor,
                                   fs=1)
      
      
      
      processed_generated <- preprocess_spike_train(generated_spike_train, stim_trace, verbose=F)
      generated_cells_active_st <- processed_generated$working_cells_spike_train
      generated_firing_rate <- processed_generated$working_firing_rate
      generated_time_bins <- processed_generated$working_time_bins_per_cells
      
      tmp <- compute_tuning(generated_cells_active_st, stim_trace)
      gen_stim_prob <- tmp[[1]]
      gen_tuning_curve <- tmp[[2]]
      rm(tmp)
      
      generated_peaks <- unlist(apply(gen_tuning_curve, 1, 
                                      function(n){return(len(get_peaks(n, threshold_from_max = peaks_significance_threshold)))}))
      generated_active_spatial_bins_per_cell <- apply(gen_tuning_curve, 1, function(n){sum(n>0)})
      generated_SI <- compute_SI(gen_stim_prob, gen_tuning_curve, generated_firing_rate)
      
      for (metric_name in names(metrics_functions)) {
        metric <- metrics_functions[[metric_name]] 
        
        results <- rbind(results,
                         data.frame(metric=metric_name,param=sprintf("Param %d", idx),
                                    pct=pct,property="SI",
                                    value=metric(true_SI[[1]], generated_SI[[1]])))
        
        
        results <- rbind(results,
                         data.frame(metric=metric_name,param=sprintf("Param %d", idx),
                                    pct=pct,property="Peaks",
                                    value=metric(true_peaks, generated_peaks)))
        
        results <- rbind(results,
                         data.frame(metric=metric_name,param=sprintf("Param %d", idx),
                                    pct=pct,property="Sbins",
                                    value=metric(true_active_spatial_bins_per_cell, 
                                                 generated_active_spatial_bins_per_cell)))
        results <- rbind(results,
                         data.frame(metric=metric_name,param=sprintf("Param %d", idx),
                                    pct=pct,property="Time bins",
                                    value=metric(true_time_bins_per_cells, 
                                                 generated_time_bins)))   
        
        
      }
    }
  }
  
  return(results)
}




supp_figure_bayesian_optimizer_param_space_video <- function() {
  write_path <- sprintf("%s\\supp_figure_simulations\\",figures_path)
  dir.create(write_path)
  
  write_path <- sprintf("%s\\param_space_video\\", write_path)
  dir.create(write_path)
  dir.create(sprintf("%s\\big\\",write_path))
  dir.create(sprintf("%s\\medium\\",write_path))
  dir.create(sprintf("%s\\medium_1\\",write_path))
  dir.create(sprintf("%s\\small\\", write_path))
  
  sizes=c(big=3,
          medium=2.5,
          medium_1=2.2,
          small=2)
  

  for (iter in 1:30) {
    
    pred_df <- 
    read.csv(file=sprintf("C:\\Users\\itayta\\Desktop\\allcells_figures\\supp_figure_simulations\\for_movie\\iter_%d.csv", 
                 iter))
    
    
    
    
    
    
    #predicted <- t(pred_df)
    #predicted <- (predicted - min(predicted)) / (max(predicted)-min(predicted))
    predicted <- as.matrix(pred_df[,-1])
    values <- list(`Double Peak`=seq(0, 0.5,length.out=7),
                        `Width (mu)`=seq(0.01, 0.17, by=0.01),
                        `Noise`=seq(0.01, 0.08, by=0.01),
                        `Width (sigma)`=seq(0.0001, 0.015, length.out=10))
    
    print(param_space)
    
    X_star <- as.matrix(cross_df(values))
    
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
              plot.margin = unit(c(0,0,0,0), "cm"),
              text=element_text(size=13)) 
      
      thme_leg <-  theme_light() +     
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              plot.margin = unit(c(0,0,0,0), "cm"),
              text=element_text(size=13)) 
      
      
      # ylab(names(values[ind[2]])) +
      #   xlab(names(values[ind[1]]))
      
      
      gmin <- 
        ggplot(melt_min, aes(x=x,y=y,fill=Cost)) + 
        geom_tile(color="NA", size=0) + 
        scale_fill_distiller(palette="Spectral") + 
        thme +
        ylab(names(values[ind[2]])) + xlab(names(values[ind[1]]))
      
      gmin_leg <- 
        ggplot(melt_min, aes(x=x,y=y,fill=Cost)) + 
        geom_tile(color="NA", size=0) + 
        scale_fill_distiller(palette="Spectral") + 
        thme_leg +
        ylab(names(values[ind[2]])) + xlab(names(values[ind[1]]))
      
      gmean <- 
        ggplot(melt_mean, aes(x=x,y=y,fill=Cost)) + 
        geom_tile(color="NA", size=0) + 
        scale_fill_distiller(palette="Spectral") + 
        thme +
        ylab(names(values[ind[2]])) + xlab(names(values[ind[1]]))
      
      
      
      min_plots <- append(min_plots, list(gmin))
      mean_plots <- append(mean_plots, list(gmean))
      min_plots_legend <- append(min_plots_legend, list(gmin_leg))
    }
    
    min_plots$nrow =2
    mean_plots$nrow = 2
    #min_plots_legend$nrow = 2
    mean_p <- do.call(plot_grid, mean_plots)
    min_p <- do.call(plot_grid, min_plots)
    #min_leg_p <- do.call(plot_grid, min_plots_legend)
    
    
    for (size_name in names(sizes)) {
      size <- sizes[[size_name]]
      
      dir.create(sprintf("%s\\%s\\",
                         write_path,
                         size_name))
      
      png(file=sprintf("%s\\%s\\min_iter_%d.png",
                       write_path,
                       size_name,
                       iter),
          width=size * 2,
          height=size * 1.2,
          units="in",
          res=300)
      plot(min_p)
      dev.off()
      
      png(file=sprintf("%s\\%s\\mean_iter_%d.png",
                       write_path,
                       size_name,
                       iter),
          width=size * 2,
          height=size * 1.2,
          units="in",
          res=300)
      plot(mean_p)
      dev.off()
      
    
    }
  }
  
  
  png_files <- sprintf("%s\\%s\\min_iter_%d.png",
                       write_path,
                       "big",
                       1:30)
  
  av::av_encode_video(png_files, 'output.mp4', framerate = 3)
  
  utils::browseURL('output.mp4')
}