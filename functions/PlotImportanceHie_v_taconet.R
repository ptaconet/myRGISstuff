PlotImportanceHie_v_taconet<-function (input.data, X.data = 2, Y.data = 3, imp.data = 4, explanatory.variable = 5, plot.type = "Tile", 
          imp.title = colnames(input.data)[imp.data], X.Title = colnames(input.data)[X.data], 
          Y.Title = colnames(input.data)[Y.data], low.col = "blue", 
          high.col = "red", geom.tile.bor.col = "white", pos.col = "green", 
          zero.col = "white", neg.col = "red", supp.warn = TRUE, ...) 
{
  X.var <- Y.var <- NULL
  localenv <- environment()
  if (plot.type != "Tile" && plot.type != "Bubble") {
    stop("\nplot.type can be either 'Tile' or 'Bubble'\n")
  }
  work.data <- input.data[, c(X.data, Y.data, imp.data, explanatory.variable)]
  colnames(work.data) <- c("X.var", "Y.var", "imp.var", "type")
  # BUG IN THE SOURCE FUNCTION WITH THE FOLLOWING LINES : 
  #work.data$X.var <- suppressWarnings(factor(work.data$X.var, 
   #                                          as.character(work.data$X.var)))
  #work.data$Y.var <- suppressWarnings(factor(work.data$Y.var, 
  #                                           as.character(work.data$Y.var)))
  
  # PTACONET FIX BUG :
  work.data$X.var <- suppressWarnings(as.factor(as.character(work.data$X.var)))
  work.data$Y.var <- suppressWarnings(as.factor(as.character(work.data$Y.var)))
  work.data$type <- suppressWarnings(as.factor(as.character(work.data$type)))
  work.data=work.data[order(work.data$type),]
  
  if (plot.type == "Tile") {
    p <- ggplot(work.data, aes(X.var, Y.var), environment = localenv)
    p <- p + geom_tile(aes(fill = work.data$imp.var), color = geom.tile.bor.col) + facet_grid(type~.,scales="free_y",space="free_y")
    p <- p + xlab(X.Title)
    p <- p + ylab(Y.Title)
    p <- p + scale_fill_continuous(low = low.col, high = high.col, 
                                   guide = "colourbar", guide_colourbar(title = imp.title))
  }
  if (plot.type == "Bubble") {
    #for (k3 in 1:dim(work.data)[1]) {
    #  if (work.data[k3, 3] > 0) {
    #    work.data[k3, 4] <- "Positive"
    #  }
    #  if (work.data[k3, 3] == 0) {
    #    work.data[k3, 4] <- "Zero"
    #  }
    #  if (work.data[k3, 3] <= 0) {
    #    work.data[k3, 4] <- "Negative"
    #  }
    #}
    #colnames(work.data)[4] <- "sign.imp"
    #fillScaleValues <- c(Positive = pos.col, Zero = zero.col, 
    #                     Negative = neg.col)
    p <- ggplot(work.data, aes(X.var, Y.var), environment = localenv)
    p <- p + geom_point(aes(size = abs(work.data$imp.var), 
                            fill = work.data[,4]), shape = 21) + facet_grid(type~.,scales="free_y",space="free_y")
    p <- p + labs(y="var. used",fill = names(input.data[explanatory.variable]),title = paste0("Variables importance by classifier and by ",names(input.data[explanatory.variable]))) 
    p <- p + xlab(X.Title)
    #p <- p + ylab(Y.Title)
    p <- p + scale_size_continuous(guide = "legend", guide_legend(title = imp.title))
    #p <- p + scale_fill_manual(values = fillScaleValues, 
    #                           name = explanatory.variable.name)
  }
  if (supp.warn) {
    return(suppressWarnings(print(p)))
  }
  else {
    return(print(p))
  }
}
