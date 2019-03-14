PerformanceHRF<-function (hie.RF, per.index = c("flat.measures", "hie.F.measure"), 
    crisp.rule = c("stepwise.majority", "multiplicative.majority", 
        "multiplicative.permutation"), perm.num = 500, by.node = TRUE, 
    div.logical = TRUE, div.print = 25, beta.h.F = 1, ...) 
{
    cat("######################")
    cat(paste("\n", "-->  Veryfing Call", "\n", sep = ""))
    if (class(hie.RF) != "HRF") {
        stop(paste("\n", "In PerformanceHRF:  hie.RF should be of class HRF", 
            "\n", sep = ""))
    }
    temp.vec <- c("flat.measures", "hie.F.measure")
    if (length(intersect(per.index, temp.vec)) < 1) {
        stop(paste("\n", "In PerformanceHRF, at least one valid option for per.index is required", 
            "\n", sep = ""))
    }
    if (length(intersect(per.index, temp.vec)) > 0) {
        cat(paste("\n", "Performance will be evaluated using: ", 
            per.index, "\n", sep = ""))
    }
    rm(temp.vec)
    temp.vec <- c("multiplicative.majority", "multiplicative.permutation", 
        "stepwise.majority")
    if (length(intersect(crisp.rule, temp.vec)) < 1) {
        stop(paste("\n", "In PerformanceHRF, at least one valid option for crisp.rule is required", 
            "\n", sep = ""))
    }
    if (length(intersect(crisp.rule, temp.vec)) > 0) {
        cat(paste("\n", "Performance will be evaluated using crisp classification based on: ", 
            crisp.rule, "\n", sep = ""))
    }
    if (!is.numeric(perm.num)) {
        stop(paste("\n", "In PerformanceHRF, perm.num should be a positive integer", 
            "\n", sep = ""))
    }
    if (perm.num < 1) {
        stop(paste("\n", "In PerformanceHRF, perm.num should be a positive integer", 
            "\n", sep = ""))
    }
    if (round(perm.num, 0) != perm.num) {
        stop(paste("\n", "In PerformanceHRF, perm.num should be a positive integer", 
            "\n", sep = ""))
    }
    if (!is.logical(by.node)) {
        cat(paste("\n", "In PerformanceHRF, by.node should be Logical. default of by.node=TRUE is used", 
            "\n", sep = ""))
        by.node <- TRUE
    }
    if (!is.logical(div.logical)) {
        cat(paste("\n", "In PerformanceHRF, div.logical should be Logical. default of div.logical=TRUE is used", 
            "\n", sep = ""))
        div.logical <- TRUE
    }
    if (!is.numeric(div.print)) {
        cat(paste("\n", "In PerformanceHRF, div.print should be a positive integer", 
            "\n", "Default of div.print=25 is used", "\n", sep = ""))
        div.print <- 25
    }
    if (div.print < 1) {
        cat(paste("\n", "In PerformanceHRF, div.print should be a positive integer", 
            "\n", "Default of div.print=25 is used", "\n", sep = ""))
        div.print <- 25
    }
    if (round(div.print, 0) != div.print) {
        cat(paste("\n", "In PerformanceHRF, div.print should be a positive integer", 
            "\n", "Default of div.print=25 is used", "\n", sep = ""))
        div.print <- 25
    }
    if (!is.numeric(beta.h.F)) {
        cat(paste("\n", "In PerformanceHRF, beta.h.F should be a non-negative number", 
            "\n", "Default of beta.h.F=1 is used", "\n", sep = ""))
        beta.h.F <- 1
    }
    if (beta.h.F < 0) {
        cat(paste("\n", "In PerformanceHRF, beta.h.F should be a non-negative number", 
            "\n", "Default of beta.h.F=1 is used", "\n", sep = ""))
        beta.h.F <- 1
    }
    cat(paste("\n", "-->  Call Verified", "\n", sep = ""))
    if (!is.na(match("multiplicative.majority", crisp.rule))) {
        mult.maj.rule.logical <- TRUE
    }
    if (is.na(match("multiplicative.majority", crisp.rule))) {
        mult.maj.rule.logical <- FALSE
    }
    if (!is.na(match("multiplicative.permutation", crisp.rule))) {
        mult.perm.logical <- TRUE
    }
    if (is.na(match("multiplicative.permutation", crisp.rule))) {
        mult.perm.logical <- FALSE
    }
    if (!is.na(match("stepwise.majority", crisp.rule))) {
        step.maj.rule.logical <- TRUE
    }
    if (is.na(match("stepwise.majority", crisp.rule))) {
        step.maj.rule.logical <- FALSE
    }
    if (!is.na(match("flat.measures", per.index))) {
        flat.measures.logical <- TRUE
    }
    if (is.na(match("flat.measures", per.index))) {
        flat.measures.logical <- FALSE
    }
    if (!is.na(match("hie.F.measure", per.index))) {
        hie.F.measure.logical <- TRUE
    }
    if (is.na(match("hie.F.measure", per.index))) {
        hie.F.measure.logical <- FALSE
    }
    cat("######################")
    cat(paste("\n", "-->  Preparing the observed terminal nodes ", 
        "\n", sep = ""))
    cat(paste("\n", "Extracting info from hie.RF ", "\n", sep = ""))
    train.data.ready <- hie.RF$train.data.ready
    case.ID <- hie.RF$case.ID
    hie.levels <- hie.RF$hie.levels
    end.path.name <- hie.RF$call$end.path.name
    if (is.null(end.path.name)) {
        end.path.name <- "END.PATH"
    }
    unique.path <- hie.RF$hier.struc$unique.path
    train.data.acc <- subset(train.data.ready, select = case.ID)
    train.data.acc$obs.term.node <- NA
    cat(paste("\n", "Evaluating terminal node for each case ", 
        "\n", sep = ""))
    for (K2 in 1:nrow(train.data.ready)) {
        focal.path <- train.data.ready[K2, hie.levels]
        focal.term.node <- GetTerminalNode(end.path.name = end.path.name, 
            unique.path = focal.path, level.depth = ncol(focal.path) - 
                2)
        train.data.acc[K2, "obs.term.node"] <- focal.term.node$term.node.name
    }
    rm(list = c("K2", "focal.path", "focal.term.node"))
    cat(paste("\n", "Ordering according to case.ID ", "\n", sep = ""))
    train.data.acc <- train.data.acc[order(train.data.acc[, 1]), 
        ]
    unique.nodes <- unique(train.data.acc$obs.term.node)
    cat(paste("\n", "-->  Observed terminal nodes ready,  ", 
        "\n", sep = ""))
    cat("######################")
    cat(paste("\n", "-->  Createing the crisp.rule by per.index data frames", 
        "\n", sep = ""))
    cat(paste("\n", "Identifying indices called by user", "\n", 
        sep = ""))
    hie.performance <- data.frame(crisp.rule = c(NA))
    if (flat.measures.logical) {
        hie.performance$Accuracy <- NA
        hie.performance$Kappa <- NA
        hie.performance$AccuracyLower <- NA
        hie.performance$AccuracyUpper <- NA
        hie.performance$AccuracyNull <- NA
        hie.performance$AccuracyPValue <- NA
        hie.performance$McnemarPValue <- NA
    }
    if (hie.F.measure.logical) {
        hie.performance$h.precision <- NA
        hie.performance$h.recall <- NA
        hie.performance$h.F.measure <- NA
    }
    if (by.node) {
        if (flat.measures.logical) {
            #nodes.acc.ind <- c("Sensitivity", "Specificity", 
            #    "Pos Pred Value", "Neg Pred Value", "Prevalence", 
            #    "Detection Rate", "Detection Prevalence", "Balanced Accuracy")
            nodes.acc.ind <- c("Sensitivity", "Specificity","Precision", "Recall", "F1", "Pos Pred Value", "Neg Pred Value", "Prevalence", "Detection Rate", "Detection Prevalence", "Balanced Accuracy")
            nodes.measures <- expand.grid(unique.nodes, nodes.acc.ind)
            colnames(nodes.measures) <- c("term.node", "index")
            nodes.measures$term.node.index <- paste(nodes.measures$term.node, 
                nodes.measures$index, sep = ";")
            hie.performance[, nodes.measures$term.node.index] <- NA
            nodes.measures.columns <- nodes.measures
        }
        if (hie.F.measure.logical) {
            nodes.h.acc.ind <- c("n.h.precision", "n.h.recall", 
                "n.h.F.measure")
            nodes.h.measures <- expand.grid(unique.nodes, nodes.h.acc.ind)
            colnames(nodes.h.measures) <- c("term.node", "index")
            nodes.h.measures$term.node.index <- paste(nodes.h.measures$term.node, 
                nodes.h.measures$index, sep = ";")
            hie.performance[, nodes.h.measures$term.node.index] <- NA
            nodes.measures.columns <- nodes.h.measures
        }
        if (flat.measures.logical && hie.F.measure.logical) {
            nodes.measures.columns <- rbind(nodes.measures, nodes.h.measures)
        }
    }
    cat(paste("\n", "Starting output data frames for each crisp rule", 
        "\n", sep = ""))
    if (mult.maj.rule.logical) {
        hie.perf.mult.maj <- hie.performance
        hie.perf.mult.maj[1, "crisp.rule"] <- "multiplicative.majority.rule"
    }
    if (mult.perm.logical) {
        hie.perf.mult.perm <- hie.performance
        hie.perf.mult.perm[c(1:perm.num), "crisp.rule"] <- c(paste("perm.", 
            1:perm.num, sep = ""))
    }
    if (step.maj.rule.logical) {
        hie.perf.step.maj <- hie.performance
        hie.perf.step.maj[1, "crisp.rule"] <- "stepwise.majority"
    }
    cat(paste("\n", "-->  crisp.rule by per.index data frames created", 
        "\n", sep = ""))
    cat("######################")
    cat(paste("\n", "-->  Preparing the predicted terminal nodes from the hierarchical RandomForest ", 
        "\n", sep = ""))
    raw.votes <- predict.HRF(object = hie.RF)
    focal.votes <- raw.votes$prop.vote.train
    crisp.case.class <- focal.votes[, c(1, 2)]
    if (mult.maj.rule.logical || mult.perm.logical) {
        cat(paste("\n", "Estimating the Multiplicative probabilities down the hierarchical tree ", 
            "\n", sep = ""))
        multiplicative.prop.last.level <- GetMultPropVotes(prop.vote = focal.votes, 
            unique.path = unique.path, all.levels = FALSE)
        multiplicative.prop <- multiplicative.prop.last.level[[1]]
        if (mult.maj.rule.logical) {
            cat(paste("\n", "Adding the multiplicative majority Rule", 
                "\n", sep = ""))
            mult.maj.rule <- GetMultMajRule(prop.multiplicative.votes = multiplicative.prop, 
                bind.prop.mul.maj = FALSE)
            crisp.case.class <- cbind(crisp.case.class, mult.maj.rule[3])
            cat(paste("\n", "Multiplicative majority Rule added", 
                "\n", sep = ""))
        }
        if (mult.perm.logical) {
            cat(paste("\n", " Adding: ", perm.num, " permutations", 
                "\n", sep = ""))
            mult.perm.rule <- GetPermMultTermNode(multiplicative.prop.votes = multiplicative.prop, 
                perm.num = perm.num, div.logical = div.logical, 
                div.print = div.print, bind.prop.perm = FALSE)
            crisp.case.class <- cbind(crisp.case.class, mult.perm.rule[3:ncol(mult.perm.rule)])
            cat(paste("\n", perm.num, " Permutations added", 
                "\n", sep = ""))
        }
        multiplicative.prop <- multiplicative.prop[order(multiplicative.prop[, 
            2]), ]
    }
    if (step.maj.rule.logical) {
        cat(paste("\n", "Adding the stepwise majority Rule", 
            "\n", sep = ""))
        step.maj.rule <- GetStepMajRule(hie.RF = hie.RF, prop.vote = focal.votes, 
            bind.prop.step.maj = FALSE)
        crisp.case.class <- cbind(crisp.case.class, step.maj.rule[3])
        cat(paste("\n", "Stepwise majority Rule added", "\n", 
            sep = ""))
    }
    cat(paste("\n", "Ordering according to case.ID ", "\n", sep = ""))
    crisp.case.class <- crisp.case.class[order(crisp.case.class[, 
        2]), ]
    if (!all(train.data.acc[, 1] == crisp.case.class[, 2])) {
        stop(paste("\n", "Error in function: PerformanceHRF", 
            "\n", "case.ID names (rows) of the training data and predicted data do not match perfectly", 
            "\n", sep = ""))
    }
    cat(paste("\n", "-->  Predicted terminal nodes ready ", "\n", 
        sep = ""))
    cat("######################")
    cat(paste("\n", "-->  Start estimating performance measures ", 
        "\n", sep = ""))
    if (mult.maj.rule.logical) {
        cat(paste("\n", "Estimating performance measures for the multiplicative majority rule ", 
            "\n", sep = ""))
        focal.detail <- "multiplicative.majority.rule"
        pred.nodes <- crisp.case.class[, match(focal.detail, 
            colnames(crisp.case.class))]
        joined.levels <- JoinLevels(vector.1 = pred.nodes, vector.2 = train.data.acc$obs.term.node)
        pred.nodes <- joined.levels$vector.1
        obse.nodes <- joined.levels$vector.2
        rm(joined.levels)
        conf.matr <- confusionMatrix(data = pred.nodes, reference = obse.nodes, 
            dnn = c("Prediction", "Observed"))
        conf.matr.table <- conf.matr$table
        conf.matr.overall <- as.data.frame(t(conf.matr$overall))
        conf.matr.byclass <- conf.matr$byClass
        rownames(conf.matr.byclass) <- sub("Class: ", "", rownames(conf.matr.byclass))
        if (flat.measures.logical) {
            hie.perf.mult.maj$Accuracy[1] <- conf.matr.overall$Accuracy
            hie.perf.mult.maj$Kappa[1] <- conf.matr.overall$Kappa
            hie.perf.mult.maj$AccuracyLower[1] <- conf.matr.overall$AccuracyLower
            hie.perf.mult.maj$AccuracyUpper[1] <- conf.matr.overall$AccuracyUpper
            hie.perf.mult.maj$AccuracyNull[1] <- conf.matr.overall$AccuracyNull
            hie.perf.mult.maj$AccuracyPValue[1] <- conf.matr.overall$AccuracyPValue
            hie.perf.mult.maj$McnemarPValue[1] <- conf.matr.overall$McnemarPValue
        }
        if (flat.measures.logical && by.node) {
            melt.byclass <- melt(conf.matr.byclass)
            colnames(melt.byclass) <- c("term.node", "index", 
                "value")
            melt.byclass$term.node.index <- paste(melt.byclass$term.node, 
                melt.byclass$index, sep = ";")
            for (i in 1:nrow(melt.byclass)) {
                col.num <- match(melt.byclass$term.node.index[i], 
                  colnames(hie.perf.mult.maj))
                hie.perf.mult.maj[1, col.num] <- melt.byclass$value[i]
            }
        }
        if (hie.F.measure.logical) {
            results.hie.F.measure <- HieFMeasure(conf.matr = conf.matr, 
                unique.path = unique.path, beta.h.F = beta.h.F, 
                by.node = by.node)
            for (count.measure in 1:nrow(results.hie.F.measure)) {
                col.num.2 <- match(results.hie.F.measure[count.measure, 
                  "measure"], colnames(hie.perf.mult.maj))
                hie.perf.mult.maj[1, col.num.2] <- results.hie.F.measure[count.measure, 
                  "values"]
            }
        }
    }
    if (mult.perm.logical) {
        cat(paste("\n", "Estimating performance measures for the multiplicative permutations rule ", 
            "\n", sep = ""))
        for (count.perm in 1:nrow(hie.perf.mult.perm)) {
            if (div.logical && round(count.perm/div.print, 0) == 
                count.perm/div.print) {
                cat(paste("\n", "  Estimating performance measures for permutation number: ", 
                  count.perm, sep = ""))
            }
            focal.detail <- hie.perf.mult.perm[count.perm, "crisp.rule"]
            pred.nodes <- crisp.case.class[, match(focal.detail, 
                colnames(crisp.case.class))]
            joined.levels <- JoinLevels(vector.1 = pred.nodes, 
                vector.2 = train.data.acc$obs.term.node)
            pred.nodes <- joined.levels$vector.1
            obse.nodes <- joined.levels$vector.2
            rm(joined.levels)
            conf.matr <- confusionMatrix(data = pred.nodes, reference = obse.nodes, 
                dnn = c("Prediction", "Observed"))
            conf.matr.table <- conf.matr$table
            conf.matr.overall <- as.data.frame(t(conf.matr$overall))
            conf.matr.byclass <- conf.matr$byClass
            rownames(conf.matr.byclass) <- sub("Class: ", "", 
                rownames(conf.matr.byclass))
            if (flat.measures.logical) {
                hie.perf.mult.perm$Accuracy[count.perm] <- conf.matr.overall$Accuracy
                hie.perf.mult.perm$Kappa[count.perm] <- conf.matr.overall$Kappa
                hie.perf.mult.perm$AccuracyLower[count.perm] <- conf.matr.overall$AccuracyLower
                hie.perf.mult.perm$AccuracyUpper[count.perm] <- conf.matr.overall$AccuracyUpper
                hie.perf.mult.perm$AccuracyNull[count.perm] <- conf.matr.overall$AccuracyNull
                hie.perf.mult.perm$AccuracyPValue[count.perm] <- conf.matr.overall$AccuracyPValue
                hie.perf.mult.perm$McnemarPValue[count.perm] <- conf.matr.overall$McnemarPValue
            }
            if (flat.measures.logical && by.node) {
                melt.byclass <- melt(conf.matr.byclass)
                colnames(melt.byclass) <- c("term.node", "index", 
                  "value")
                melt.byclass$term.node.index <- paste(melt.byclass$term.node, 
                  melt.byclass$index, sep = ";")
                for (i in 1:nrow(melt.byclass)) {
                  col.num <- match(melt.byclass$term.node.index[i], 
                    colnames(hie.perf.mult.perm))
                  hie.perf.mult.perm[count.perm, col.num] <- melt.byclass$value[i]
                }
            }
            if (hie.F.measure.logical) {
                results.hie.F.measure <- HieFMeasure(conf.matr = conf.matr, 
                  unique.path = unique.path, beta.h.F = beta.h.F, 
                  by.node = by.node)
                for (count.measure in 1:nrow(results.hie.F.measure)) {
                  col.num.2 <- match(results.hie.F.measure[count.measure, 
                    "measure"], colnames(hie.perf.mult.perm))
                  hie.perf.mult.perm[count.perm, col.num.2] <- results.hie.F.measure[count.measure, 
                    "values"]
                }
            }
        }
    }
    if (step.maj.rule.logical) {
        cat(paste("\n", "Estimating performance measures for the stepwise majority rule ", 
            "\n", sep = ""))
        focal.detail <- "stepwise.majority.rule"
        pred.nodes <- crisp.case.class[, match(focal.detail, 
            colnames(crisp.case.class))]
        joined.levels <- JoinLevels(vector.1 = pred.nodes, vector.2 = train.data.acc$obs.term.node)
        pred.nodes <- joined.levels$vector.1
        obse.nodes <- joined.levels$vector.2
        rm(joined.levels)
        conf.matr <- confusionMatrix(data = pred.nodes, reference = obse.nodes, 
            dnn = c("Prediction", "Observed"))
        conf.matr.table <- conf.matr$table
        conf.matr.overall <- as.data.frame(t(conf.matr$overall))
        conf.matr.byclass <- conf.matr$byClass
        rownames(conf.matr.byclass) <- sub("Class: ", "", rownames(conf.matr.byclass))
        if (flat.measures.logical) {
            hie.perf.step.maj$Accuracy[1] <- conf.matr.overall$Accuracy
            hie.perf.step.maj$Kappa[1] <- conf.matr.overall$Kappa
            hie.perf.step.maj$AccuracyLower[1] <- conf.matr.overall$AccuracyLower
            hie.perf.step.maj$AccuracyUpper[1] <- conf.matr.overall$AccuracyUpper
            hie.perf.step.maj$AccuracyNull[1] <- conf.matr.overall$AccuracyNull
            hie.perf.step.maj$AccuracyPValue[1] <- conf.matr.overall$AccuracyPValue
            hie.perf.step.maj$McnemarPValue[1] <- conf.matr.overall$McnemarPValue
        }
        if (flat.measures.logical && by.node) {
            melt.byclass <- melt(conf.matr.byclass)
            colnames(melt.byclass) <- c("term.node", "index", 
                "value")
            melt.byclass$term.node.index <- paste(melt.byclass$term.node, 
                melt.byclass$index, sep = ";")
            for (i in 1:nrow(melt.byclass)) {
                col.num <- match(melt.byclass$term.node.index[i], 
                  colnames(hie.perf.step.maj))
                hie.perf.step.maj[1, col.num] <- melt.byclass$value[i]
            }
        }
        if (hie.F.measure.logical) {
            results.hie.F.measure <- HieFMeasure(conf.matr = conf.matr, 
                unique.path = unique.path, beta.h.F = beta.h.F, 
                by.node = by.node)
            for (count.measure in 1:nrow(results.hie.F.measure)) {
                col.num.2 <- match(results.hie.F.measure[count.measure, 
                  "measure"], colnames(hie.perf.step.maj))
                hie.perf.step.maj[1, col.num.2] <- results.hie.F.measure[count.measure, 
                  "values"]
            }
        }
    }
    cat(paste("\n", "-->  Performance measures estimated  ", 
        "\n", sep = ""))
    cat("######################")
    for (i in 1:nrow(crisp.case.class)) {
        crisp.case.class[i, "obs.term.node"] <- train.data.acc[match(train.data.acc[i, 
            1], crisp.case.class[, 2]), "obs.term.node"]
    }
    hie.performance = rbind(if (step.maj.rule.logical) 
        hie.perf.step.maj, if (mult.maj.rule.logical) 
        hie.perf.mult.maj, if (mult.perm.logical) 
        hie.perf.mult.perm)
    return.list <- list(raw.votes = focal.votes, crisp.case.class = crisp.case.class, 
        hie.performance = hie.performance)
    if (mult.maj.rule.logical || mult.perm.logical) {
        new.list <- list(multiplicative.prop = multiplicative.prop)
        return.list <- c(return.list, new.list)
    }
    if (by.node) {
        new.list <- list(nodes.measures.columns = nodes.measures.columns)
        return.list <- c(return.list, new.list)
    }
    new.list <- list(call = match.call())
    return.list <- c(return.list, new.list)
    return(return.list)
}

