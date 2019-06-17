#' k-fold cross validation for logistic (fixed or mixed-effects) model
#' 
#' \code{runGLM} is to run k-fold cross validation for logistic (fixed or
#' mixed-effects) model.
#'
#' @param tse A TreeSummarizedExperiment object
#' @param assayNum A numeric vector. It specifies which matrix-like elements
#'   in assays will be used to do analysis. If NULL, the first one is used.
#' @param onRow TRUE or FALSE. If TRUE, entities are on rows and samples on
#'   columns.
#' @param kFold A numeric value to specify k-fold cross validation.
#' @param confounder The confounder variables. The names of columns in
#'   \code{colData} if \code{onRow=TRUE} that store the confounders.
#' @param respY The response variables. The names of columns in
#'   \code{colData} if \code{onRow=TRUE} that store the response variable.
#' @param lmeTip TRUE or FALSE. Default is FALSE. If TRUE, logistic
#'   mixed-effects model is used on the leaf nodes.
#' @param lmeInNode TRUE or FALSE. Default is FALSE. If TRUE, logistic
#'   mixed-effects model is used on the internal nodes.
#' @param familyTip "binomial" or "quasibinomial".
#' @param familyInNode "binomial" or "quasibinomial"
#' @param message TRUE or FALSE
#' 
#' @import TreeSummarizedExperiment
#' @importFrom methods is
#' @importFrom stats as.formula glm predict coef
#' @importFrom dplyr filter
#' @importFrom lme4 glmer 
#' 
#' @export
#' @return a data frame
#' @author Ruizhu Huang
#' @examples
#' library(ape)
#' library(ggtree)
#' 
#' n <- 20
#' nSam <- c(50, 50)
#' lib <- 500
#' set.seed(3)
#' exTree <- rtree(n)
#' ggtree(exTree) + geom_text2(aes(label = node))
#' prob <- rep(1/n, n)
#' names(prob) <- exTree$tip.label
#' lse <- simMult(pr = prob, libSize = lib, 
#'                tree = exTree,
#'                minTip.A = 3, maxTip.A = 6,
#'                nSam = nSam, scenario = "S1")
#' 
#' tse2 <- aggData(x = lse, onRow = TRUE)
#' out2 <- cvGLM(tse = tse2, onRow = TRUE, kFold = 10,
#'               lmeTip = FALSE, lmeInNode = FALSE,
#'               familyTip = "quasibinomial", 
#'               familyInNode = "quasibinomial")
#' 
cvGLM <- function(tse, assayNum = NULL,
                  onRow = TRUE, kFold = NULL, 
                  confounder = NULL, respY = "group",
                  lmeTip = FALSE, lmeInNode = TRUE,
                  familyTip = "binomial", 
                  familyInNode = "binomial",
                  message = FALSE) {


    # ================ check inputs ===================================
    # The input data should be a TreeSummarizedExperiment object
    # It should includes only data on the leaf nodes
    if (!is(tse, "TreeSummarizedExperiment")) {
        stop("tse should be a TreeSummarizedExperiment.")
    }
    
    # If not specified, the first table in the assays is used
    if (is.null(assayNum)) {
        assayNum <- 1
    }
    
    if (onRow) {
        rY <- colData(tse)[, respY]
        uY <- unique(rY)
        sp <- ncol(tse)
        
    } else {
        rY <- rowData(tse)[, respY]
        uY <- unique(rY)    
        sp <- nrow(tse)
    }
    
    if (length(uY) != 2) {
        stop(length(uY), "levels detetected in the response variable;
             Two levels are required to do logistic regression.")}
    
    if (kFold > sp) {
        stop("kFold should be less than the number of samples.")
    }
    
    
    if (is.null(kFold)) {
        kFold <- sp
    }
    if (sp < 10 & kFold <10) {
        warning("leave-one-out cross validation is recommended ")
    }
    
    
    # ================create level ===================================
    if (onRow) {
        lev <- unique(rowLinks(tse)$levelNum)
        if (is.null(lev)) {
            lev <- unique(rowLinks(tse)$nodeNum)
        }
    } else {
        lev <- unique(colLinks(tse)$levelNum)
        if (is.null(lev)) {
            lev <- unique(colLinks(tse)$nodeNum)
        }
    }
    
    # ================ run k fold cross validation ===========================
    if (message) {
        message("Preparing train data and test data ...")
    }
    ind <- .createFolds(m = sp, k = kFold)
    errDat <- matrix(NA, nrow = length(lev), ncol = kFold)
    colnames(errDat) <- paste0("err.cv", seq_len(kFold))
    
    if (message) {
        message("Running ", kFold, " cross validation ...")
    }
    for (i in seq_len(kFold)) {
        if (message) {
            message("----------------- ", i)
        }
        
        ts <- unlist(ind[i])
        tr <- unlist(ind[-i])

        if (onRow) {
            test <- convertDf(tse = tse[, ts], onRow = onRow,
                               assayNum = assayNum)
            train <- convertDf(tse = tse[, tr], onRow = onRow,
                               assayNum = assayNum)
        } else {
            test <- convertDf(tse = tse[ts, ], onRow = onRow,
                              assayNum = assayNum)
            train <- convertDf(tse = tse[tr, ], onRow = onRow,
                               assayNum = assayNum)
        }

        test[[respY]] <- factor(test[[respY]])
        train[[respY]] <- factor(train[[respY]])
        
        errDat[, i] <- .singleCV(tse = tse, train, test,
                                 onRow = onRow, 
                                 level = lev,
                                 confounder = confounder, respY = respY,
                                 lmeTip = lmeTip, lmeInNode = lmeInNode,
                                 familyTip = familyTip, 
                                 familyInNode = familyInNode,
                                 message = message) 
    }

    out <- data.frame(errDat)
    out$node <- lev
    return(out)
}


.singleCV <- function(tse, train, test, onRow = TRUE, level,
                     confounder = NULL, respY = "group",
                     lmeTip = FALSE, lmeInNode = TRUE,
                     familyTip = "binomial", 
                     familyInNode = "binomial", 
                     message = message) {
    
    # ================create formula=======================================
    # if (is.null(confounder)) {
    #     if (onRow) {
    #         confounder <- setdiff(colnames(colData(tse)), respY)
    #     } else {
    #         confounder <- setdiff(colnames(rowData(tse)), respY)
    #     }
    # }
    
    formula.glm <- as.formula(paste(respY, " ~",
                                    paste(c("value", confounder),
                                          collapse = " + ")))
    
    formula.glmer <- as.formula(paste(respY, " ~",
                                      paste(c("value", confounder, "(1|node)"),
                                            collapse = " + "),
                                      sep = ""))
    
   
    
    # on leaf nodes
    trainL <- train[train$isLeaf, ]
    levL <- intersect(level, unique(trainL$level))
    lenL <- length(levL)
    trainList <- split(trainL, f = factor(trainL$level, levels = levL))
    
    testL <- test[test$isLeaf, ]
    testList <- split(testL, f = factor(testL$level, levels = levL))
    
    if (message) {
        message("   working on ", lenL, " leaves...",
                "\r", appendLF = FALSE)
        flush.console()
    }
    if (lmeTip) {
        errL <- .cvMix(train = trainList, test = testList, 
                      family = familyTip, fm = formula.glmer,
                      respY = respY, message = message)  
        
    } else {
        errL <- .cvFix(train = trainList, test = testList, 
                      family = familyTip, fm = formula.glm,
                      respY = respY, message = message)
    }

    
    # on internal nodes
    trainI <- train[!train$isLeaf, ]
    levI <- intersect(level, unique(trainI$level))
    lenI <- length(levI)
    trainList <- split(trainI, f = factor(trainI$level, levels = levI))
    
    testI <- test[!test$isLeaf, ]
    testList <- split(testI, f = factor(testI$level, levels = levI))
    
    if (message) {
        message("   working on ", lenI, " internal nodes...",
                "\r", appendLF = FALSE)
        flush.console()
    }
    if (lmeInNode) {
        errI <- .cvMix(train = trainList, test = testList, 
                       family = familyInNode, fm = formula.glmer,
                       respY = respY, message = message)  
        
    } else {
        errI <- .cvFix(train = trainList, test = testList, 
                       family = familyInNode, fm = formula.glm,
                       respY = respY, message = message)
    }
    
    err <- c(errL, errI)[match(c(levL, levI), level)]
    return(err)
    
    # for (i in seq_along(lev)) {
    #     train.i <- train %>%
    #         filter(level == lev[i])
    #     test.i <- test %>%
    #         filter(level == lev[i])
    #     
    #     isTip <- unique(train.i$isLeaf)
    #     
    #     if (isTip) {
    #         if (lmeTip) {
    #             # logistic mixed effects model on inner nodes
    #             mod <- glmer(formula.glmer, family = familyTip, 
    #                          data = train.i)
    #         } else {
    #             # logistic model on leaf
    #             mod <- glm(formula.glm, family = familyTip, 
    #                        data = train.i, control = list(maxit = 50))
    #         }
    #         
    #     } else {
    #         
    #         if (lmeInNode) {
    #             # logistic mixed effects model on inner nodes
    #             mod <- glmer(formula.glmer, family = familyInNode, 
    #                          data = train.i)
    #         } else {
    #             # logistic model on inner nodes
    #             mod <- glm(formula.glm, family = familyInNode, 
    #                        data = train.i, control = list(maxit = 50)) 
    #         }
    #         
    #     }
    #     predG <- predict(mod, newdata = test.i, type = "response")
    #     predGv <- ifelse(predG > 0.5, levels(train.i$group)[2], 
    #                      levels(train.i$group)[1])
    #     err[i] <- sum(predGv != test.i$group)/nrow(test.i)
    # }
    
    
}


.cvFix <- function(train, test, family, fm, message = TRUE, respY) {
    
    lh <- levels(train[[1]][[respY]])[2]
    ll <- levels(train[[1]][[respY]])[1]
    len <- length(train)
    
    out <- lapply(seq_along(train), FUN = function(x) {
        if (message) {
            message("   ", x, " out of ", len,
                    " finished", "\r", appendLF = FALSE)
            flush.console()
        }
        
        mod <- glm(fm, family = family, 
                   data = train[[x]],
                   control = list(maxit = 50))
        predG <- predict(mod, newdata = test[[x]], type = "response")
        
        predGv <- ifelse(predG > 0.5, lh, ll)
        sum(predGv != test[[x]][[respY]])/nrow(test[[x]])
    })
    
    out <- unlist(out)
    return(out)
    
}

.cvMix <- function(train, test, family, fm, message = TRUE, respY) {
    
    lh <- levels(train[[1]][[respY]])[2]
    ll <- levels(train[[1]][[respY]])[1]
    len <- length(train)
    
    out <- lapply(seq_along(train), FUN = function(x) {
        if (message) {
            message("    ", x, " out of ", len,
                    " finished", "\r", appendLF = FALSE)
            flush.console()
        }
        
        mod <- glmer(fm, family = family, 
                     data = train[[x]])
        predG <- predict(mod, newdata = test[[x]], type = "response")
        predGv <- ifelse(predG > 0.5, lh, ll)
        sum(predGv != test[[x]][[respY]])/nrow(test[[x]])
    })
    
    out <- unlist(out)
    return(out)
    
}
