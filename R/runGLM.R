#' run logistic (mixed-effects) model
#' 
#' \code{runGLM} is to run logistic model analysis for hierarchical data. For
#' data on the leaf level, the logistic fixed-effect model is used. However, for
#' data on the internal nodes, uses could decide whether the mixed-effects or
#' the fixed effect model should be used.
#'
#' @param tse A TreeSummarizedExperiment object
#' @param onRow TRUE or FALSE
#' @param confounder The confounder variables.
#' @param assayNum A numeric vector. It specifies which matrix-like elements
#'   in assays will be used to do analysis. If NULL, the first one is used.
#' @param respY The response variables.
#' @param lmeTip TRUE or FALSE. Default is FALSE. If TRUE, mixed-effects model
#'   is used on the leaf nodes.
#' @param lmeInNode TRUE or FALSE. Default is FALSE. If TRUE, mixed-effects model
#'   is used on the internal nodes.
#' @param familyTip "binomial" or "quasibinomial".
#' @param familyInNode "binomial" or "quasibinomial"
#' 
#' @importFrom stats as.formula glm coef
#' @importFrom dplyr filter
#' @importFrom lme4 glmer
#' @importFrom utils flush.console
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
#' tse <- collectData(x = lse, onRow = TRUE)
#' 
#' out <- runGLM(tse = tse, onRow = TRUE)
#' 

# runGLM <- function(tse, onRow = TRUE, confounder = NULL, 
#                    assayNum = NULL, respY = "group",
#                    lmeTip = FALSE, lmeInNode = TRUE,
#                    familyTip = "binomial", 
#                    familyInNode = "binomial",
#                    message = TRUE) {
#     if (message) {
#         message("Preparing data... ")}
#     
#     # A data frame
#     df <- convertDf(tse = tse, onRow = onRow, assayNum = assayNum)
#     if (!respY %in% colnames(df)) {
#         stop(respY, "doesn't exist.")
#     }
#     
#     df[[respY]] <- factor(df[[respY]])
#     
#     # if explanatoryX is null, use all columns.
#     # if (is.null(confounder)) {
#     #     if (onRow) {
#     #         confounder <- setdiff(colnames(colData(tse)), respY)
#     #     } else {
#     #         confounder <- setdiff(colnames(rowData(tse)), respY)
#     #     }
#     # } 
#     if (message) {
#         message("Generating model formula ... ")}
#     formula.glm <- as.formula(paste(respY, " ~", 
#                                     paste(c("value", confounder),
#                                           collapse = " + ")))
#     
#     formula.glmer <- as.formula(paste(respY, " ~", 
#                                       paste(c("value", confounder, "(1|node)"),
#                                             collapse = " + "),
#                                       sep = ""))
#     
#     
#     lev <- unique(df$level)
#     if (is.null(lev)) {
#         lev <- unique(df$node)
#     }
#     
#     out <- data.frame("node" = lev, 
#                       "Estimate" = rep(NA, length(lev)),
#                       "Std. Error" = rep(NA, length(lev)),
#                       "Pr(>|z|)" = rep(NA, length(lev)),
#                       check.names = FALSE)
#     
#     for (i in seq_along(lev)) {
#         if (message) {
#         message(i, " out of ", length(lev),
#                 " finished", "\r", appendLF = FALSE)
#         flush.console()
#         }
#         
#          dfi <- df %>%
#             dplyr::filter(level == lev[i]) 
#         
#         isTip <- unique(dfi$isLeaf)
#         
#         if (isTip) {
#             if (lmeTip) {
#                 # logistic mixed effects model on inner nodes
#                 mod <- glmer(formula.glmer, family = familyTip, 
#                              data = dfi)
#             } else {
#                 # logistic model on leaf
#                 mod <- glm(formula.glm, family = familyTip, 
#                            data = dfi, control = list(maxit = 50))
#             }
#             
#         } else {
#             
#             if (lmeInNode) {
#                 # logistic mixed effects model on inner nodes
#                 mod <- glmer(formula.glmer, family = familyInNode, 
#                              data = dfi)
#             } else {
#                 # logistic model on inner nodes
#                 mod <- glm(formula.glm, family = familyInNode, 
#                            data = dfi, control = list(maxit = 50)) 
#             }
#             
#         }
#         if ("value" %in% rownames(coef(summary(mod)))) {
#             if ("Pr(>|t|)" %in% colnames(coef(summary(mod)))) {
#                 out$`Pr(>|z|)`[i] <- coef(summary(mod))["value", "Pr(>|t|)"]
#             } else {
#                 out$`Pr(>|z|)`[i] <- coef(summary(mod))["value", "Pr(>|z|)"]
#             }
#             
#             out$Estimate[i] <- coef(summary(mod))["value", "Estimate"]
#             out$`Std. Error`[i] <- coef(summary(mod))["value", "Std. Error"]
#         } else {
#             warning("Exact the same value in two groups")
#             out$`Pr(>|z|)`[i] <- NA
#             out$Estimate[i] <- NA
#             out$`Std. Error`[i] <- NA
#         }
#         
#         
#         # prediction
#         # predG <- predict(mod, newdata = dfi, type = "response")
#         # err[[i]] <- sum(round(predG) == dfi$group)/nrow(dfi)
#         # cross validation
#         
#     }
#     
#    return(out)
# }

runGLM <- function(tse, onRow = TRUE, confounder = NULL, 
                    assayNum = NULL, respY = "group",
                    lmeTip = FALSE, lmeInNode = TRUE,
                    familyTip = "binomial", 
                    familyInNode = "binomial",
                    message = TRUE) {
    if (message) {
        message("Preparing data... ")}
    
    # A data frame
    df <- convertDf(tse = tse, onRow = onRow, assayNum = assayNum)
    if (!respY %in% colnames(df)) {
        stop(respY, "doesn't exist.")
    }
    
    df[[respY]] <- factor(df[[respY]])
    # if (familyTip %in% c("binomial", "quasibinomial")) {
    #     if (!is.numeric(df[[respY]])) {
    #         stop("The column", respY, "should be binary numeric value")
    #     }
    # }
    
    # ====================== Formula =============================
    if (message) {
        message("Generating model formula ... ")}
    
    formula_glm <- as.formula(paste(respY, " ~", 
                                    paste(c("value", confounder),
                                          collapse = " + ")))
    
    formula_glmer <- as.formula(paste(respY, " ~", 
                                      paste(c("value", confounder, "(1|node)"),
                                            collapse = " + "),
                                      sep = ""))
    if (message) {
        message("running analysis on leaf nodes... ")}
    
    lev <- unique(df$level)
    if (is.null(lev)) {
        lev <- unique(df$node)
    }
    
    # ================== analysis on leaf nodes ==============
    out <- data.frame("node" = lev, 
                      "Estimate" = rep(NA, length(lev)),
                      "Std. Error" = rep(NA, length(lev)),
                      "Pr(>|z|)" = rep(NA, length(lev)),
                      check.names = FALSE)
    tg <- c("Estimate", "Std. Error", "Pr(>|z|)" )
    
    
    # on leaf nodes
    dL <- df[df$isLeaf, ]
    nL <- unique(dL$level)
    lenL <- length(nL)
    dLL <- split(dL, f = factor(dL$level, levels = nL))
    
    if (lmeTip) {
        for (i in seq_len(lenL)) {
            if (message) {
                message(i, " out of ", lenL,
                        " finished", "\r", appendLF = FALSE)
                flush.console()
            }
            out[i, tg] <- .runMix(df = dLL[[i]], family = familyTip, 
                                  fm = formula_glmer)
        }
        
    } else {
        for (i in seq_len(lenL)) {
            if (message) {
                message(i, " out of ", lenL,
                        " finished", "\r", appendLF = FALSE)
                flush.console()
            }
            out[i, tg] <- .runFix(df = dLL[[i]], family = familyTip, 
                                  fm = formula_glm)
        }
    }
    
    # ================== analysis on internal nodes ==============
    if (message) {
        message("running analysis on internal nodes... ")}
    # on internal nodes
    dI <- df[!df$isLeaf, ]
    nI <- unique(dI$level)
    lenI <- length(nI)
    dII <- split(dI, f = factor(dI$level, levels = nI))
    
    
    if (lmeInNode) {
        for (i in seq_len(lenI)) {
            if (message) {
                message(i, " out of ", lenI,
                        " finished", "\r", appendLF = FALSE)
                flush.console()
            }
            out[i + lenL, tg] <- .runMix(df = dII[[i]], family = familyInNode,
                                         fm = formula_glmer)
        }
        
    } else {
        for (i in seq_len(lenI)) {
            if (message) {
                message(i, " out of ", lenI,
                        " finished", "\r", appendLF = FALSE)
                flush.console()
            }
            out[i + lenL, tg] <- .runFix(df = dII[[i]], family = familyInNode,
                                         fm = formula_glm)
        }
        
    }
    
    return(out)
}



.runMix <- function(df, family, fm) {
    mod <- glmer(fm, family = family, 
                 data = df)
    res <- coef(summary(mod))
    rn <- rownames(res)
    cn <- colnames(res)
    sn <- cn[grepl(pattern = "Pr", cn)]
    
    if ("value" %in% rn) {
        vv <- res["value", c("Estimate", "Std. Error", sn)]
    } else {
        vv <- rep(NA, 3)
    }
    
}

.runFix <- function(df, family, fm) {
    mod <- glm(fm, family = family, 
               data = df, control = list(maxit = 50))
    res <- coef(summary(mod))
    rn <- rownames(res)
    cn <- colnames(res)
    sn <- cn[grepl(pattern = "Pr", cn)]
    
    if ("value" %in% rn) {
        vv <- res["value", c("Estimate", "Std. Error", sn)]
    } else {
        vv <- rep(NA, 3)
    }
}



