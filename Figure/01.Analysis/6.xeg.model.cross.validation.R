### HISTORY #####################################################################
# This script fits the XGBoost model that predicts XEG status from gene, tumour
# and somatic features (Figure 2k, Supplementary Figure 12, Supplementary Note 3).
# Date: 2026-08-04

### DESCRIPTION #################################################################
# Repeated stratified five-fold outer cross-validation with a training-only
# inner validation split for feature selection. The fold-averaged summaries
# written to the output directory are the xeg.model.* objects described in
# 'File information.txt' and read by Figure2k.R.
#
# Usage: Rscript 6.xeg.model.cross.validation.R <input.rda> <output.dir> \
#            <n.repeats> <n.folds> <nrounds> <drop.threshold> <max.iter> \
#            <max.workers> <inner.frac> <nthread> <save.predictions>
# The published run used n.repeats = 10, n.folds = 5 (50 folds).


# Design: repeated stratified five-fold outer cross-validation with a
# training-only inner validation split. Within each outer training fold, a
# patient-level 75/25 holdout is carved out and iterative ablation-based feature
# selection is scored on that inner holdout. The outer test fold is never used
# for feature selection or tuning; it is scored once, after the feature set is
# frozen, to give the outer-test AUROC and the reporting-only SHAP, ablation,
# single-feature and block-ablation values that feed back into nothing.
#
# Fold membership is seeded by 100 * repeat_id, so folds are reproducible.
#
# Memory: each worker holds several multi-million-row subsets. Inner objects are
# released before the outer refit allocates. Measure worker RSS before scaling
# up the worker count.

### PREAMBLE ###################################################################
options(stringsAsFactors = FALSE);

args <- commandArgs(trailingOnly = TRUE);

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x) || x == '') y else x;

get.arg <- function(i, default = NULL) {
    if (length(args) >= i) args[[i]] else default;
    }

input.rda      <- get.arg(1, '2026-03-17_nested_cv_inputs.rda');
output.dir     <- get.arg(2, file.path(getwd(), paste0('innerselect_run_', format(Sys.time(), '%Y%m%d_%H%M%S'))));
n.repeats      <- as.integer(get.arg(3, 1));
n.folds        <- as.integer(get.arg(4, 5));
fixed.nrounds  <- as.integer(get.arg(5, 100));
drop.threshold <- as.numeric(get.arg(6, 0.001));
max.iter       <- as.integer(get.arg(7, 10));
max.workers    <- as.integer(get.arg(8, 5));
inner.frac     <- as.numeric(get.arg(9, 0.25));   # fraction of outer-train patients held out for selection
task.nthread   <- as.integer(get.arg(10, 4));
# 0 = AUROC only (leanest), 1 = keep sample, 2 = keep sample+gene
save.predictions <- as.integer(get.arg(11, 1));
# Shard bounds for multi-node runs: this process handles only tasks
# task_from..task_to (1-based, inclusive) out of n_repeats*n_folds. Defaults
# cover everything, so single-node behaviour is unchanged. Without these, every
# node would pick up every unfinished fold and duplicate the work.
task.from <- as.integer(get.arg(12, 1));
task.to   <- as.integer(get.arg(13, n.repeats * n.folds));
# Skip the aggregation step (set 0 on worker nodes; run once at the end on one
# machine after every shard has finished).
do.aggregate <- as.integer(get.arg(14, 1));

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE);
task.dir <- file.path(output.dir, 'task_results');
dir.create(task.dir, recursive = TRUE, showWarnings = FALSE);
log.file <- file.path(output.dir, 'run.log');

append.log <- function(...) {
    msg <- paste0(...);
    cat(sprintf("[%s] %s\n", format(Sys.time(), '%Y-%m-%d %H:%M:%S'), msg),
        file = log.file, append = TRUE);
    }

log.task.step <- function(repeat.id, fold.id, iter.idx, msg) {
    append.log(sprintf('STEP repeat=%s fold=%s iter=%s :: %s',
                       repeat.id, fold.id, iter.idx, msg));
    }

# Load necessary libraries
suppressPackageStartupMessages({
    library(Matrix);
    library(xgboost);
    library(pROC);
    library(foreach);
    library(doParallel);
    });

### DATA PREPARATION ###########################################################
load(input.rda);

required.objects <- c(
    'rf.combined.df.v2.new',
    'base.vars.v2.new',
    'params_iter_base',
    'run_refit_iter',
    'esc_re_it',
    'map_onehot_to_base_it'
    );
missing.objects <- setdiff(required.objects, ls());
if (length(missing.objects) > 0) {
    stop('Missing required objects in input_rda: ', paste(missing.objects, collapse = ', '));
    }

max.workers <- max(1L, max.workers);
task.nthread <- max(1L, task.nthread);

msg.lines <- c(
    paste0('Input RDA: ', normalizePath(input.rda, mustWork = FALSE)),
    paste0('Output dir: ', normalizePath(output.dir, mustWork = FALSE)),
    paste0('Repeats: ', n.repeats),
    paste0('Folds: ', n.folds),
    paste0('Fixed nrounds: ', fixed.nrounds),
    paste0('Drop threshold: ', drop.threshold),
    paste0('Max iterations: ', max.iter),
    paste0('Workers: ', max.workers),
    paste0('Inner holdout fraction: ', inner.frac),
    paste0('Threads per task: ', task.nthread),
    'Selection scored on: INNER HOLDOUT (outer test used once, for final refit only)'
    );
for (ln in msg.lines) { message(ln);; append.log(ln); }

### FEATURE BLOCKS #############################################################
# Feature blocks for Supplementary Figure 12c, validated before the fork.
#
# A typo here would not error: intersect() would return empty, the block would
# record AUC_Drop = 0 in every fold, and the figure would show a zero bar that
# reads as "this category contributes nothing". Validate against the real
# feature list up front so that failure mode is impossible.
BLOCK.TEMPLATE <- list(
    'Gene Intrinsic' = c('gene_length', 'exon_num', 'GC_content', 'chromosome'),
    'Gene Extrinsic' = c('CNA_raw', 'Meth_z', 'ecDNA', 'SV'),
    'Tumor'          = c('subtype', 'stage', 'age')
    );
local({
    somatic <- grep('^(driver_|xeg_snv)', base.vars.v2.new, value = TRUE);
    unknown <- setdiff(unlist(BLOCK.TEMPLATE, use.names = FALSE), base.vars.v2.new);
    if (length(unknown) > 0) {
        stop('Block features absent from base.vars.v2.new: ',
             paste(unknown, collapse = ', '));
        }
    if (length(somatic) == 0) stop('Somatic Landscape regex matches no feature');
    unassigned <- setdiff(base.vars.v2.new,
                          c(unlist(BLOCK.TEMPLATE, use.names = FALSE), somatic));
    if (length(unassigned) > 0) {
        append.log('WARNING features in no block: ',
                   paste(unassigned, collapse = ', '));
        message('WARNING features in no block: ',
                paste(unassigned, collapse = ', '));
        } else {
        append.log(sprintf('Block validation OK: all %d features assigned',
                           length(base.vars.v2.new)));
        }
    });

### CROSS-VALIDATION FOLDS #####################################################
# Patient-level outer folds, stratified by cohort and XEG status.
make.stratified.folds <- function(df, n.folds, seed) {
    patient.info <- unique(df[, c('sample', 'cohort'), drop = FALSE]);
    xeg.by.patient <- aggregate(
        XEG ~ sample,
        data = df,
        FUN = function(x) as.integer(any(as.character(x) %in% c('1', 'yes')))
        );
    colnames(xeg.by.patient)[2] <- 'has_XEG';
    patient.info <- merge(patient.info, xeg.by.patient, by = 'sample');
    patient.info$strat_key <- paste(patient.info$cohort, patient.info$has_XEG, sep = '_');
    patient.info$fold <- NA_integer_;

    set.seed(seed);
    for (key in unique(patient.info$strat_key)) {
        idx <- which(patient.info$strat_key == key);
        shuffled <- sample(idx);
        patient.info$fold[shuffled] <- rep(seq_len(n.folds), length.out = length(idx));
        }

    patient.info;
    }

# Patient-level, cohort+XEG-stratified inner holdout carved out of the outer train set.
make.inner.split <- function(df, train.samples, frac, seed) {
    info <- unique(df[df$sample %in% train.samples, c('sample', 'cohort'), drop = FALSE]);
    xeg.by.patient <- aggregate(
        XEG ~ sample,
        data = df[df$sample %in% train.samples, , drop = FALSE],
        FUN = function(x) as.integer(any(as.character(x) %in% c('1', 'yes')))
        );
    colnames(xeg.by.patient)[2] <- 'has_XEG';
    info <- merge(info, xeg.by.patient, by = 'sample');
    info$strat_key <- paste(info$cohort, info$has_XEG, sep = '_');

    set.seed(seed);
    holdout <- character(0);
    for (key in unique(info$strat_key)) {
        idx <- which(info$strat_key == key);
        n.hold <- max(1L, round(length(idx) * frac));
        if (n.hold >= length(idx)) n.hold <- length(idx) - 1L;
        if (n.hold < 1L) next;
        holdout <- c(holdout, info$sample[sample(idx, n.hold)]);
        }
    list(
        inner_train = setdiff(info$sample, holdout),
        inner_hold  = holdout
        );
    }

# Meth_z: mean and sd estimated on the training rows only, applied to the eval set.
add.meth.z.from.training <- function(train.df, test.df) {
    train.df$Meth_z <- NA_real_;
    test.df$Meth_z  <- NA_real_;

    for (coh in unique(train.df$cohort)) {
        idx.coh <- train.df$cohort == coh;
        genes.in.cohort <- unique(train.df$gene[idx.coh & !is.na(train.df$Meth_raw)]);
        for (g in genes.in.cohort) {
            idx.tr <- idx.coh & train.df$gene == g & !is.na(train.df$Meth_raw);
            vals <- train.df$Meth_raw[idx.tr];
            if (length(vals) <= 1) next;
            m <- mean(vals);
            s <- sd(vals);
            if (is.na(s) || s <= 0) next;

            train.df$Meth_z[idx.tr] <- (vals - m) / s;

            idx.te <- test.df$cohort == coh & test.df$gene == g & !is.na(test.df$Meth_raw);
            if (any(idx.te)) {
                test.df$Meth_z[idx.te] <- (test.df$Meth_raw[idx.te] - m) / s;
                }
            }
        }

    list(train = train.df, test = test.df);
    }

strip.refit.result <- function(res) {
    list(auc = res$auc, best_nrounds = res$best_nrounds,
         vars = res$vars, feat_eval = res$feat_eval);
    }

### DESIGN MATRIX ##############################################################
# Shared design-matrix builder.
#
# Categorical levels in the evaluation set that were not observed in training
# are mapped to the reference level rather than NA, so that
# Matrix::sparse.model.matrix() (na.action = na.omit) cannot drop the row and
# every indicator column for that feature is 0 ("none of the known
# categories"). See the inline comment for why not an "unseen" level or the mode.
build.xy <- function(train.df, test.df, variables) {
    keep.cols <- c('XEG', intersect(variables, names(train.df)));
    train.sel <- train.df[, keep.cols, drop = FALSE];
    test.sel  <- test.df[, keep.cols, drop = FALSE];

    for (cc in setdiff(names(train.sel), 'XEG')) {
        if (is.factor(train.sel[[cc]]) || is.character(train.sel[[cc]])) {
            tr <- as.character(train.sel[[cc]]);
            te <- as.character(test.sel[[cc]]);
            tr[is.na(tr)] <- 'missing';
            te[is.na(te)] <- 'missing';
            # Only levels actually observed in training. Adding "missing"/
            # "unseen" unconditionally would emit all-zero indicator columns,
            # and colsample_bytree counts columns.
            lv <- unique(tr);
            unseen <- !(te %in% lv);
            if (any(unseen)) {
                # Map to the REFERENCE level (lv[1]), which under treatment
                # contrasts sets every indicator column for this feature to 0 --
                # i.e. "none of the known categories". Alternatives not used:
                #   - adding an "unseen" level: it appears in NO training row, so
                #     an all-zero training column is emitted, and since
                #     colsample_bytree counts columns that perturbs the fit.
                #   - mapping to the most frequent level: that asserts the row IS
                #     that category.
                te[unseen] <- lv[1];
                append.log(sprintf("NOTE unseen level in %s: %d rows -> '%s'",
                                   cc, sum(unseen), lv[1]));
                }
            train.sel[[cc]] <- factor(tr, levels = lv);
            test.sel[[cc]]  <- factor(te, levels = lv);
            } else {
            train.sel[[cc]] <- suppressWarnings(as.numeric(train.sel[[cc]]));
            test.sel[[cc]]  <- suppressWarnings(as.numeric(test.sel[[cc]]));
            train.sel[[cc]][is.na(train.sel[[cc]])] <- -999;
            test.sel[[cc]][is.na(test.sel[[cc]])]   <- -999;
            }
        }

    y.tr <- as.numeric(as.character(train.sel$XEG) %in% c('1', 'yes'));
    y.te <- as.numeric(as.character(test.sel$XEG) %in% c('1', 'yes'));
    x.tr <- Matrix::sparse.model.matrix(~ . - XEG, data = train.sel);
    x.te <- Matrix::sparse.model.matrix(~ . - XEG, data = test.sel);
    miss <- setdiff(colnames(x.tr), colnames(x.te));
    if (length(miss) > 0) {
        z <- Matrix::Matrix(0, nrow = nrow(x.te), ncol = length(miss), sparse = TRUE);
        colnames(z) <- miss;; x.te <- cbind(x.te, z);
        }
    x.te <- x.te[, colnames(x.tr), drop = FALSE];

    # No row may be silently dropped.
    stopifnot(nrow(x.tr) == length(y.tr), nrow(x.te) == length(y.te),
              nrow(x.tr) == nrow(train.sel), nrow(x.te) == nrow(test.sel));

    list(x.tr = x.tr, x.te = x.te, y.tr = y.tr, y.te = y.te);
    }

### SELECTION REFIT ############################################################
# Refit used during feature selection only.
#
# Selection consumes Ablation_Drop, computed as auc - auc_without_feature on the
# supplied eval set (same seed, same params). SHAP is not computed here: selection
# never reads it, and on a 2.1M-row holdout that dense matrix costs ~390 MB per
# iteration. Gain is still returned for reporting.
run.selection.iter <- function(sel.vars, train.df, test.df,
                               params.base, nrounds.base) {
    m <- build.xy(train.df, test.df, sel.vars);
    x.tr <- m$x.tr;; x.te <- m$x.te;; y.tr <- m$y.tr;; y.te <- m$y.te;
    rm(m);

    spw <- sum(y.tr == 0) / max(1, sum(y.tr == 1));
    params.r <- params.base;
    params.r$scale_pos_weight <- spw;

    d.tr <- xgboost::xgb.DMatrix(x.tr, label = y.tr, missing = -999);
    d.te <- xgboost::xgb.DMatrix(x.te, label = y.te, missing = -999);

    set.seed(123);
    mdl <- xgboost::xgb.train(params = params.r, data = d.tr,
                              nrounds = nrounds.base, verbose = 0);
    auc <- as.numeric(pROC::auc(pROC::roc(y.te, predict(mdl, d.te), quiet = TRUE)));

    imp <- xgboost::xgb.importance(model = mdl);
    gain.v <- setNames(rep(0, length(sel.vars)), sel.vars);
    gmap <- map_onehot_to_base_it(as.character(imp$Feature), sel.vars);
    for (i in seq_len(nrow(imp))) {
        b <- gmap[i];; if (!is.na(b)) gain.v[b] <- gain.v[b] + imp$Gain[i];
        }

    # Ablation: retrain without each feature, rescore on the SAME eval set.
    abl.v <- setNames(rep(NA_real_, length(sel.vars)), sel.vars);
    for (vi in seq_along(sel.vars)) {
        v <- sel.vars[vi];
        kp <- setdiff(sel.vars, v);
        ma <- build.xy(train.df, test.df, kp);
        d.tr.a <- xgboost::xgb.DMatrix(ma$x.tr, label = ma$y.tr, missing = -999);
        d.te.a <- xgboost::xgb.DMatrix(ma$x.te, label = ma$y.te, missing = -999);
        set.seed(123);
        m.a <- xgboost::xgb.train(params = params.r, data = d.tr.a,
                                  nrounds = nrounds.base, verbose = 0);
        auc.a <- as.numeric(pROC::auc(pROC::roc(ma$y.te, predict(m.a, d.te.a),
                                                quiet = TRUE)));
        abl.v[v] <- auc - auc.a;
        rm(ma, d.tr.a, d.te.a, m.a);
        if (vi %% 4 == 0) gc(FALSE);
        }

    feat.eval <- data.frame(
        Variable = sel.vars,
        Gain_sum = as.numeric(gain.v[sel.vars]),
        SHAP_sum = NA_real_,          # not computed during selection
        Ablation_Drop = as.numeric(abl.v[sel.vars]),
        stringsAsFactors = FALSE
        );
    feat.eval <- feat.eval[order(-feat.eval$Ablation_Drop), ];

    rm(x.tr, x.te, d.tr, d.te, mdl);; gc(FALSE);
    list(auc = auc, best_nrounds = nrounds.base,
         feat_eval = feat.eval, vars = sel.vars);
    }

### FINAL FIT ##################################################################
# Final-stage fit: train on the outer training fold with the already-fixed
# feature set, predict once on outer test, score once. Preprocessing is identical
# to the selection stage (same factor handling, same -999 missing sentinel, same
# scale_pos_weight from the training rows only).
fit.selected.model <- function(train.df, test.df, variables, params.base, nrounds,
                               want.shap = TRUE, want.ablation = TRUE,
                               want.uni = TRUE, want.block = TRUE) {
    m <- build.xy(train.df, test.df, variables);
    x.tr <- m$x.tr;; x.te <- m$x.te;; y.tr <- m$y.tr;; y.te <- m$y.te;
    rm(m);

    # class weight from the outer TRAINING rows only
    spw <- sum(y.tr == 0) / max(1, sum(y.tr == 1));
    params.r <- params.base;
    params.r$scale_pos_weight <- spw;

    d.tr <- xgboost::xgb.DMatrix(x.tr, label = y.tr, missing = -999);
    d.te <- xgboost::xgb.DMatrix(x.te, label = y.te, missing = -999);

    set.seed(123);
    mdl <- xgboost::xgb.train(params = params.r, data = d.tr,
                              nrounds = nrounds, verbose = 0);

    pred <- as.numeric(predict(mdl, d.te));
    auc  <- as.numeric(pROC::auc(pROC::roc(y.te, pred, quiet = TRUE)));

    imp <- xgboost::xgb.importance(model = mdl);
    gain.v <- setNames(rep(0, length(variables)), variables);
    gmap <- map_onehot_to_base_it(as.character(imp$Feature), variables);
    for (i in seq_len(nrow(imp))) {
        b <- gmap[i];; if (!is.na(b)) gain.v[b] <- gain.v[b] + imp$Gain[i];
        }

    # SHAP on the outer test fold, for Figure 2k's mean|SHAP| panel. Computed
    # once per fold (not per selection iteration) and reduced to column means
    # immediately so the dense matrix is not retained.
    shap.v <- setNames(rep(NA_real_, length(variables)), variables);
    if (isTRUE(want.shap)) {
        sm <- predict(mdl, d.te, predcontrib = TRUE);
        sn <- colnames(sm)[-ncol(sm)];
        sa <- colMeans(abs(sm[, -ncol(sm), drop = FALSE]));
        rm(sm);; gc(FALSE);
        shap.v[] <- 0;
        smap <- map_onehot_to_base_it(sn, variables);
        for (i in seq_along(sa)) {
            b <- smap[i];; if (!is.na(b)) shap.v[b] <- shap.v[b] + sa[i];
            }
        }

    # Ablation drop on the outer test fold, for Figure 2k's AUROC-Drop panel.
    # This is REPORTING ONLY: selected_vars is already frozen by the caller and
    # nothing downstream re-reads this to change the feature set.
    abl.v <- setNames(rep(NA_real_, length(variables)), variables);
    if (isTRUE(want.ablation) && length(variables) > 1) {
        for (vi in seq_along(variables)) {
            kp <- setdiff(variables, variables[vi]);
            ma <- build.xy(train.df, test.df, kp);
            dta <- xgboost::xgb.DMatrix(ma$x.tr, label = ma$y.tr, missing = -999);
            dea <- xgboost::xgb.DMatrix(ma$x.te, label = ma$y.te, missing = -999);
            set.seed(123);
            mza <- xgboost::xgb.train(params = params.r, data = dta,
                                      nrounds = nrounds, verbose = 0);
            abl.v[variables[vi]] <- auc - as.numeric(
                pROC::auc(pROC::roc(ma$y.te, predict(mza, dea), quiet = TRUE)));
            rm(ma, dta, dea, mza);
            if (vi %% 4 == 0) gc(FALSE);
            }
        }

    rm(x.tr, x.te, d.tr, d.te, mdl);; gc(FALSE);

    # Single-feature AUROC (Figure 2k, left panel): each feature alone, trained
    # on outer-train and scored on outer test.
    uni <- data.frame(Variable = character(0), AUC = numeric(0),
                      stringsAsFactors = FALSE);
    if (isTRUE(want.uni)) {
        for (vi in seq_along(variables)) {
            mu <- build.xy(train.df, test.df, variables[vi]);
            dtu <- xgboost::xgb.DMatrix(mu$x.tr, label = mu$y.tr, missing = -999);
            deu <- xgboost::xgb.DMatrix(mu$x.te, label = mu$y.te, missing = -999);
            set.seed(123);
            mzu <- xgboost::xgb.train(params = params.r, data = dtu,
                                      nrounds = nrounds, verbose = 0);
            uni <- rbind(uni, data.frame(
                Variable = variables[vi],
                AUC = as.numeric(pROC::auc(pROC::roc(mu$y.te, predict(mzu, deu),
                                                     quiet = TRUE))),
                stringsAsFactors = FALSE));
            rm(mu, dtu, deu, mzu);
            if (vi %% 4 == 0) gc(FALSE);
            }
        }

    # Block ablation (Supplementary Figure 12c): drop a whole feature category at
    # once.
    blocks <- c(BLOCK.TEMPLATE,
                list('Somatic Landscape' =
                     grep('^(driver_|xeg_snv)', variables, value = TRUE)));
    blk <- data.frame(Block = character(0), AUC = numeric(0), AUC_Drop = numeric(0),
                      Removed_Features = character(0), Removed_Count = integer(0),
                      stringsAsFactors = FALSE);
    if (isTRUE(want.block)) {
        for (bn in names(blocks)) {
            bv <- intersect(blocks[[bn]], variables);
            if (length(bv) == 0) {
                blk <- rbind(blk, data.frame(Block = bn, AUC = auc, AUC_Drop = 0,
                    Removed_Features = '', Removed_Count = 0L,
                    stringsAsFactors = FALSE));
                next;
                }
            kp <- setdiff(variables, bv);
            if (length(kp) == 0) {
                blk <- rbind(blk, data.frame(Block = bn, AUC = NA_real_,
                    AUC_Drop = NA_real_,
                    Removed_Features = paste(bv, collapse = ', '),
                    Removed_Count = length(bv), stringsAsFactors = FALSE));
                next;
                }
            mb <- build.xy(train.df, test.df, kp);
            dtb <- xgboost::xgb.DMatrix(mb$x.tr, label = mb$y.tr, missing = -999);
            deb <- xgboost::xgb.DMatrix(mb$x.te, label = mb$y.te, missing = -999);
            set.seed(123);
            mzb <- xgboost::xgb.train(params = params.r, data = dtb,
                                      nrounds = nrounds, verbose = 0);
            ab <- as.numeric(pROC::auc(pROC::roc(mb$y.te, predict(mzb, deb),
                                                 quiet = TRUE)));
            blk <- rbind(blk, data.frame(Block = bn, AUC = ab, AUC_Drop = auc - ab,
                Removed_Features = paste(bv, collapse = ', '),
                Removed_Count = length(bv), stringsAsFactors = FALSE));
            rm(mb, dtb, deb, mzb);; gc(FALSE);
            }
        }

    list(
        auc = auc,
        pred = pred,
        y_te = y.te,
        vars = variables,
        scale_pos_weight = spw,
        uni_auc = uni,
        block_ablation = blk,
        feat_eval = data.frame(
            Variable      = variables,
            Gain_sum      = as.numeric(gain.v[variables]),
            SHAP_sum      = as.numeric(shap.v[variables]),
            Ablation_Drop = as.numeric(abl.v[variables]),
            stringsAsFactors = FALSE
            )
        );
    }

### ONE OUTER FOLD #############################################################
run.innerselect.fold <- function(df, train.samples, test.samples, repeat.id, fold.id) {
    inner <- make.inner.split(df, train.samples, inner.frac,
                              seed = 100L * repeat.id + fold.id);

    inner.train.df <- df[df$sample %in% inner$inner_train, , drop = FALSE];
    inner.hold.df  <- df[df$sample %in% inner$inner_hold,  , drop = FALSE];

    log.task.step(
        repeat.id, fold.id, 'setup',
        sprintf('outer_train=%d outer_test=%d inner_train=%d inner_hold=%d',
                length(train.samples), length(test.samples),
                length(inner$inner_train), length(inner$inner_hold))
        );

    # Meth_z for the selection stage: computed from inner-train only.
    mz.in <- add.meth.z.from.training(inner.train.df, inner.hold.df);
    inner.train.df <- mz.in$train;
    inner.hold.df  <- mz.in$test;
    rm(mz.in);; gc(FALSE);   # list still referenced both frames; free it now

    params.base <- params_iter_base;
    params.base$nthread <- task.nthread;

    cur.vars <- base.vars.v2.new;
    iter.results <- list();
    removed.by.iter <- list();
    final.iter.key <- NA_character_;

    # Selection loop: scored on the inner holdout, never on the outer test fold
    for (iter.idx in 0:max.iter) {
        log.task.step(repeat.id, fold.id, iter.idx,
                      sprintf('start n_features=%d', length(cur.vars)));

        res <- run.selection.iter(
            sel.vars     = cur.vars,
            train.df     = inner.train.df,
            test.df      = inner.hold.df,
            params.base  = params.base,
            nrounds.base = fixed.nrounds
            );

        key <- paste0('iter', iter.idx);
        iter.results[[key]] <- strip.refit.result(res);

        neg.vars <- res$feat_eval$Variable[res$feat_eval$Ablation_Drop < drop.threshold];
        removed.by.iter[[key]] <- neg.vars;

        log.task.step(
            repeat.id, fold.id, iter.idx,
            sprintf('done inner_auc=%.4f kept=%d removed=%d',
                    res$auc, length(res$vars), length(neg.vars))
            );

        if (length(neg.vars) == 0 || iter.idx == max.iter) {
            final.iter.key <- key;
            log.task.step(repeat.id, fold.id, iter.idx,
                          sprintf('stop final_iter=%s selected_features=%d',
                                  key, length(res$vars)));
            break;
            }

        cur.vars <- setdiff(cur.vars, neg.vars);
        if (length(cur.vars) == 0) {
            stop('All features were removed at repeat ', repeat.id, ', fold ', fold.id);
            }
        }

    selected.vars <- cur.vars;
    inner.auc <- iter.results[[final.iter.key]]$auc;

    # Release the inner-stage frames BEFORE allocating the outer ones. Without
    # this the worker holds inner_train + inner_hold + train + test
    # simultaneously (~18.9M rows); dropping the inner pair first cuts the peak
    # to ~10.5M rows. `res` holds the last selection model, also dead here.
    rm(inner.train.df, inner.hold.df, res);
    gc(FALSE);

    # Final refit: full outer train -> outer test, touched exactly once
    train.df <- df[df$sample %in% train.samples, , drop = FALSE];
    test.df  <- df[df$sample %in% test.samples,  , drop = FALSE];
    mz.out <- add.meth.z.from.training(train.df, test.df);
    train.df <- mz.out$train;
    test.df  <- mz.out$test;
    rm(mz.out);; gc(FALSE);

    res.final <- fit.selected.model(
        train.df  = train.df,
        test.df   = test.df,
        variables = selected.vars,
        params.base = params.base,
        nrounds   = fixed.nrounds,
        want.shap = TRUE,
        want.ablation = TRUE,
        want.uni = TRUE,
        want.block = TRUE
        );

    pred <- res.final$pred;
    y.te <- res.final$y_te;

    log.task.step(repeat.id, fold.id, 'final',
                  sprintf('outer_test_auc=%.4f n_features=%d',
                          res.final$auc, length(selected.vars)));

    list(
        repeat_id = repeat.id,
        fold_id = fold.id,
        n_train_rows = nrow(train.df),
        n_test_rows = nrow(test.df),
        n_train_samples = length(unique(train.df$sample)),
        n_test_samples = length(unique(test.df$sample)),
        n_inner_train_samples = length(inner$inner_train),
        n_inner_hold_samples = length(inner$inner_hold),
        final_iter = final.iter.key,
        inner_holdout_auc = inner.auc,
        final_auc = res.final$auc,
        final_n_features = length(selected.vars),
        final_vars = selected.vars,
        final_feat_eval = res.final$feat_eval,
        final_uni_auc = res.final$uni_auc,
        final_block_ablation = res.final$block_ablation,
        final_scale_pos_weight = res.final$scale_pos_weight,
        inner_feat_eval = iter.results[[final.iter.key]]$feat_eval,
        removed_by_iter = removed.by.iter,
        iter_results = iter.results,
        # Per-row predictions are ~2.1M rows/fold. save_predictions=0 keeps only
        # what AUROC needs; 1 keeps sample (for patient-level follow-up);
        # 2 keeps sample+gene.
        predictions = if (save.predictions >= 2) {
            data.frame(repeat_id = repeat.id, fold_id = fold.id,
                       sample = test.df$sample, gene = test.df$gene,
                       y_true = y.te, y_pred = pred, stringsAsFactors = FALSE);
            } else if (save.predictions == 1) {
            data.frame(repeat_id = repeat.id, fold_id = fold.id,
                       sample = test.df$sample,
                       y_true = y.te, y_pred = pred, stringsAsFactors = FALSE);
            } else {
            NULL;
            }
        );
    }

### RUN ########################################################################
tasks <- do.call(rbind, lapply(seq_len(n.repeats), function(rep.id) {
    fold.info <- make.stratified.folds(rf.combined.df.v2.new, n.folds = n.folds, seed = 100L * rep.id);
    do.call(rbind, lapply(seq_len(n.folds), function(fold.id) {
        data.frame(
            repeat_id = rep.id,
            fold_id = fold.id,
            train_samples = I(list(fold.info$sample[fold.info$fold != fold.id])),
            test_samples  = I(list(fold.info$sample[fold.info$fold == fold.id])),
            stringsAsFactors = FALSE
            );
        }));
    }));

get.task.output.path <- function(repeat.id, fold.id) {
    file.path(task.dir, sprintf('repeat_%02d_fold_%02d.rds', repeat.id, fold.id));
    }

append.log('Total tasks: ', nrow(tasks));
message('Total tasks: ', nrow(tasks));

pending <- which(!file.exists(vapply(
    seq_len(nrow(tasks)),
    function(i) get.task.output.path(tasks$repeat_id[[i]], tasks$fold_id[[i]]),
    character(1)
    )));

# Restrict to this node's shard. Tasks are ordered repeat-major
# (r1f1, r1f2 ... r1f5, r2f1 ...), so a contiguous index range is a clean split.
task.from <- max(1L, task.from);
task.to   <- min(nrow(tasks), task.to);
if (task.from > task.to) stop('Empty shard: task_from > task_to');
pending <- pending[pending >= task.from & pending <= task.to];
append.log(sprintf('Shard: tasks %d..%d of %d', task.from, task.to, nrow(tasks)));
message(sprintf('Shard: tasks %d..%d of %d', task.from, task.to, nrow(tasks)));

append.log('Pending tasks: ', length(pending));
message('Pending tasks: ', length(pending));

if (length(pending) > 0) {
    workers <- min(max.workers, length(pending));
    cl <- parallel::makeCluster(workers, type = 'FORK');
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE);
    doParallel::registerDoParallel(cl);
    append.log(sprintf('Starting %d task(s) with %d worker(s)', length(pending), workers));
    message(sprintf('Starting %d task(s) with %d worker(s)', length(pending), workers));

    invisible(foreach(i = pending, .packages = c('Matrix', 'xgboost', 'pROC')) %dopar% {
        rep.id  <- tasks$repeat_id[[i]];
        fold.id <- tasks$fold_id[[i]];
        t0 <- Sys.time();
        append.log(sprintf('START repeat=%s fold=%s', rep.id, fold.id));
        out <- tryCatch({
            res <- run.innerselect.fold(
                df = rf.combined.df.v2.new,
                train.samples = tasks$train_samples[[i]],
                test.samples  = tasks$test_samples[[i]],
                repeat.id = rep.id,
                fold.id = fold.id
                );
            res$elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = 'secs'));
            saveRDS(res, file = get.task.output.path(rep.id, fold.id));
            append.log(sprintf(
                'DONE repeat=%s fold=%s :: outer_test_auc=%.4f inner_auc=%.4f n_features=%d elapsed_sec=%.1f',
                rep.id, fold.id, res$final_auc, res$inner_holdout_auc,
                res$final_n_features, res$elapsed_sec));
            TRUE;
            }, error = function(e) {
            append.log(sprintf('ERROR repeat=%s fold=%s :: %s',
                               rep.id, fold.id, conditionMessage(e)));
            FALSE;
            });
        out;
        });

    parallel::stopCluster(cl);
    }

### COMBINE ####################################################################
# Worker nodes (do_aggregate = 0) stop here: their fold rds files are on shared
# storage and one designated run aggregates everything afterwards. Without this
# a shard would hit the completeness guard below and exit non-zero even though
# it did its own work correctly.
if (!identical(do.aggregate, 1L)) {
    append.log(sprintf('Shard %d..%d finished; aggregation skipped.',
                       task.from, task.to));
    message(sprintf('Shard %d..%d finished; aggregation skipped.',
                    task.from, task.to));
    quit(save = 'no', status = 0);
    }

result.files <- vapply(
    seq_len(nrow(tasks)),
    function(i) get.task.output.path(tasks$repeat_id[[i]], tasks$fold_id[[i]]),
    character(1)
    );
result.files <- result.files[file.exists(result.files)];
if (length(result.files) == 0) stop('No task results produced.');

all.res <- lapply(result.files, readRDS);

summary.df <- do.call(rbind, lapply(all.res, function(r) data.frame(
    repeat_id = r$repeat_id,
    fold_id = r$fold_id,
    final_iter = r$final_iter,
    inner_holdout_auc = r$inner_holdout_auc,
    final_auc = r$final_auc,
    final_n_features = r$final_n_features,
    n_train_samples = r$n_train_samples,
    n_test_samples = r$n_test_samples,
    n_inner_train_samples = r$n_inner_train_samples,
    n_inner_hold_samples = r$n_inner_hold_samples,
    elapsed_sec = r$elapsed_sec %||% NA_real_,
    stringsAsFactors = FALSE
    )));
summary.df <- summary.df[order(summary.df$repeat_id, summary.df$fold_id), ];

# A fold that died (OOM, error) writes no rds. Its absence must not pass
# silently as a mean over fewer folds than claimed. Write what we have FIRST so
# a single dead fold does not destroy hours of completed work, then stop.
# Completed rds files persist; re-running the same command resumes only the
# missing folds.
if (length(result.files) < nrow(tasks)) {
    write.csv(summary.df,
              file.path(output.dir, 'innerselect_summary_PARTIAL.csv'),
              row.names = FALSE);
    stop(sprintf(
        paste0('Incomplete run: %d of %d fold files present. Partial summary ',
               'written to innerselect_summary_PARTIAL.csv. Check run.log for ',
               'ERROR lines, then re-run the same command to resume.'),
        length(result.files), nrow(tasks)));
    }

write.csv(summary.df, file.path(output.dir, 'innerselect_summary.csv'), row.names = FALSE);

feat.counts <- table(unlist(lapply(all.res, function(r) r$final_vars)));
freq.df <- data.frame(
    Variable = names(feat.counts),
    Selected_Count = as.integer(feat.counts),
    Selected_Fraction = as.numeric(feat.counts) / length(all.res),
    stringsAsFactors = FALSE
    );
freq.df <- freq.df[order(-freq.df$Selected_Fraction), ];
write.csv(freq.df, file.path(output.dir, 'innerselect_feature_frequency.csv'), row.names = FALSE);

saveRDS(all.res, file.path(output.dir, 'innerselect_combined_results.rds'));

### FIGURE INPUTS ##############################################################
# Fold-averaged summaries, with the column names the plotting scripts read.
# Each quantity is averaged over the folds in which that feature was selected.
avg.by <- function(extract, value.cols, key = 'Variable') {
    d <- do.call(rbind, lapply(all.res, extract));
    if (is.null(d) || nrow(d) == 0) return(NULL);
    agg <- aggregate(d[, value.cols, drop = FALSE],
                     by = setNames(list(d[[key]]), key),
                     FUN = function(x) mean(x, na.rm = TRUE));
    n <- as.data.frame(table(d[[key]]), stringsAsFactors = FALSE);
    names(n) <- c(key, 'N_folds');
    merge(agg, n, by = key, all.x = TRUE);
    }

# Fig 2k middle/right: mean|SHAP|, Gain, per-feature AUROC drop
fe <- avg.by(function(r) r$final_feat_eval,
             c('Gain_sum', 'SHAP_sum', 'Ablation_Drop'));
if (!is.null(fe)) {
    fe <- fe[order(-fe$Ablation_Drop), ];
    write.csv(fe, file.path(output.dir, 'innerselect_feat_eval_foldmean.csv'),
              row.names = FALSE);
    }

# Fig 2k left: single-feature AUROC
ua <- avg.by(function(r) r$final_uni_auc, 'AUC');
if (!is.null(ua)) {
    ua <- ua[order(ua$AUC), ];
    write.csv(ua, file.path(output.dir, 'innerselect_uni_auc_foldmean.csv'),
              row.names = FALSE);
    }

# Supp Fig 12c: per-category AUROC drop
ba <- avg.by(function(r) r$final_block_ablation,
             c('AUC', 'AUC_Drop', 'Removed_Count'), key = 'Block');
if (!is.null(ba)) {
    write.csv(ba, file.path(output.dir, 'innerselect_block_ablation_foldmean.csv'),
              row.names = FALSE);
    }

# Fig 2k "Full Model" point: mean and SD of the outer-test AUROC.
write.csv(data.frame(
    Metric = c('mean_outer_test_auc', 'sd_outer_test_auc', 'n_folds',
               'mean_inner_holdout_auc', 'mean_n_features'),
    Value = c(mean(summary.df$final_auc), sd(summary.df$final_auc),
              nrow(summary.df), mean(summary.df$inner_holdout_auc),
              mean(summary.df$final_n_features))),
    file.path(output.dir, 'innerselect_headline.csv'), row.names = FALSE);

append.log('Completed inner-selection run.');
append.log(sprintf('Mean outer-test AUC: %.4f', mean(summary.df$final_auc)));
append.log(sprintf('SD outer-test AUC: %.4f', sd(summary.df$final_auc)));
append.log(sprintf('Mean inner-holdout AUC: %.4f', mean(summary.df$inner_holdout_auc)));
message(sprintf('Completed. Mean outer-test AUC: %.4f (SD %.4f)',
                mean(summary.df$final_auc), sd(summary.df$final_auc)));
