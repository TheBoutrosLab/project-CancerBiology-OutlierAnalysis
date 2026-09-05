### HISTORY #####################################################################
# This script visualizes the XGBoost model that predicts XEG status from gene,
# tumour and somatic features.
# Date: 2026-08-05

### DESCRIPTION #################################################################
# This script visualizes the feature importance of the XGBoost model predicting
# XEG status. The cross-validation summaries are loaded from the input data, and
# the single-feature AUROC, mean absolute SHAP value and AUROC drop of the
# selected features are combined into a multipanel plot.

### PREAMBLE ####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

stable.feature <- xeg.model.feature.frequency$Variable[
    xeg.model.feature.frequency$Selected_Fraction >= 0.5
    ];

full.model.auroc <- xeg.model.summary$Value[
    xeg.model.summary$Metric == 'mean_outer_test_auc'
    ];

feature.name <- c(
    'CNA_raw' = 'CNA',
    'subtype' = 'Subtype',
    'GC_content' = 'GC Content',
    'Meth_z' = 'DNA Methylation',
    'chromosome' = 'Chromosome',
    'gene_length' = 'Gene Length',
    'exon_num' = 'Exon Number',
    'driver_cna_del_FANCA' = 'FANCA Del'
    );

# Feature categories and their colours
feature.category <- c(
    'gene_length' = 'Gene Intrinsic',
    'exon_num' = 'Gene Intrinsic',
    'GC_content' = 'Gene Intrinsic',
    'chromosome' = 'Gene Intrinsic',
    'CNA_raw' = 'Gene Extrinsic',
    'Meth_z' = 'Gene Extrinsic',
    'subtype' = 'Tumour',
    'driver_cna_del_FANCA' = 'Somatic'
    );

category.colour <- c(
    'Gene Intrinsic' = '#7DA84F',
    'Gene Extrinsic' = '#4E7FA5',
    'Tumour' = '#C0504D',
    'Somatic' = '#F0CE63'
    );

univariate.auc <- xeg.model.univariate.auc[
    xeg.model.univariate.auc$Variable %in% stable.feature,
    ];
univariate.auc <- univariate.auc[order(univariate.auc$AUC), ];

feature.order <- univariate.auc$Variable;
feature.label <- feature.name[feature.order];

feature.evaluation <- xeg.model.feature.evaluation[
    match(feature.order, xeg.model.feature.evaluation$Variable),
    ];

univariate.auc$label <- factor(feature.label, levels = feature.label);
feature.evaluation$label <- factor(feature.label, levels = feature.label);

### COVARIATE STRIP #############################################################
category.index <- as.numeric(factor(
    feature.category[feature.order],
    levels = names(category.colour)
    ));

category.strip <- create.heatmap(
    x = matrix(category.index, nrow = length(feature.order), ncol = 1),
    clustering.method = 'none',
    colour.scheme = unname(category.colour),
    total.colours = length(category.colour) + 1,
    at = seq(0.5, length(category.colour) + 0.5, 1),
    yaxis.lab = feature.label,
    yat = seq_along(feature.order),
    xaxis.lab = '',
    ylab.label = '',
    xlab.label = '',
    yaxis.tck = c(0.2, 0),
    xaxis.tck = 0,
    yaxis.cex = 0.9,
    xaxis.cex = 0,
    yaxis.rot = 0,
    print.colour.key = FALSE,
    filename = NULL
    );

### AUROC PANEL #################################################################
auroc.limits <- c(0.45, max(c(univariate.auc$AUC, full.model.auroc)) + 0.05);

full.model.point <- data.frame(
    label = factor('Full Model', levels = 'Full Model'),
    AUC = full.model.auroc,
    ci.min = full.model.auroc,
    ci.max = full.model.auroc
    );

full.model.segplot <- create.segplot(
    formula = label ~ ci.min + ci.max,
    data = full.model.point,
    centers = full.model.point$AUC,
    main = '',
    xlab.label = '',
    ylab.label = '',
    xlab.cex = 0,
    ylab.cex = 0,
    yaxis.cex = 0,
    xaxis.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xlimits = auroc.limits,
    segments.col = 'firebrick3',
    symbol.cex = 1.2,
    abline.v = 0.5,
    abline.lty = 3,
    abline.col = 'black',
    disable.factor.sorting = TRUE,
    filename = NULL
    );

univariate.auc$ci.min <- univariate.auc$AUC;
univariate.auc$ci.max <- univariate.auc$AUC;

stripe.bottom <- seq(0.5, nrow(univariate.auc) - 1.5, 2);
stripe.top <- seq(1.5, nrow(univariate.auc) - 0.5, 2);

univariate.segplot <- create.segplot(
    formula = label ~ ci.min + ci.max,
    data = univariate.auc,
    centers = univariate.auc$AUC,
    main = '',
    xlab.label = expression('AUROC'),
    ylab.label = '',
    xlab.cex = 1.2,
    ylab.cex = 0,
    yaxis.cex = 0,
    xaxis.cex = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    xlimits = auroc.limits,
    segments.col = 'grey50',
    symbol.cex = 1.0,
    abline.v = 0.5,
    abline.lty = 3,
    abline.col = 'black',
    add.rectangle = TRUE,
    xleft.rectangle = 0,
    xright.rectangle = 1,
    ybottom.rectangle = stripe.bottom,
    ytop.rectangle = stripe.top,
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    disable.factor.sorting = TRUE,
    filename = NULL
    );

### SHAP PANEL ##################################################################
shap.barplot <- create.barplot(
    formula = label ~ SHAP_sum,
    data = feature.evaluation,
    main = NULL,
    xlab.label = expression('Mean |SHAP|'),
    ylab.label = '',
    xlab.cex = 1.2,
    ylab.cex = 0,
    yaxis.cex = 0,
    xaxis.cex = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0, 0),
    xaxis.tck = c(0.2, 0),
    xlimits = c(-0.02, max(feature.evaluation$SHAP_sum) * 1.15),
    col = '#4C7A4F',
    plot.horizontal = TRUE,
    filename = NULL
    );

### ABLATION PANEL ##############################################################
ablation.barplot <- create.barplot(
    formula = label ~ Ablation_Drop,
    data = feature.evaluation,
    main = NULL,
    xlab.label = expression('AUROC Drop'),
    ylab.label = '',
    xlab.cex = 1.2,
    ylab.cex = 0,
    yaxis.cex = 0,
    xaxis.cex = 1,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0, 0),
    xaxis.tck = c(0.2, 0),
    xlimits = c(
        min(c(feature.evaluation$Ablation_Drop, 0)) * 1.3 - 0.0013,
        max(feature.evaluation$Ablation_Drop) * 1.15
        ),
    col = '#A23B2E',
    plot.horizontal = TRUE,
    filename = NULL
    );

### MULTIPANEL ##################################################################
category.legend <- legend.grob(
    legends = list(
        legend = list(
            colours = unname(category.colour),
            labels = names(category.colour),
            title = 'Feature'
            )
        ),
    label.cex = 0.9,
    title.cex = 1.0,
    size = 2,
    title.just = 'left'
    );

model.multipanel <- create.multipanelplot(
    plot.objects = list(
        full.model.segplot,
        category.strip,
        univariate.segplot,
        shap.barplot,
        ablation.barplot
        ),
    filename = NULL,
    layout.skip = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
    layout.width = 4,
    layout.height = 2,
    plot.objects.width = c(0.29, 0.30, 0.24, 0.24),
    plot.objects.height = c(0.185, 1),
    legend = list(right = list(fun = category.legend)),
    main = '',
    y.spacing = -1,
    x.spacing = -0.2
    );

save.outlier.figure(
    model.multipanel,
    c('Figure2k', 'xgboost', 'feature', 'importance', 'multipanel'),
    width = 9,
    height = 5.5
    );

save.session.profile(file.path('output', 'Figure2k.txt'));
