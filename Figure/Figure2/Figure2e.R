### HISTORY #####################################################################
# This script tests whether patients carrying a copy number aberration of a
# cancer driver gene are more likely to harbour an XEG, within PAM50 subtype.
# Date: 2026-08-05

### DESCRIPTION #################################################################
# This script analyses the enrichment of XEGs in patients carrying a copy number
# amplification or deletion of a cancer driver gene, within each PAM50 subtype.
# Odds ratios are calculated for each of METABRIC, TCGA-BRCA and ICGC BRCA-EU and
# combined by meta-analysis. The significantly enriched driver genes are shown as
# dot maps.

### PREAMBLE ####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);
library(metafor);

# Source the helper library
library(outlierAnalysisSupport);

### DATA PREPARATION ############################################################
attach(get.outlier.data.path());

load.multiple.computed.variables(c(
    'icgc.cnv.chr.new.gis.fpkm.order.match'
    ));

# Analysis settings
amplification.cut <- 1;
deletion.cut <- -1;
prefilter.proportion <- 0.05;
fdr.cutoff <- 0.10;

# Cancer driver genes recurrently amplified or deleted in breast cancer
cna.drivers.amp <- c(
    'AKT2', 'ALK', 'CCNE1', 'DROSHA', 'EGFR', 'ERBB2', 'ERG', 'FLT4', 'HGF',
    'JUN', 'KAT6A', 'LMO1', 'MDM2', 'MDM4', 'MITF', 'MYC', 'MYCL', 'MYCN',
    'NKX2-1', 'WHSC1L1', 'NTRK1', 'PPM1D', 'RAF1', 'REL', 'RICTOR', 'SOX2',
    'YAP1'
    );

cna.drivers.del <- c(
    'ACVR1B', 'AMER1', 'APC', 'ARID1A', 'ATM', 'AXIN1', 'BIRC3', 'BRCA1',
    'BRCA2', 'CDKN2A', 'CDKN2C', 'FANCA', 'FANCC', 'FANCD2', 'FBXO11', 'FBXW7',
    'GPC3', 'IKZF1', 'IKZF3', 'KDM6A', 'LATS2', 'LRP1B', 'MAP2K4', 'MEN1',
    'MLH1', 'MSH2', 'MUTYH', 'NF1', 'NF2', 'PAX5', 'PBRM1', 'PRDM1', 'PTEN',
    'RAD50', 'RB1', 'SMAD4', 'SMARCB1', 'STK11', 'SUFU', 'TNFAIP3', 'TP53',
    'TSC1', 'TSC2', 'VHL', 'WT1'
    );

all.drivers <- c(cna.drivers.amp, cna.drivers.del);

# PAM50 subtypes, in the order used by subtype.total.outlier.num
subtype.key <- c('basal', 'her2', 'luma', 'lumb', 'normal');
subtype.label <- c('Basal', 'Her2', 'LumA', 'LumB', 'Normal');

### XEG COUNT PER PATIENT #######################################################
# Number of XEGs in each patient
xeg.per.patient <- c(
    apply(outlier.patient.tag.01.brca, 2, sum),
    apply(outlier.patient.tag.01.meta, 2, sum),
    apply(outlier.patient.tag.01.icgc, 2, sum)
    );

build.subtype.xeg <- function(subtype.id) {
    subtype.patient <- subtype.total.outlier.num[
        subtype.total.outlier.num$subtype == subtype.id,
        ,
        drop = FALSE
        ];
    subtype.xeg <- xeg.per.patient[rownames(na.omit(subtype.patient))];

    brca.xeg <- subtype.xeg[substr(names(subtype.xeg), 1, 4) == 'TCGA'];
    # TCGA-BRCA copy number barcodes carry no trailing A
    names(brca.xeg) <- gsub('A$', '', names(brca.xeg));

    return(list(
        meta = subtype.xeg[substr(names(subtype.xeg), 1, 2) == 'MB'],
        brca = brca.xeg,
        icgc = subtype.xeg[substr(names(subtype.xeg), 1, 2) == 'PR']
        ));
    }

### DRIVER COPY NUMBER STATUS ###################################################
extract.driver.cna <- function(cna.data, gene.symbol, drivers) {
    driver.cna <- NULL;
    driver.kept <- NULL;
    for (gene in drivers) {
        index <- which(gene.symbol == gene);
        if (length(index) > 0) {
            driver.cna <- rbind(driver.cna, cna.data[index[1], ]);
            driver.kept <- c(driver.kept, gene);
            }
        }
    driver.cna <- data.frame(driver.cna, check.names = FALSE);
    rownames(driver.cna) <- driver.kept;
    return(driver.cna);
    }

# METABRIC and TCGA-BRCA GISTIC tables carry two annotation columns
annotation.column <- c('Hugo_Symbol', 'Entrez_Gene_Id');

driver.cna.meta <- extract.driver.cna(
    meta.cnv.chr.new.gis,
    meta.cnv.chr.new.gis$Hugo_Symbol,
    all.drivers
    );
driver.cna.meta <- driver.cna.meta[
    ,
    !colnames(driver.cna.meta) %in% annotation.column
    ];

driver.cna.brca <- extract.driver.cna(
    brca.cnv.chr.new.gis,
    brca.cnv.chr.new.gis$Hugo_Symbol,
    all.drivers
    );
driver.cna.brca <- driver.cna.brca[
    ,
    !colnames(driver.cna.brca) %in% annotation.column
    ];

# The ICGC table uses two older gene symbols
icgc.gene.symbol <- icgc.cnv.chr.new.gis.raw$Gene.Symbol;
icgc.gene.symbol[icgc.gene.symbol == 'MYCL1'] <- 'MYCL';
icgc.gene.symbol[icgc.gene.symbol == 'MLL3'] <- 'KMT2C';

driver.cna.icgc <- extract.driver.cna(
    icgc.cnv.chr.new.gis.fpkm.order.match,
    icgc.gene.symbol,
    all.drivers
    );

define.cna.status <- function(cna.row, mode) {
    values <- suppressWarnings(as.numeric(cna.row));
    if ('amp' == mode) {
        status <- ifelse(values >= amplification.cut, 1, 0);
        } else {
        status <- ifelse(values <= deletion.cut, 1, 0);
        }
    names(status) <- names(cna.row);
    return(status);
    }

driver.status <- list();
for (gene in all.drivers) {
    mode <- if (gene %in% cna.drivers.amp) 'amp' else 'del';
    driver.status[[gene]] <- list(
        meta = if (gene %in% rownames(driver.cna.meta)) define.cna.status(driver.cna.meta[gene, ], mode) else NULL,
        brca = if (gene %in% rownames(driver.cna.brca)) define.cna.status(driver.cna.brca[gene, ], mode) else NULL,
        icgc = if (gene %in% rownames(driver.cna.icgc)) define.cna.status(driver.cna.icgc[gene, ], mode) else NULL
        );
    }

### ENRICHMENT ##################################################################
# Odds ratio and p-value from the 2 x 2 table of driver status against XEG status
compute.driver.enrichment <- function(status, xeg.count) {
    if (is.null(status)) {
        return(NULL);
        }
    common <- intersect(names(status), names(xeg.count));
    if (length(common) < 10) {
        return(NULL);
        }

    driver <- suppressWarnings(as.numeric(status[common]));
    xeg <- ifelse(xeg.count[common] > 0, 1, 0);
    usable <- is.finite(driver) & is.finite(xeg);
    driver <- driver[usable];
    xeg <- xeg[usable];
    if (length(driver) < 10) {
        return(NULL);
        }

    a <- sum(driver == 1 & xeg == 1) + 0.5;
    b <- sum(driver == 1 & xeg == 0) + 0.5;
    c <- sum(driver == 0 & xeg == 1) + 0.5;
    d <- sum(driver == 0 & xeg == 0) + 0.5;

    return(list(
        or = (a * d) / (b * c),
        se.log.or = sqrt(1 / a + 1 / b + 1 / c + 1 / d),
        p.value = tryCatch(
            fisher.test(matrix(c(a - 0.5, b - 0.5, c - 0.5, d - 0.5), nrow = 2, byrow = TRUE))$p.value,
            error = function(e) NA_real_
            ),
        exposed = sum(driver == 1),
        n = length(driver)
        ));
    }

# Pooled across datasets by DerSimonian-Laird random-effects meta-analysis
combine.driver.enrichment <- function(results) {
    valid <- results[!sapply(results, function(x) is.null(x) || !is.finite(x$or) || !is.finite(x$se.log.or))];
    if (length(valid) <= 1) {
        return(NULL);
        }

    pooled <- tryCatch(
        rma.uni(
            yi = sapply(valid, function(x) log(x$or)),
            sei = sapply(valid, function(x) x$se.log.or),
            method = 'DL'
            ),
        error = function(e) NULL
        );
    if (is.null(pooled)) {
        return(NULL);
        }

    return(list(
        or = as.numeric(exp(pooled$beta)),
        p.value = as.numeric(pooled$pval)
        ));
    }

dataset.key <- c('meta', 'brca', 'icgc');
dataset.name <- c('METABRIC', 'TCGA-BRCA', 'ICGC BRCA-EU');

test.subtype.drivers <- function(subtype.id, drivers) {
    subtype.xeg <- build.subtype.xeg(subtype.id);

    per.dataset <- list();
    passed <- NULL;
    for (gene in drivers) {
        per.dataset[[gene]] <- lapply(
            dataset.key,
            function(key) compute.driver.enrichment(driver.status[[gene]][[key]], subtype.xeg[[key]])
            );
        names(per.dataset[[gene]]) <- dataset.key;

        supporting <- 0;
        for (key in dataset.key) {
            result <- per.dataset[[gene]][[key]];
            if (!is.null(result) && result$exposed > length(subtype.xeg[[key]]) * prefilter.proportion) {
                supporting <- supporting + 1;
                }
            }
        if (supporting >= 2) {
            passed <- c(passed, gene);
            }
        }

    if (is.null(passed)) {
        return(NULL);
        }
    common.genes <- unique(passed);
    if (0 == length(common.genes)) {
        return(NULL);
        }

    combined <- lapply(common.genes, function(gene) combine.driver.enrichment(per.dataset[[gene]]));
    names(combined) <- common.genes;

    pooled.or <- sapply(common.genes, function(gene) combined[[gene]]$or);
    pooled.p <- sapply(common.genes, function(gene) combined[[gene]]$p.value);
    pooled.fdr <- p.adjust(pooled.p, method = 'BH');

    tested <- data.frame(
        gene = common.genes,
        odds.ratio = as.numeric(pooled.or),
        p.value = as.numeric(pooled.p),
        fdr = as.numeric(pooled.fdr),
        stringsAsFactors = FALSE
        );
    tested <- tested[order(tested$fdr), ];

    retained <- common.genes[is.finite(pooled.fdr) & pooled.fdr < fdr.cutoff];
    if (0 == length(retained)) {
        return(list(tested = tested, genes = character(0)));
        }

    odds.ratio <- matrix(
        NA_real_,
        nrow = length(dataset.key),
        ncol = length(retained),
        dimnames = list(dataset.name, retained)
        );
    p.value <- odds.ratio;
    for (i in seq_along(retained)) {
        for (j in seq_along(dataset.key)) {
            result <- per.dataset[[retained[i]]][[dataset.key[j]]];
            if (!is.null(result)) {
                odds.ratio[j, i] <- result$or;
                p.value[j, i] <- result$p.value;
                }
            }
        }

    fdr <- p.value;
    for (j in seq_len(nrow(p.value))) {
        fdr[j, ] <- p.adjust(p.value[j, ], method = 'BH');
        }

    return(list(
        tested = tested,
        genes = retained,
        odds.ratio = odds.ratio,
        fdr = fdr,
        pooled.or = pooled.or[retained],
        pooled.fdr = pooled.fdr[retained]
        ));
    }

### DOTMAP ######################################################################
spot.size.function <- function(x) {
    0.1 + 2 * abs(x);
    }

spot.colour.function <- function(x) {
    colour <- rep('white', length(x));
    colour[sign(x) == -1] <- default.colours(2, palette.type = 'dotmap')[1];
    colour[sign(x) == 1] <- default.colours(2, palette.type = 'dotmap')[2];
    return(colour);
    }

background.cutoff <- 6;
colourkey.labels.at <- seq(0, background.cutoff, by = 1);
colourkey.labels <- sapply(
    X = colourkey.labels.at,
    FUN = function(x) {
        if (x == 0) {
            return(expression('10'^'0'));
            } else if (x != background.cutoff) {
            return(as.expression(bquote('10'^-.(as.character(x)))));
            } else {
            return(as.expression(bquote('<10'^-.(as.character(x)))));
            }
        }
    );

fdr.legend <- legend.grob(
    list(
        legend = list(
            title = expression(underline('FDR')),
            continuous = TRUE,
            colours = c('white', 'black'),
            total.colours = 100,
            labels = colourkey.labels,
            cex = 0.9,
            at = seq(0, 100, length.out = length(colourkey.labels)),
            height = 3
            )
        ),
    label.cex = 1,
    title.cex = 1,
    title.just = 'left',
    title.fontface = 'plain',
    between.row = 4
    );

spot.size.key <- list(
    space = 'right',
    points = list(
        cex = spot.size.function(seq(-2, 2, 1)),
        col = spot.colour.function(seq(-2, 2, 1)),
        pch = 19
        ),
    text = list(
        lab = c('0.25', '0.5', '1', '2', '4'),
        cex = 1,
        adj = 1
        ),
    padding.text = 8
    );

create.driver.dotmap <- function(result, panel.title) {
    dot.each <- create.dotmap(
        x = log2(result$odds.ratio),
        main = panel.title,
        main.cex = 1.4,
        yaxis.lab = dataset.name,
        yaxis.cex = 1.2,
        xaxis.rot = 90,
        xaxis.cex = 0,
        xaxis.fontface = 1,
        yaxis.fontface = 1,
        yaxis.tck = c(0.2, 0),
        xaxis.tck = c(0.2, 0),
        spot.size.function = spot.size.function,
        spot.colour.function = spot.colour.function,
        legend = list(inside = list(fun = fdr.legend, x = 1.07, y = 0)),
        key = spot.size.key,
        key.top = 1,
        right.padding = 2,
        pch = 21,
        pch.border.col = 'white',
        bg.data = -log10(result$fdr),
        colourkey = FALSE,
        bg.alpha = 1,
        colour.scheme = c('white', 'black'),
        at = seq(0, background.cutoff, 0.01),
        row.colour = 'white',
        col.colour = 'white',
        row.lwd = 1.3,
        col.lwd = 1.3,
        filename = NULL
        );

    dot.all <- create.dotmap(
        x = t(log2(result$pooled.or)),
        xaxis.lab = result$genes,
        yaxis.lab = 'All patients',
        top.padding = 10,
        yaxis.cex = 1.2,
        xaxis.rot = 90,
        xaxis.cex = 1,
        xaxis.fontface = 1,
        yaxis.fontface = 1,
        yaxis.tck = c(0.2, 0),
        xaxis.tck = c(0.2, 0),
        spot.size.function = spot.size.function,
        spot.colour.function = spot.colour.function,
        key.top = 1,
        right.padding = 2,
        pch = 21,
        pch.border.col = 'white',
        bg.data = t(-log10(result$pooled.fdr)),
        colourkey = FALSE,
        bg.alpha = 1,
        colour.scheme = c('white', 'black'),
        at = seq(0, background.cutoff, 0.01),
        row.colour = 'white',
        col.colour = 'white',
        row.lwd = 1.3,
        col.lwd = 1.3,
        filename = NULL
        );

    return(create.multipanelplot(
        list(dot.each, dot.all),
        layout.height = 2,
        layout.width = 1,
        plot.objects.heights = c(10.5, 5),
        x.spacing = -1,
        y.spacing = -10,
        bottom.padding = 0,
        top.padding = 2,
        right.padding = 0
        ));
    }

### OUTPUT ######################################################################
# Every subtype and both directions are tested
tested.summary <- NULL;

for (i in seq_along(subtype.key)) {
    for (mode in c('amp', 'del')) {
        drivers <- if ('amp' == mode) cna.drivers.amp else cna.drivers.del;
        direction <- if ('amp' == mode) 'Amplification' else 'Deletion';

        result <- test.subtype.drivers(i, drivers);
        if (is.null(result)) {
            next;
            }

        tested.summary <- rbind(
            tested.summary,
            data.frame(
                subtype = subtype.label[i],
                direction = direction,
                result$tested,
                retained = result$tested$gene %in% result$genes,
                stringsAsFactors = FALSE
                )
            );

        if (0 == length(result$genes)) {
            next;
            }

        driver.dotmap <- create.driver.dotmap(
            result,
            paste0('CNA Driver Genes - ', direction, ' (', subtype.label[i], ')')
            );

        save.outlier.figure(
            driver.dotmap,
            c('Figure2e', 'drivergene', 'cna', mode, subtype.key[i], 'multipanel'),
            width = 6,
            height = 5.5
            );
        }
    }

write.csv(
    tested.summary,
    file.path(here::here('output'), 'Figure2e_driver_cna_enrichment.csv'),
    row.names = FALSE
    );

save.session.profile(file.path('output', 'Figure2e.txt'));
