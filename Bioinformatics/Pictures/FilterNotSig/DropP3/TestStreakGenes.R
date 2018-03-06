load("../soybean_ir_noFilt_metrics.rda")
metrics <- soybean_ir_noFilt_metrics[[1]]

streakGenes <- c("Glyma.18G057100.Wm82.a2.v1", "Glyma.18G092200.Wm82.a2.v1", "Glyma.09G164000.Wm82.a2.v1","Glyma.08G269900.Wm82.a2.v1", "Glyma.07G079300.Wm82.a2.v1")

streakMetrics <- metrics[which(metrics$ID %in% streakGenes),]

streakMetrics$FDR
# [1] 0.5087831 0.6025199 0.8859584 0.9118598 0.9224019