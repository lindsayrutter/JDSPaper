load("soybean_ir_noFiltnoP3_metrics.rda")
metrics <- soybean_ir_noFiltnoP3_metrics[[1]]

streakGenes <- c("Glyma.18G057100.Wm82.a2.v1", "Glyma.18G092200.Wm82.a2.v1", "Glyma.09G164000.Wm82.a2.v1","Glyma.08G269900.Wm82.a2.v1", "Glyma.07G079300.Wm82.a2.v1")

streakMetrics <- metrics[which(metrics$ID %in% streakGenes),]

streakMetrics$FDR
# [1] 0.004365995 0.005296390 0.008645290 0.010265304 0.010453529