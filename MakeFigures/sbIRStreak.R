data(soybean_ir)
soybean_ir <- soybean_ir
soybean_ir[,-1] <- log(soybean_ir[,-1]+1)
p = plotScatterStatic(soybean_ir, option ="point", saveFile = FALSE, pointSize = 0.1)
p