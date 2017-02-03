library("humarray")
rd <- rranges(9,GRanges=FALSE, fakeids=TRUE)
rd[["fakecol"]] <- sample(nrow(rd))
rd[["rs.id"]] <- paste0("rs",sample(10000,9))
ranged.to.data.frame(rd)
ranged.to.data.frame(rd,,FALSE)
ranged.to.data.frame(rd,TRUE) # keep all the columns
df.to.GRanges(ranged.to.data.frame(rd,TRUE)) 
