# point to the stock downloaded and unzipped from stockassessment.org
load(url(paste0("https://www.stockassessment.org/datadisk/stockassessment/userdirs/user3/", name, "/run/model.RData")))
load(url(paste0("https://www.stockassessment.org/datadisk/stockassessment/userdirs/user3/", name, "/run/data.RData")))
load(url(paste0("https://www.stockassessment.org/datadisk/stockassessment/userdirs/user3/", name, "/run/forecast.RData")))
load(url(paste0("https://www.stockassessment.org/datadisk/stockassessment/userdirs/user3/", name, "/run/retro.RData")))


basename <- paste0('her27.nirs_',name,'_dataYear=',dataYear)

#Change the year range
year.range <- 1980:dataYear
maxyear <- year.range[length(year.range)]
ages <- 1:8
Fcv       = 0.16
Fphi      = 0.423
bio.years = c((maxyear-9),maxyear)
sel.years = c((maxyear-9),maxyear)

Fcv       = 0.212 # default 0.212, value from 2017: 0.231
Fphi      = 0.423 # default value

fit$data$natMor <- fit$data$natMor[ac(1980:dataYear),]

# build the FLStock obj slot by slot
stk <- FLStock(stock.n = FLQuant(as.matrix(t(ntable(fit))), dimnames=list(age=ages,year=year.range)))
stock.wt(stk) <- FLQuant(as.matrix(t(fit$data$stockMeanWeight)), dimnames=list(age=ages,year=year.range))
catch.wt(stk) <- FLQuant(as.matrix(t(drop(fit$data$catchMeanWeight))), dimnames=list(age=ages,year=year.range))
mat(stk) <- FLQuant(as.matrix(t(fit$data$propMat)), dimnames=list(age=ages,year=year.range))
harvest(stk) <- FLQuant(as.matrix(t(faytable(fit))), dimnames=list(age=ages,year=year.range))
m(stk) <- FLQuant(as.matrix(t(fit$data$natMor[rownames(fit$data$natMor) %in% year.range,])), dimnames=list(age=ages,year=year.range))
catch(stk) <- landings(stk) <- FLQuant(catchtable(fit)[,1], dimnames=list(age="all",year=year.range))
m.spwn(stk) <- FLQuant(as.matrix(t(fit$data$propM)), dimnames=list(age=ages,year=year.range))
harvest.spwn(stk) <- FLQuant(as.matrix(t(drop(fit$data$propF))), dimnames=list(age=ages,year=year.range))
catch.n(stk) <- FLQuant(as.matrix(t(caytable(fit))), dimnames=list(age=ages,year=year.range))
stk@discards[,,,,,] <- 0
stk@discards.n[,,,,,] <- 0
stk@discards.wt[,,,,,] <- 0
stk@landings.n <- stk@catch.n
stk@landings.wt <- stk@catch.wt
stk@harvest@units <- "f"
range(stk)['minfbar'] <- 2
range(stk)['maxfbar'] <- 4

# --------------------------------------------------------------------------------------------------------------  
# 1. Get estimate of Blim at breakpoint and infer Bpa
# --------------------------------------------------------------------------------------------------------------  
setwd(file.path("report"))

prefix <- basename
# fit HS S-R
srFit <- msy::eqsr_fit(stk, nsamp=2000, models="Segreg", rshift=1)

taf.png(paste0(prefix,"_SR_segReg"))
print(eqsr_plot(srFit, n=2e4, ggPlot=TRUE) +
        theme(legend.position="none"))
dev.off()

# extract Blim
Blim <- round(srFit$sr.det$b)
Blim

idx <- names(fit$sdrep$value) == "logssb"
sigmaSSB <- fit$sdrep$sd[idx][fit$data$years==max(fit$data$years)]

Bpa       <- round(Blim * exp(1.645*sigmaSSB))
Bpa

# --------------------------------------------------------------------------------------------------------------  
# 2. fit the stock recruitment models
# --------------------------------------------------------------------------------------------------------------  

SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blim, 
                                            ab$a * Blim, 
                                            ab$a * ssb))

FIT <- eqsr_fit(stk, nsamp = 2000, models = c("SegregBlim","Ricker", "Bevholt"), rshift=1)

taf.png(paste0(prefix,"_SR_all_models"))
eqsr_plot(FIT, n=2e4) +
  theme(legend.position="none")
dev.off()

# --------------------------------------------------------------------------------------------------------------  
# 3. Get Flim and thereby Fpa. Run EqSim with no MSY Btrigger (i.e. run EqSim with Btrigger=0), and Fcv=Fphi=0
# --------------------------------------------------------------------------------------------------------------  
SIM1 <- eqsim_run(FIT,
                  bio.years        = bio.years,
                  bio.const        = FALSE,
                  sel.years        = sel.years,
                  sel.const        = FALSE,
                  recruitment.trim = c(3, -3),
                  Fcv              = 0,
                  Fphi             = 0,
                  Blim             = Blim,
                  Bpa              = Bpa,
                  Btrigger         = 0,
                  Fscan            = seq(0,1,by=0.01),#seq(0,1.5,len=40),### 0 to 1 in steps in 0.01
                  verbose          = TRUE,
                  extreme.trim     = c(0.05,0.95))

Flim      <- SIM1$Refs2["catF","F50"]

# --------------------------------------------------------------------------------------------------------------  
# 4. Run EqSim with assessment error but no MSY Btrigger (i.e. run EqSim with Btrigger=0) to get initial FMSY
# --------------------------------------------------------------------------------------------------------------  
SIM2 <- eqsim_run(FIT,
                  bio.years = bio.years,
                  bio.const = FALSE,
                  sel.years = sel.years,
                  sel.const = FALSE,
                  recruitment.trim = c(3, -3),
                  Fcv       = Fcv,            
                  Fphi      = Fphi,            
                  Blim      = Blim,
                  Bpa       = Bpa,
                  Btrigger  = 0,
                  Fscan     = seq(0,1,by=0.01),#seq(0,1.5,len=40),### 0 to 1 in steps in 0.01
                  verbose          = TRUE,
                  extreme.trim     = c(0.05,0.95))

Fmsy      <- SIM2$Refs2["lanF","medianMSY"]

MSYBtrigger <- SIM2$Refs2["catB","F05"]
MSYBtrigger <- round(MSYBtrigger) # rounding

# --------------------------------------------------------------------------------------------------------------  
# 5. Check if FMSY is precautionary, so do a scan on Fp05. If Fmsy is larger than Fp05, reduce to Fp05
# --------------------------------------------------------------------------------------------------------------  
SIM3 <- eqsim_run(FIT,
                  bio.years = bio.years,
                  bio.const = FALSE,
                  sel.years = sel.years,
                  sel.const = FALSE,
                  recruitment.trim = c(3, -3),
                  Fcv       = Fcv,
                  Fphi      = Fphi,
                  Blim      = Blim,
                  Bpa       = Bpa,
                  Btrigger  = MSYBtrigger,
                  Fscan     = seq(0,1,by=0.01),#seq(0,1.5,len=40),##
                  verbose   = TRUE,
                  extreme.trim=c(0.05,0.95))

# If the precautionary criterion (FMSY < Fp.05) evaluated is not met, then FMSY should be reduced to  Fp.05. 
Fp05      <- SIM3$Refs2["catF","F05"]
#DM: define new Fpa here
Fpa <- Fp05
#DM: if Fpa > Flim, then Flim will be undefined
if (Fpa>Flim) Flim <- NA

propFmsy  <- subset(SIM3$pProfile, round(Ftarget, 2)==round(Fmsy,2) & variable=="Blim")$value
if (Fmsy > Fp05) {Fmsy <- Fp05}

if (Flim<0.2) Flim   <- round(Flim,3) else Flim   <- round(Flim,2)
if (Fpa<0.2) Fpa     <- round(Fpa,3) else Fpa   <- round(Fpa,2)
if (Fmsy<0.2) Fmsy   <- round(Fmsy,3) else Fmsy   <- round(Fmsy,2)

#Flim   <- round(Flim,2)
#Fpa    <- round(Fpa,2)
#Fmsy   <- round(Fmsy, 2)
refpts <- data.frame(Flim       = Flim,
                     Fpa        = Fpa,
                     Fmsy       = Fmsy,
                     Fp05       = Fp05,
                     Blim       = Blim,
                     Bpa        = Bpa,
                     MSYBtrigger= MSYBtrigger,
                     Fcv        = Fcv,
                     Fphi       = Fphi)

# print(refpts)
pander::pandoc.table(refpts, 
                     style        = "simple",
                     split.tables = 200, 
                     split.cells  = c(rep(7,10)),
                     justify      = "right",
                     missing      =" ",
                     big.mark     = '', 
                     round        = c(2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0))

write.taf(refpts,file = file.path(paste0(basename,'_output.csv')))

save(stk,srFit, FIT, 
     SIM1, SIM2, SIM3,
     refpts, 
     file=file.path(paste0(basename,"_refpoints.RData")))

setwd('..')


# -------------------------------
# plots
# -------------------------------
# taf.png(paste0(prefix,"_yield"))
# print(eqsim_plot_range(SIM3, type="mean"))
# dev.off()
# 
# taf.png(paste0(prefix,"_summary"))
# print(eqsim_plot(SIM2))
# dev.off()
# 
# # -------------------------------
# ssbOut <- ssbtable(fit) # at spawning time (Aut)
# recOut <- rectable(fit) # Aut spw, recr is 1-wr
# 
# ssbOut <- ssbOut %>%
#     as.data.frame() %>%
#     mutate(Year=as.numeric(fit$data$years))
# recOut <- recOut %>%
#     as.data.frame() %>%
#     mutate(Year=as.numeric(fit$data$years) - 2) # -2 to get the right yearclass 
# 
# tmp <- inner_join(ssbOut %>% select(Year,SSB=Estimate),
#                   recOut %>% select(Year,Rec=Estimate))
# ggplot(tmp) +
#     geom_text(aes(SSB,Rec,label=substring(Year,3,4)))
