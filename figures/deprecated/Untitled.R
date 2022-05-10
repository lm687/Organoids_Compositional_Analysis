readRDS("../../../../GlobalDA/data/assessing_models_simulation/datasets//multiple_GenerationDMFE1_80_180_9_6_0_betaintercept2_betaslope2_1_dataset99.RDS")$beta
readRDS("../../../../GlobalDA/data/assessing_models_simulation/datasets//multiple_GenerationDMFE1_80_180_9_6_0_betaintercept1_betaslope1_1_dataset99.RDS")$beta
readRDS("../../../../GlobalDA/data/assessing_models_simulation/additional_files/multiple_fixed_sdRE1.RDS")


saveRDS(object = runif(2, -4, 4), file = "../../../../data/assessing_models_simulation/additional_files/multiple_fixed_betaintercept4.RDS")
saveRDS(object = runif(2, -4, 4), file = "../../../../data/assessing_models_simulation/additional_files/multiple_fixed_betaslope4.RDS")
