

source("functions/aux_functions.r")

res_ucell = run_ucell_trajectory(
  SC_obj_sub = SC_obj_sub,
  cds = cds
)

tiff("figures/hypoxia_trajectory.tiff", width = 6, height = 5, units = "in", res = 300)
res_ucell$plots$HALLMARK_HYPOXIA
dev.off()

tiff("figures/TNFa_trajectory.tiff", width = 6, height = 5, units = "in", res = 300)
res_ucell$plots$HALLMARK_TNFA_SIGNALING_VIA_NFKB
dev.off()



