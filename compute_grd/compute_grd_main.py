# Run all processes to compute GRD data
import compute_load_tws
import compute_grd_tws
import compute_grd_glacier_regions
import compute_masks
import save_love_npy
import compute_grd_ensemble
import compute_grace_annual

save_love_npy.main()
compute_masks.main()
compute_load_tws.main()
compute_grd_tws.main()
compute_grd_glacier_regions.main()
compute_grd_ensemble.main()
compute_grace_annual.main()
