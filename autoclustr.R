# Automated cluster analysis of single molecule localisation data
#
# Automate the analysis, exploration and visualisation of clustering in single
# and multi-color single molecule localisation data, as described by
# Owen et al. (2010) and Rossy et al. (2014).
#
# Owen, D. M., Rentero, C., Rossy, J., Magenau, A., Williamson, D.,
# Rodriguez, M., & Gaus, K. (2010). PALM imaging and cluster analysis of protein
# heterogeneity at the cell surface. Journal of Biophotonics, 3(7), 446–454.
# http://doi.org/10.1002/jbio.200900089
#
# Rossy, J., Cohen, E., Gaus, K., & Owen, D. M. (2014). Method for co-cluster
# analysis in multichannel single-molecule localisation data. Histochemistry and
# Cell Biology, 141(6), 605–612. http://doi.org/10.1007/s00418-014-1208-z

library(ggplot2)
library(plyr)
library(reshape2)
library(spatstat)
library(doParallel)
library(RColorBrewer)
library(scales)
library(stringr)
library(grid)
library(gridExtra)
library(viridis)
library(imager)
library(supr)

# ------------------------------------------------------------------------
## Utility functions

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

accept <- function(ext) {
  function(name) {
    substring(name, nchar(name)-nchar(ext)+1, nchar(name)) == ext
  }
}

peaks <- function(series, span = 3) {
  z <- embed(series, span)
  s <- span%/%2
  v<- max.col(z) == 1 + s
  result <- c(rep(FALSE,s),v)
  result <- result[1:(length(result)-s)]
  result
}

any_objects <- function(label_img) {
  if (is.im(label_img)) {
    any(as.logical(label_img$v), na.rm = TRUE)
  } else if(is.mask(label_img)) {
    any(label_img$m, na.rm = TRUE)
  } else {
    stop(paste(label_img, " is not an im or mask."))
  }
}

aspect_ratio <- function(win) {
  if (!is.owin(win)) stop("points must be a planar point pattern")
  rect <- as.rectangle(win)
  width <- diff(rect$xrange)
  height <- diff(rect$yrange)
  height/width
}

# ------------------------------------------------------------------------
# Analysis Functions

#' Load RapidSTORM localisation files
#'
#' Loads all RapidSTORM localisations files from a specified directory into a
#' list of point patterns (\code{\link{ppp.object}}s). Each localisation file
#' becomes a unique \code{ppp} object with \code{x} and \code{y} coordinates
#' for each localisation and other localisation parameters store in \code{ppp}
#' marks.
#'
#' @param localisations_dir path pointing to a directory of RapidSTORM
#'    localisations files
#' @return list of point patterns
load_localisations <- function(localisations_dir, hdrs) {
  localisation_fns <- Filter(accept(".txt"), list.files(localisations_dir))
  sample_names <- tools::file_path_sans_ext(localisation_fns)
  ppps <- read.moleculelists(file.path(localisations_dir, localisation_fns),
                             rdr = read.moleculelist_rapidstorm, hdrs = hdrs)
  names(ppps) <- sample_names
  ppps
}


load_rois <- function(rois_dir) {
  roi_fns <- Filter(accept(".txt"), list.files(rois_dir))
  rois <- read.rois(file.path(rois_dir, roi_fns), xsize = roi_xy_pix_size,
                    ysize = roi_xy_pix_size)
  names(rois) <- roi_fns
  rois
}

apply_rois <- function(points, rois) {
  p <- mapply(function(sample, points, rois) {
    roi_indexes <- grep(sample, do.call(rbind, str_split(names(rois), "_"))[,1])
    lapply(roi_indexes, function(index, points, rois) {
      points[rois[[index]]]
    }, points, rois)
  }, names(points), points, MoreArgs=list(rois = rois), SIMPLIFY = FALSE)

  p <- if (length(p) > 1) do.call(c, p) else p[[1]]
  names(p) <- names(rois)
  mapply(
    function(name, points) {
      if(npoints(points) < 100) warning(paste("Low number of localisations in ROI: ", name, ". Are you should you configured the correct xy pixel size?", sep = ""))
    },
    names(rois),
    p)
  p
}

calculate_local_l <- function(points, rvalue, local_l_label, ...) {
  ll <- local_l(points, rvalue = rvalue, ...)
  marx <- cbind(marks(points), ll)
  if (spatstat::markformat(points) == "vector") {
    colnames(marx) <- c(spatstat::short.deparse(substitute(points)),
                        local_l_label)
  } else {
    colnames(marx) <- c(colnames(marx)[-length(colnames(marx))],
                        local_l_label)
  }
  marks(points) <- marx
  points
}

interpolate_surfaces <- function(points, matlab_path, which.marks=NULL,
                                 port=9999, trials = 60, maxTries = 60,
                                 interval = 2) {
  with_matlab(
    function(matlab)
      lapply(points,
             function(points, matlab) {
               img <- interpolate_surface(
                 matlab,
                 data = subset(points, select = which.marks),
                 xsize = cluster_heatmap_xy_pix_size,
                 ysize = cluster_heatmap_xy_pix_size
               )
               img[points$window, drop = FALSE]
             }, matlab),
    matlab_path = matlab_path,
    port = matlab_port,
    trials = 60,
    maxTries = 60,
    interval = 2)
}

## Measurements
cluster_centroids <- function(img, win = NULL) {
  if (is.null(win)) win <- as.rectangle(img)

  props <- region_props(img)

  if (length(props) < 1) return(ppp(c(), c(), window = win))

  centroids <- setNames(ldply(lapply(props, function(p) p$centroid()), data.frame),
                        c("id", "x", "y"))
  ppp(centroids$x, centroids$y, window = win)
}

cluster_pairwise_distances <- function(cluster_points) {
  if (npoints(cluster_points) < 2) return(data.frame(pwd=NA))

  pd <- pairdist(cluster_points)
  data.frame(pwd=pd[upper.tri(pd)])
}

cluster_nndist <- function(cluster_points) {
  if (npoints(cluster_points) < 2) return(data.frame(nn=NA))
  data.frame(nn=nndist(cluster_points))
}

prop_points_in_clusters <- function(points, mask) {
  if(!any_objects(mask)) return(0)

  mask <- as.mask(mask)
  points_in <- points[mask]
  npoints(points_in)/npoints(points)
}

## Plotting functions
plot_overview <- function(points) {
  ggplot(points, aes(x=x, y=y, color = frame)) +
      geom_point(size=0.1) + coord_fixed() +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      scale_color_viridis() +
      theme_linedraw(base_size = 8) + theme(axis.text = element_blank(),
                                            panel.grid = element_blank(),
                                            axis.ticks = element_blank(),
                                            axis.title = element_blank())
}

plot_lr <- function(data) {
  if(!is.fv(data)) stop("plot_lr require an fv.object")
  ggplot(data) + geom_line(aes(x = r, y = iso - r), colour = "blue") +
      geom_line(aes(x = r, y = theo - r), colour = "red") +
      ylab(expression(paste("L(", italic(r), ") - ", italic(r)))) +
      xlab(expression(italic(r))) + theme_linedraw(base_size = 12)
}

plot_points <- function(points) {
  ggplot(points, aes(x=x, y=y, colour=frame)) +
    geom_point(size=0.1) +
    geom_polygon(data=points$window, aes(x=x, y=y), fill=NA, colour="black") +
    coord_fixed() +
    scale_color_viridis() +
    theme_linedraw(base_size = 8) +
    theme(axis.text = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
}

plot_heatmap <- function(im, name=NULL, colours=jet.colors(256), limits = c(0, 200)) {
  ggplot(im) + geom_tile(aes(x, y, fill=value)) + coord_fixed() +
    scale_fill_gradientn(name = if (is.null(name)) "L-value" else name,
                         colours = colours, limits = limits, oob=squish) +
    theme_linedraw(base_size = 8) + theme(axis.text = element_blank(),
                                          panel.grid = element_blank(),
                                          axis.ticks = element_blank(),
                                          axis.title = element_blank())
}

plot_mask <- function(mask, window = NULL) {
  df <- as.data.frame(mask, drop=FALSE)
  plt <- ggplot() + geom_raster(data = df, aes(x, y, fill=inside)) + coord_fixed() +
    scale_fill_manual(values=c("white", "black")) +
    theme_linedraw(base_size = 8) + theme(axis.text = element_blank(),
                                          panel.grid = element_blank(),
                                          axis.ticks = element_blank(),
                                          axis.title = element_blank(),
                                          legend.position = "none")
  if (!is.null(window)) plt + geom_polygon(data=as.polygonal(window), aes(x, y, group=id), fill=NA, colour="black") else plt
}

plot_points_masks_bw <- function(points, mask) {
  plt <- ggplot() + geom_point(data=points, aes(x, y), size=0.1) + geom_polygon(data = as.polygonal(points$window), aes(x, y, group=id), fill=NA, colour="black") + coord_fixed() + theme_linedraw() + theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank())

  if (any_objects(mask))
    plt + geom_polygon(data=as.polygonal(mask), aes(x, y, group=id), fill=NA, colour="red")
  else plt
}

plot_points_masks_colour <- function(points, mask, which.marks = NULL, colours=rev(brewer.pal(11, "Spectral"))) {
  if (markformat(points) == "dataframe") {
    if (is.null(which.marks)) stop("which.marks must be specified when markformat is data.frame")
    points <- subset(points, select = which.marks)
  }
  ggplot() + geom_point(data=points, aes(x, y, colour=marks), size=0.1) + geom_polygon(data = points$window, aes(x, y, group=id), fill=NA, colour="black") + geom_polygon(data=as.polygonal(mask), aes(x, y, group=id), fill=NA, colour="red") + coord_fixed() + scale_colour_gradientn(name = which.marks, colours = colours, oob=squish) + theme_linedraw() + theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank())
}

plot_pwd_histogram <- function(pwd) {
  ggplot(pwd, aes(x=pwd)) + geom_histogram(binwidth=50) + theme_linedraw(base_size = 12) + xlab("Pair-wise distances") +
    ylab("Count")
}

plot_nn_histogram <- function(nn) {
  ggplot(nn, aes(x=nn)) + geom_histogram(binwidth=50) + theme_linedraw(base_size = 12) + xlab("Nearest-neighbour distances") +
    ylab("Count")
}

# ------------------------------------------------------------------------------
# Analysis
process_dir <- function(input_dir, output_dir_name,
                        hdrs = c("x", "y", "frame", "amplitude", "chisq", "bkgd"),
                        roi_xy_pix_size = 5,
                        use_global_k_peak = FALSE, global_k_correction = 0.75,
                        local_l_radius = 50, cluster_heatmap_xy_pix_size = 5,
                        min_degree_clustering = 3, split_clumped_clusters = TRUE,
                        split_min_r_correction_factor = 1, split_min_threshold_correction_factor = 1,
                        exclude_small_objects = TRUE, min_object_size = 1000,
                        matlab_path = NULL, matlab_port = 9999) {
  message(paste("Processing:", input_dir))
  localisations_dir <- file.path(input_dir, "Localizations")
  rois_dir <- file.path(input_dir, "ROIs")

  output_dir <- file.path(input_dir, output_dir_name)
  if(!dir.exists(output_dir)) dir.create(output_dir)

  # Find localisation files and read them in as PPPs
  pb <- txtProgressBar(style = 3)
  setTxtProgressBar(pb, 0.1, title = "Loading localisations")
  ppps <- load_localisations(localisations_dir, hdrs)

  # Find ROIs and read them in as polygon owins
  setTxtProgressBar(pb, 0.15, title = "Loading ROIs")
  rois <- load_rois(rois_dir)

  # Match PPPs and ROIs, then apply the ROIs to each PPP.
  setTxtProgressBar(pb, 0.2, title = "Applying ROIs")
  points <- apply_rois(ppps, rois)

  # Calculate L-function and peak r values
  setTxtProgressBar(pb, 0.3, title = "Calculating global L-function")
  points_lest <- lapply(points, Lest, correction="isotropic", rmax=500)

  points_lest_peak <- lapply(points_lest, peak_r, "iso")

  # Setup for local Ripley's K analysis
  ## Parameters
  if (use_global_k_peak) {
    local_l_radius <- round(mean(unlist(points_lest_peak), na.rm = TRUE) * global_k_correction, digits = 0)
  }
  local_l_label <- paste("L", local_l_radius, sep = "")

  # Calculate localL
  setTxtProgressBar(pb, 0.5, title = "Calculating local L-function")
  ## Start cluster
  points <- lapply(points, calculate_local_l, rvalue = local_l_radius, local_l_label = local_l_label,
                   correction = "isotropic")

  # Interpolate surface using matlab
  setTxtProgressBar(pb, 0.6, title = "Generating cluster heatmap")
  imgs <- interpolate_surfaces(points = points, which.marks = local_l_label,
                               matlab_path = matlab_path, port = matlab_port)

  # Segment clusters by thresholding L50
  setTxtProgressBar(pb, 0.75, title = "Threshold cluster heatmap")
  cluster_heatmap_threshold <- sqrt(min_degree_clustering) * local_l_radius
  cluster_masks <- lapply(imgs, levelset, thresh = cluster_heatmap_threshold,
                  compare = ">=")

  setTxtProgressBar(pb, 0.8, title = "Label clusters")

  cluster_labels <- if (split_clumped_clusters) {
    split_min_r <- local_l_radius * split_min_r_correction_factor
    split_min_threshold <- cluster_heatmap_threshold * split_min_threshold_correction_factor
    lpeaks <- lapply(mapply("[.ppp", points, cluster_masks, SIMPLIFY = FALSE),
                     peak_local_max, min_distance = split_min_r,
                     threshold = split_min_threshold,
                     which.marks = local_l_label)
    mapply(watershed, mask = cluster_masks, seeds = lpeaks,
           fill_lines = FALSE, SIMPLIFY = FALSE)
  } else {
    lapply(cluster_masks, function(m) if (any_objects(m)) connected(m) else as.im(m))
  }

  if (exclude_small_objects)
    cluster_labels <- lapply(cluster_labels, remove_small_objects,
                             size = min_object_size, output_mask = FALSE)

  setTxtProgressBar(pb, 0.85, title = "Calculate and output cluster statistics")

  cluster_points <- Map(function (lbl_img, points) cluster_centroids(lbl_img, points$window),
                        cluster_labels, points)
  cluster_density <- lapply(cluster_points, function (p, s) intensity(p) * s, 1000000)

  cluster_in_prop <- Map(prop_points_in_clusters, points, cluster_labels)
  cluster_pwd <- lapply(cluster_points, cluster_pairwise_distances)
  cluster_nndist <- lapply(cluster_points, cluster_nndist)

  pwd_hist_df <- ldply(cluster_pwd, as.data.frame)
  write.csv(pwd_hist_df, file = file.path(output_dir, "cluster_pwd.csv"))

  nn_hist_df <- ldply(cluster_nndist, as.data.frame)
  write.csv(nn_hist_df, file = file.path(output_dir, "cluster_nn.csv"))

  ## Peak r
  peak_r_df <- setNames(ldply(points_lest_peak, data.frame),
                        c("Sample", "Peak r (nm)"))
  write.csv(peak_r_df, file = file.path(output_dir, "peak_r.csv"))

  ## Cluster density
  cluster_density <- setNames(ldply(cluster_density, data.frame),
                              c("Sample", "Cluster density"))
  write.csv(cluster_density, file = file.path(output_dir, "cluster_density.csv"))

  ## Cluster pair-wise distances
  mean_cluster_pwd <- setNames(ldply(cluster_pwd, summarise, mean=mean(pwd)),
                               c("Sample", "Mean PWD"))
  write.csv(mean_cluster_pwd, file = file.path(output_dir, "mean_cluster_pwd.csv"))

  ## Cluster nearest neighbour distances
  mean_cluster_nndist <- setNames(ldply(cluster_nndist, summarise, mean=mean(nn)),
                               c("Sample", "Mean NN Dist"))
  write.csv(mean_cluster_nndist, file = file.path(output_dir, "mean_cluster_nndist.csv"))

  ## Proportion of Localisations inside clusters
  localisations_in_clusters <- setNames(ldply(cluster_in_prop, data.frame), c("Sample", "Propotion localisations"))
  write.csv(localisations_in_clusters, file = file.path(output_dir, "prop_localisations_in_clusters.csv"))


  # Calculate cluster size and diameter
  cluster_labels_edges <- Map(function(labels, points) clear_border(labels, points$window),
                              cluster_labels, points)

  mean_cluster_areas <- lapply(cluster_labels_edges,
                               function (c) {
                                 if (!any_objects(c)) return(NA)
                                 props <- region_props(c)
                                 mean(sapply(props, function(cp) cp$area()))
                               })
  mean_cluster_diameters <- lapply(cluster_labels_edges,
                                   function (c) {
                                     if (!any_objects(c)) return(NA)
                                     props <- region_props(c)
                                     mean(sapply(props, function(cp) cp$diameter()))
                                   })

  mean_cluster_area <- setNames(ldply(mean_cluster_areas, data.frame),
                                c("Sample", "Mean cluster area"))
  write.csv(mean_cluster_area, file = file.path(output_dir, "mean_cluster_area.csv"))

  mean_cluster_diameter <- setNames(ldply(mean_cluster_diameters, data.frame),
                                    c("Sample", "Mean cluster diameter"))
  write.csv(mean_cluster_diameter, file = file.path(output_dir, "mean_cluster_diameter.csv"))

  # Output plots
  if (!file.exists(file.path(output_dir, "reconstructed")))
    dir.create(file.path(output_dir, "reconstructed"))
  Map(function(n, p) ggsave(paste("Reconstructed ", n, ".png", sep=""), plot_overview(p),
                            path = file.path(output_dir, "reconstructed"),
                            width = 3+1, height = max(3 * aspect_ratio(p$window), 3),
                            units="in", dpi = 300), names(ppps), ppps)

  if (!file.exists(file.path(output_dir, "points")))
    dir.create(file.path(output_dir, "points"))
  Map(function(n, p) ggsave(paste(tools::file_path_sans_ext(n), "_points.pdf", sep=""), plot_points(p),
                            path = file.path(output_dir, "points"),
                            width = 3+1, height = max(3 * aspect_ratio(p$window), 2),
                            units="in", dpi = 300), names(points), points)

  if (!file.exists(file.path(output_dir, "L(r)-r")))
    dir.create(file.path(output_dir, "L(r)-r"))
  Map(function(n, p) ggsave(paste(tools::file_path_sans_ext(n), "_lr.pdf", sep=""), plot_lr(p),
                            path = file.path(output_dir, "L(r)-r"),
                            width = 3+1, height = 2,
                            units="in", dpi = 300), names(points_lest), points_lest)

  if (!file.exists(file.path(output_dir, "heatmaps")))
    dir.create(file.path(output_dir, "heatmaps"))
  Map(function(n, im) ggsave(paste(tools::file_path_sans_ext(n), "_hm.pdf", sep=""), plot_heatmap(im),
                            path = file.path(output_dir, "heatmaps"),
                            width = 3+1, height = max(3 * aspect_ratio(as.rectangle(im)), 2),
                            units="in", dpi = 300), names(imgs), imgs)

  if (!file.exists(file.path(output_dir, "masks")))
    dir.create(file.path(output_dir, "masks"))
  Map(function(n, mask) ggsave(paste(tools::file_path_sans_ext(n), "_mask.pdf", sep=""), plot_mask(mask),
                             path = file.path(output_dir, "masks"),
                             width = 3+1, height = max(3 * aspect_ratio(as.rectangle(mask)), 2),
                             units="in", dpi = 300), names(cluster_masks), cluster_masks)

  if (!file.exists(file.path(output_dir, "masks")))
    dir.create(file.path(output_dir, "masks"))
  Map(function(n, points, mask) ggsave(paste(tools::file_path_sans_ext(n), "_points_mask.pdf", sep=""),
                                       plot_points_masks_bw(points, mask),
                               path = file.path(output_dir, "masks"),
                               width = 3+1, height = max(3 * aspect_ratio(as.rectangle(mask)), 2),
                               units="in", dpi = 300), names(cluster_masks), points, cluster_masks)

  if (!file.exists(file.path(output_dir, "filtered_masks")))
    dir.create(file.path(output_dir, "filtered_masks"))
  Map(function(n, mask) ggsave(paste(tools::file_path_sans_ext(n), "_fmask.pdf", sep=""), plot_mask(as.mask(mask)),
                               path = file.path(output_dir, "filtered_masks"),
                               width = 3+1, height = max(3 * aspect_ratio(as.rectangle(mask)), 2),
                               units="in", dpi = 300), names(cluster_labels), cluster_labels)

  if (!file.exists(file.path(output_dir, "filtered_masks")))
    dir.create(file.path(output_dir, "filtered_masks"))
  Map(function(n, points, labels) ggsave(paste(tools::file_path_sans_ext(n), "_points_filtered_mask.pdf", sep=""),
                                       plot_points_masks_bw(points, as.mask(labels)),
                                       path = file.path(output_dir, "filtered_masks"),
                                       width = 3+1, height = max(3 * aspect_ratio(as.rectangle(labels)), 2),
                                       units="in", dpi = 300), names(cluster_masks), points, cluster_labels)

  if (!file.exists(file.path(output_dir, "edge_filtered_masks")))
    dir.create(file.path(output_dir, "edge_filtered_masks"))
  Map(function(n, mask) ggsave(paste(tools::file_path_sans_ext(n), "_efmask.pdf", sep=""), plot_mask(as.mask(mask)),
                               path = file.path(output_dir, "edge_filtered_masks"),
                               width = 3+1, height = max(3 * aspect_ratio(as.rectangle(mask)), 2),
                               units="in", dpi = 300), names(cluster_labels_edges), cluster_labels_edges)

  if (!file.exists(file.path(output_dir, "edge_filtered_masks")))
    dir.create(file.path(output_dir, "edge_filtered_masks"))
  Map(function(n, points, labels) ggsave(paste(tools::file_path_sans_ext(n), "_points_efmask.pdf", sep=""),
                                       plot_points_masks_bw(points, as.mask(labels)),
                                       path = file.path(output_dir, "edge_filtered_masks"),
                                       width = 3+1, height = max(3 * aspect_ratio(as.rectangle(labels)), 2),
                                       units="in", dpi = 300), names(cluster_labels_edges), points, cluster_labels_edges)

  if (!file.exists(file.path(output_dir, "pwd")))
    dir.create(file.path(output_dir, "pwd"))
  Map(function(n, pwd) {
    if(length(pwd) > 1 && !is.na(pwd)) {
      ggsave(paste(tools::file_path_sans_ext(n), "_pwd.pdf", sep=""), plot_pwd_histogram(pwd),
             path = file.path(output_dir, "pwd"), width = 3+1, height = 2, units="in", dpi = 300)
    }
  }, names(cluster_pwd), cluster_pwd)

  if (all(!is.na(pwd_hist_df$pwd)))
    ggsave("all_pwd.pdf", plot_pwd_histogram(pwd_hist_df),
           path = file.path(output_dir, "pwd"), width = 3+1, height = 2, units="in", dpi = 300)

  if (!file.exists(file.path(output_dir, "nn")))
    dir.create(file.path(output_dir, "nn"))
  Map(function(n, nn) {
    if (length(nn) > 1 && !is.na(nn)) {
      ggsave(paste(tools::file_path_sans_ext(n), "_nn.pdf", sep=""), plot_nn_histogram(nn),
             path = file.path(output_dir, "nn"),
             width = 3+1, height = 2,
             units="in", dpi = 300)
    }
  }, names(cluster_nndist), cluster_nndist)

  if (all(!is.na(nn_hist_df$nn)))
    ggsave("all_nn.pdf", plot_nn_histogram(nn_hist_df), path = file.path(output_dir, "nn"),
           width = 3+1, height = 2, units="in", dpi = 300)


  # Finally save R workspace
  setTxtProgressBar(pb, 0.95, title = "Saving R environment.")
  base::save(ppps, points, points_lest, points_lest_peak, imgs, cluster_masks,
             cluster_labels, cluster_labels_edges, cluster_points, cluster_density,
             cluster_in_prop, cluster_pwd, cluster_nndist, mean_cluster_areas,
             mean_cluster_diameters, lpeaks, cluster_heatmap_xy_pix_size,
             cluster_heatmap_xy_pix_size, exclude_small_objects, global_k_correction,
             local_l_radius, min_degree_clustering, min_object_size, roi_xy_pix_size,
             split_min_r, split_min_r_correction_factor, split_min_threshold,
             split_min_r_correction_factor, use_global_k_peak, file = file.path(output_dir, "output_data.RData"))
  setTxtProgressBar(pb, 1)
}




# Settings
#-------------------------------------------------------------------------
# Input Settings
project_dir <- "./test/test_data"

localisation_file_hdrs <- c("x", "y", "frame", "amplitude", "chisq", "bkgd")
# localisation_file_hdrs <- c("x", "precision_x", "y", "precision_y", "frame", "amplitude", "chisq", "bkgd")

## ROI pixel size calibration
roi_xy_pix_size <- 10

#-------------------------------------------------------------------------
# Output Settings
## Output directory name
output_dir_name <- "output"

## Local L and heatmap settings
localL_radius <- 40
use_global_k_peak <- FALSE
global_k_correction <- 0.75
min_degree_clustering <- 3
cluster_heatmap_xy_pix_size <- 5

## Split clusters that are stuck together
split_clumped_clusters <- TRUE
split_min_r_correction_factor <- 0.6
split_min_threshold_correction_factor <- 1.1

## Filter clusters
exclude_small_objects <- TRUE
min_object_size <- 1000

#------------------------------------------------------------------------
# Miscellaneous Settings
matlab_path <- "/Applications/MATLAB_R2017a.app/bin/matlab"
matlab_port <- 9999

parallel = TRUE
ncores <- if (detectCores() > 2) detectCores()-2 else 1

if (parallel) {
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
}

tryCatch(
  process_dir(
    input_dir = project_dir,
    output_dir_name = output_dir_name,
    hdrs = localisation_file_hdrs,
    roi_xy_pix_size = roi_xy_pix_size,
    use_global_k_peak = use_global_k_peak,
    local_l_radius = localL_radius,
    cluster_heatmap_xy_pix_size = cluster_heatmap_xy_pix_size,
    min_degree_clustering = min_degree_clustering,
    split_clumped_clusters = split_clumped_clusters,
    split_min_r_correction_factor = split_min_r_correction_factor,
    split_min_threshold_correction_factor = split_min_threshold_correction_factor,
    exclude_small_objects = exclude_small_objects,
    min_object_size = min_object_size,
    matlab_path = matlab_path
  ),
  finally = function() {
    if (parallel) {
      registerDoSEQ()
      stopCluster(cl)
    }
  }
)
