#' Combine several MizerParams objects
#'
#' @description `r lifecycle::badge("experimental")`
#'
#'   Takes several \linkS4class{MizerParams} objects and fuses them into a
#'   single \linkS4class{MizerParams} object containing all species. The
#'   interaction matrix of the combined object is set to all zeros.
#'
#' @param ... Two or more \linkS4class{MizerParams} objects to combine.
#'   Alternatively a single list of \linkS4class{MizerParams} objects.
#'
#' @return A \linkS4class{MizerParams} object containing all species from all
#'   input params objects.
#'
#' @details All input params objects must have identical `w` and `w_full` slots
#'   and identical `resource_params`. No duplicate species names are allowed.
#'
#'   The `species_params`, `given_species_params` and `gear_params` data frames
#'   are combined with `rbind()`. The per-species arrays from the various slots
#'   are stacked along the species dimension. The interaction matrix is set to
#'   all zeros. All other slots (resource spectrum, rates functions, etc.) are
#'   taken from the first params object.
#'
#' @seealso [addSpecies()], [removeSpecies()]
#' @export
bindParams <- function(...) {
    # Collect params into a list ----
    params_list <- list(...)
    if (length(params_list) == 1 && is.list(params_list[[1]]) &&
        !is(params_list[[1]], "MizerParams")) {
        params_list <- params_list[[1]]
    }
    if (length(params_list) < 2) {
        stop("bindParams() requires at least two MizerParams objects.")
    }
    params_list <- lapply(params_list, validParams)

    # Check identical w and w_full ----
    w_ref <- params_list[[1]]@w
    w_full_ref <- params_list[[1]]@w_full
    for (i in seq_along(params_list)[-1]) {
        if (!isTRUE(all.equal(params_list[[i]]@w, w_ref))) {
            stop("All MizerParams objects must have identical `w` slots.")
        }
        if (!isTRUE(all.equal(params_list[[i]]@w_full, w_full_ref))) {
            stop("All MizerParams objects must have identical `w_full` slots.")
        }
    }

    # Check identical resource_params ----
    rp_ref <- params_list[[1]]@resource_params
    for (i in seq_along(params_list)[-1]) {
        if (!isTRUE(all.equal(params_list[[i]]@resource_params, rp_ref))) {
            stop("All MizerParams objects must have identical `resource_params`.")
        }
    }

    # Check no duplicate species ----
    all_species <- unlist(lapply(params_list, function(p) p@species_params$species))
    if (anyDuplicated(all_species)) {
        dupes <- unique(all_species[duplicated(all_species)])
        stop("Duplicate species found: ", toString(dupes))
    }
    no_sp <- length(all_species)
    no_w <- length(w_ref)

    # Start with first params as template ----
    p <- params_list[[1]]

    # Helper: align columns across data frames and rbind ----
    align_rbind <- function(dfs) {
        all_cols <- Reduce(union, lapply(dfs, names))
        dfs <- lapply(dfs, function(df) {
            df[setdiff(all_cols, names(df))] <- NA
            df[all_cols]
        })
        do.call(rbind, c(dfs, list(stringsAsFactors = FALSE)))
    }

    # Combine data frames ----
    p@species_params <-
        align_rbind(lapply(params_list, function(x) x@species_params))
    p@given_species_params <-
        align_rbind(lapply(params_list, function(x) x@given_species_params))
    p@gear_params <- align_rbind(lapply(params_list, function(x) x@gear_params))
    p@gear_params <- validGearParams(p@gear_params, p@species_params)

    # Combine 2D [species x w] arrays ----
    stack_rows <- function(slot_name) {
        arrays <- lapply(params_list, function(x) slot(x, slot_name))
        result <- do.call(rbind, arrays)
        names(dimnames(result)) <- names(dimnames(arrays[[1]]))
        result
    }
    p@psi           <- stack_rows("psi")
    p@maturity      <- stack_rows("maturity")
    p@initial_n     <- stack_rows("initial_n")
    p@intake_max    <- stack_rows("intake_max")
    p@search_vol    <- stack_rows("search_vol")
    p@metab         <- stack_rows("metab")
    p@mu_b          <- stack_rows("mu_b")
    p@ext_encounter <- stack_rows("ext_encounter")
    p@diffusion     <- stack_rows("diffusion")

    # Combine [species x w_full] arrays ----
    p@ft_pred_kernel_e <- stack_rows("ft_pred_kernel_e")
    p@ft_pred_kernel_p <- stack_rows("ft_pred_kernel_p")
    p@ft_mask          <- stack_rows("ft_mask")

    # Combine pred_kernel (may be 3D [species x sp_size x prey_size] or empty) ----
    has_pk <- sapply(params_list, function(x) length(dim(x@pred_kernel)) == 3)
    if (all(has_pk)) {
        sp_dims <- dim(params_list[[1]]@pred_kernel)[2:3]
        combined_pk <- array(0, dim = c(no_sp, sp_dims),
                             dimnames = c(list(sp = all_species),
                                          dimnames(params_list[[1]]@pred_kernel)[2:3]))
        row_start <- 1L
        for (pp in params_list) {
            nsp_pp <- nrow(pp@species_params)
            combined_pk[row_start:(row_start + nsp_pp - 1L), , ] <- pp@pred_kernel
            row_start <- row_start + nsp_pp
        }
        p@pred_kernel <- combined_pk
    } else if (any(has_pk)) {
        stop("Either all or none of the params objects must have a pred_kernel.")
    }

    # Combine selectivity [gear x species x w] ----
    all_gears <- unique(unlist(lapply(params_list, function(x) dimnames(x@selectivity)[[1]])))
    w_dimnames <- dimnames(params_list[[1]]@selectivity)[[3]]
    combined_sel <- array(0, dim = c(length(all_gears), no_sp, no_w),
                          dimnames = list(gear = all_gears, sp = all_species,
                                          w = w_dimnames))
    sp_start <- 1L
    for (pp in params_list) {
        nsp_pp <- nrow(pp@species_params)
        gears_pp <- dimnames(pp@selectivity)[[1]]
        if (length(gears_pp) > 0) {
            combined_sel[gears_pp, sp_start:(sp_start + nsp_pp - 1L), ] <- pp@selectivity
        }
        sp_start <- sp_start + nsp_pp
    }
    p@selectivity <- combined_sel

    # Combine catchability [gear x species] ----
    combined_catch <- array(0, dim = c(length(all_gears), no_sp),
                            dimnames = list(gear = all_gears, sp = all_species))
    sp_start <- 1L
    for (pp in params_list) {
        nsp_pp <- nrow(pp@species_params)
        gears_pp <- dimnames(pp@catchability)[[1]]
        if (length(gears_pp) > 0) {
            combined_catch[gears_pp, sp_start:(sp_start + nsp_pp - 1L)] <- pp@catchability
        }
        sp_start <- sp_start + nsp_pp
    }
    p@catchability <- combined_catch

    # Combine initial_effort (union of gears, max for shared gears) ----
    combined_effort <- setNames(numeric(length(all_gears)), all_gears)
    for (pp in params_list) {
        eff <- pp@initial_effort
        combined_effort[names(eff)] <- pmax(combined_effort[names(eff)], eff)
    }
    p@initial_effort <- combined_effort

    # Set interaction matrix to 0 ----
    p@interaction <- array(0, dim = c(no_sp, no_sp),
                           dimnames = list(predator = all_species, prey = all_species))

    # Combine per-species named vectors ----
    p@w_min_idx <- unlist(lapply(params_list, function(x) x@w_min_idx))
    p@A <- unlist(lapply(params_list, function(x) x@A))

    # Combine linecolour and linetype (preserve special entries from first params) ----
    special <- c("Resource", "Total", "Background", "Fishing", "External")
    species_lc <- unlist(lapply(params_list, function(x) {
        lc <- x@linecolour
        lc[!(names(lc) %in% special)]
    }))
    p@linecolour <- c(species_lc, p@linecolour[intersect(special, names(p@linecolour))])
    species_lt <- unlist(lapply(params_list, function(x) {
        lt <- x@linetype
        lt[!(names(lt) %in% special)]
    }))
    p@linetype <- c(species_lt, p@linetype[intersect(special, names(p@linetype))])

    p@time_modified <- lubridate::now()
    validObject(p)
    return(p)
}
