#' Calculate catch selectivity params from catch data
#'
#' This function determines the selectivity parameter for the params object from
#' the user provided scaled catch data. Catch data is rebinned into mizer bins
#' and then divided by the total abundace. When abundance reaches 0 this means
#' selectivity values will reach Inf and a warning will be issued that the
#' steady state model does not support large enough fish according to the catch
#' data.
#'
#' @param params A MizerParams object. a and b species parameters are required
#'   in the params object
#' @param catch A data frame containing at least the following columns:
#'   'species' (common names for each species), 'gear', 'length', and 'catch'
#'   (scaled to total g/m^2/year). Additional value columns are allowed.
#' @return A MizerParams object with selectivity for catch data

catchSelectivity <- function(params,
                             catch){
    # Validate inputs
    if (!inherits(params, "MizerParams")) {
        stop("`params` must be a MizerParams object.")
    }

    if (!is.data.frame(catch)) {
        stop("`catch` must be a data frame.")
    }

    if (!"species" %in% colnames(catch)){
        stop("`species` column must be in data frame.")
    }

    if (!"gear" %in% colnames(catch)){
        stop("`gear` column must be in data frame.")
    }

    library(dplyr)
    library(purrr)
    library(reshape2)
    library(tidyr)
    library(mizer)

    #Get Ni(w)*dw
    Niwdw<-sweep(initialN(params), MARGIN = 2, STATS = dw(params), FUN = "*")

    #Catch to mizer bins
    gear_list<-unique(catch$gear)
    results <- list()
    for (g in gear_list) {
        # select catch for this gear
        catch_gear <- catch %>% filter(gear == g)
        # --- 1. Prepare catch with LWR ---
        LWR <- species_params(params) %>% select(species, a, b)

        catch_LWR <- left_join(catch_gear, LWR, by = "species") %>%
            mutate(weight = a * length^b) %>%
            group_by(species) %>%
            arrange(length, .by_group = TRUE) %>%
            mutate(bin_weight_width = lead(weight) - weight) %>%
            ungroup()

        # --- 2. Mizer bin edges ---
        bin_mids <- w(params)
        bin_widths <- dw(params)
        bin_count <- length(bin_mids)
        mizer_edges <- c(bin_mids - bin_widths / 2,
                         tail(bin_mids, 1) + tail(bin_widths, 1) / 2)

        # --- 3. Add lower and upper weight bounds ---
        catch_bins <- catch_LWR %>%
            arrange(species, length) %>%
            group_by(species) %>%
            mutate(
                lower_w = weight,
                upper_w = lead(weight),
                upper_w = ifelse(is.na(upper_w) & !is.na(bin_weight_width),
                                 lower_w + bin_weight_width, upper_w)
            ) %>%
            ungroup() %>%
            filter(!is.na(lower_w) & !is.na(upper_w) & !is.na(catch))

        # --- 4. Bin splitting (using base R lapply) ---
        all_splits <- lapply(seq_len(nrow(catch_bins)), function(i) {
            row <- catch_bins[i, ]
            species <- row$species
            lower_w <- row$lower_w
            upper_w <- row$upper_w
            catch_value <- row$catch

            if (is.na(catch_value) || is.na(lower_w) || is.na(upper_w) || upper_w <= lower_w) {
                return(data.frame(species = species, bin = seq_len(bin_count), catch = rep(NA_real_, bin_count)))
            }

            overlaps <- numeric(bin_count)
            for (j in seq_len(bin_count)) {
                m_lower <- mizer_edges[j]
                m_upper <- mizer_edges[j + 1]
                overlap <- max(0, min(upper_w, m_upper) - max(lower_w, m_lower))
                if (!is.na(overlap) && overlap > 0) {
                    prop <- overlap / (upper_w - lower_w)
                    overlaps[j] <- catch_value * prop
                }
            }

            data.frame(species = species, bin = seq_len(bin_count), catch = overlaps)
        })

        catch_mizer <- do.call(rbind, all_splits)

        # --- 5. Summarise and fill missing species-bin combos ---
        catch_mizer <- catch_mizer %>%
            group_by(species, bin) %>%
            summarise(
                catch = if (all(is.na(catch))) NA_real_ else sum(catch, na.rm = TRUE),
                .groups = "drop"
            )

        # Manually complete species-bin combinations
        all_species <- unique(catch_mizer$species)
        all_bins <- seq_len(bin_count)
        full_grid <- expand.grid(species = all_species, bin = all_bins)

        catch_mizer <- full_join(full_grid, catch_mizer, by = c("species", "bin"))

        # --- 6. Add bin metadata and calculate values ---
        w_df <- data.frame(w = signif(bin_mids,3), bin = seq_len(bin_count))
        dw_df <- data.frame(dw = bin_widths, bin = seq_len(bin_count))

        catch_mizer <- catch_mizer %>%
            left_join(w_df, by = "bin") %>%
            left_join(dw_df, by = "bin") %>%
            mutate(value = catch,
                   sp=species)

        catch_mizer$gear <- g

        # store result
        results[[g]] <- catch_mizer
    }
    # combine everything into one data frame
    final_catch <- do.call(rbind, results)

    # --- 7. Build matrix
    catch_array <- with(final_catch,
                        tapply(value, list(gear, sp, w), FUN = function(x) if(length(x)) x else NA))

    # Suppose Sgiw has dimnames like this:
    # dimnames(Sgiw)[[2]] -> species names in current order
    current_species <- dimnames(catch_array)[[2]]

    # Get desired order from selectivity(params)
    desired_species <- dimnames(Niwdw)[[1]]

    # Find the matching indices
    reorder_idx <- match(desired_species, current_species)

    # Reorder Sgiw along the species dimension (2nd dimension)
    catch_array <- catch_array[, reorder_idx, , drop = FALSE]

    #get Sg,i(w)
    Sgiw <- sweep(catch_array, MARGIN = c(2,3), STATS = Niwdw, FUN = "/")

    if (any(is.infinite(Sgiw))) {
        # Find positions of infinite values
        inf_positions <- which(is.infinite(Sgiw), arr.ind = TRUE)

        # Extract species names corresponding to those positions
        inf_species <- unique(dimnames(Sgiw)[[2]][inf_positions[, "dim2"]])

        # Build message
        msg <- paste(
            "Gear is catching individuals larger than present in the model.\n",
            "Increase individuals at large sizes for the following species:", paste(inf_species, collapse = ", ")
        )

        # Show warning
        withr::with_options(
            list(warn = 1),
            warning(msg, call. = FALSE)
        )
    }
    Sgiw[is.infinite(Sgiw)] <- 0
    Sgiw[is.nan(Sgiw)] <- 0

    sp <- species_params(params)
    sp<-sp$species
    for (i in sp){
        selectivity(params)[,i,]<-Sgiw[,i,]
    }
    initial_effort(params)<-1
    gear_params(params)$catchability<-1
    ## Change external mortality rate
    fish_mort <- getFMort(params)
    new_ext_mort <- ext_mort(params) - fish_mort

    sp <- species_params(params)
    no_sp <- nrow(sp)
    # Check that fishing mortality rate is less than total mortality rate
    # for all species up to at least their largest size
    for (i in 1:no_sp) {
        if (any(new_ext_mort[i, ] < 0) &&
            (min(params@w[new_ext_mort[i, ] < 0]) <=
             max(params@w[initialN(params)[i, ] > 0]))) {
            warning("Implied negative external mortality rate for ",
                    params@species_params$species[i],
                    ". Setting it to zero. This may alter the steady state.")
        }
    }
    # Don't allow negative mortality rates
    new_ext_mort[new_ext_mort < 0] <- 0
    ext_mort(params) <- new_ext_mort

    return(params)
}
