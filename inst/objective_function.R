catch_pdf <- function(pars) {
    gp$l50 <- pars["l50"]
    gp$l25 <- pars["ratio"] * gp$l50
    sp$M <- pars["M"]
    ext_mort(p)[] <- sp$M * w(p)^sp$d
    if (hasName(pars, "w_mat")) {
        sp$w_mat <- pars["w_mat"]
    }
    if (hasName(pars, "U")) {
        sp$w_mat25 <- sp$w_mat / 3^(1 / pars["U"])
    }
    p@species_params <- sp

    if (hasName(pars, "w_mat") || hasName(pars, "U")) {
        p <- setReproduction(p)
    }
    gear_params(p) <- gp

    p <- steadySingleSpecies(p)
    catch_dens <- initialN(p) * getFMort(p) * dw(p)
    names(catch_dens) <- model_lengths

    return(catch_dens)
}

objective_function <- function(params) {
    # Get the PDF from the model
    pdf_result <- catch_pdf(params)

    # Ensure pdf_result is a named vector
    if (is.null(names(pdf_result))) {
        stop("The result of catch_pdf(params) must be a named vector with names as pdf_lengths.")
    }

    # Extract pdf_lengths and pdf_values
    pdf_lengths <- as.numeric(names(pdf_result))
    pdf_values <- as.numeric(pdf_result)

    # Compute the negative log-likelihood
    neg_log_likelihood <- -log_likelihood_fn(pdf_lengths, pdf_values)

    # Handle cases where the log-likelihood is infinite or NaN
    if (is.infinite(neg_log_likelihood) || is.nan(neg_log_likelihood)) {
        neg_log_likelihood <- .Machine$double.xmax
    }

    return(neg_log_likelihood)
}
