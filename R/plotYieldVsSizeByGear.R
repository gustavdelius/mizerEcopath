plotYieldVsSizeByGear <- function( model, catch, species = 1, return_df = FALSE){

    params <- validParams(model)

    spname <- valid_species_arg(params, species, error_on_empty = TRUE)

    idx <- which(params@species_params$species==spname)

    s <- params@species_params[idx,'species']
    a <- as.numeric(params@species_params[idx,"a"])
    b <- as.numeric(params@species_params[idx,"b"])
    w <- params@w
    l <- wlf(w,a,b)

    icatch <- catch |> filter( species == spname)

    gears <- unique(icatch$gear)
    glength <- length(gears)

    fmortal <- getFMortGear(params)

    df <- NULL

    for( i in gears){

        idx2 <- which(rownames(fmortal)==i)

        f_mort <- fmortal[idx2,idx,]

        catch_w <- f_mort * params@initial_n[idx,]
        catch_w <- catch_w/sum(catch_w * params@dw)
        catch_l <- catch_w * b * w/l

        df <- rbind(df, data.frame( Length=l, Catch_l=catch_l,
                                    Gear=i, Type="Estimated"))

        cind <- which(catch$gear==i)

        len <- catch$length[cind]
        catch_l <- catch$catch[cind]
        catch_l <- catch_l/sum(catch_l)

        df <- rbind(df, data.frame(Length=len, Catch_l=catch_l,
                                   Gear=i, Type = "Observed"))

    }

    pl <-   ggplot( df, aes(x = Length, y = Catch_l, color = Type, fill = Type)) +
        geom_bar(data = subset(df, Type != 'Estimated'), stat = "identity", position = "dodge", alpha = 0.6) +
        geom_line(data = subset(df, Type == 'Estimated'), linewidth = 1) + theme_bw() +
        facet_wrap( ~Gear) +
        labs( x = "Size [cm]", y = "Normalised number density [1/cm]", color = NULL, fill = NULL)

    print(pl)

    if(return_df) return(df)

}
