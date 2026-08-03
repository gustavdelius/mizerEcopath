#include <TMB.hpp>

// Trapezoidal bin-average of the weight K of a summary integral
// \int N(w) K(w) dw, discretised as \sum_j N_j Kbar_j dw_j. This is the C++
// counterpart of mizer's internal `bin_average_weight()`, which the package's
// `get...()` functions use, so that the quantities the objective function
// penalises are the ones those functions report. The top bin has no right-hand
// neighbour and is left unaveraged, as in mizer.
template<class Type>
vector<Type> bin_average_weight(vector<Type> K)
{
    int n = K.size();
    if (n < 2) return K;
    vector<Type> Kbar(n);
    for (int i = 0; i < n - 1; ++i) {
        Kbar(i) = Type(0.5) * (K(i) + K(i + 1));
    }
    Kbar(n - 1) = K(n - 1);
    return Kbar;
}

template<class Type>
vector<Type> calculate_F_mort(Type sel_func, Type logit_l50, Type log_ratio_left, Type log_l50_right_offset, Type log_ratio_right,
                              Type log_catchability, vector<Type> L, Type min_len, Type max_len,
                              vector<Type> w, vector<Type> dw, Type b_lw, int second_order)
{
    Type l50 = min_len + (max_len - min_len) * invlogit(logit_l50);
    Type l25 = l50 * (1 - invlogit(log_ratio_left));
    Type sr = l50 - l25;
    Type s1 = l50 * log(Type(3.0)) / sr;
    Type s2 = s1 / l50;

    Type catchability = exp(log_catchability);

    vector<Type> F_mort(L.size());

    // Setup for double sigmoid
    Type l50_right = Type(0.0);
    Type l25_right = Type(0.0);
    Type sr_right = Type(0.0);
    Type s1_right = Type(0.0);
    Type s2_right = Type(0.0);

    if(sel_func == 1) {
        l50_right = l50 + exp(log_l50_right_offset);
        l25_right = l50_right * (1 + exp(log_ratio_right));
        sr_right = l50_right - l25_right;
        s1_right = l50_right * log(Type(3.0)) / sr_right;
        s2_right = s1_right / l50_right;
    }

    if (second_order == 1) {
        // Bin-averaged selectivity
        int Q = 100;
        for (int i = 0; i < L.size(); i++) {
            Type w_val = w(i);
            Type dw_val = dw(i);
            Type beta = Type(1.0) + dw_val / w_val;
            
            Type sum_sel_w = Type(0.0);
            Type sum_w = Type(0.0);
            
            for (int q = 1; q <= Q; q++) {
                Type q_frac = (Type(q) - Type(0.5)) / Type(Q);
                Type w_q = w_val * pow(beta, q_frac);
                Type l_q = L(i) * pow(beta, q_frac / b_lw);
                
                Type sel = Type(0.0);
                if (sel_func == 1) {
                    sel = (Type(1.0) / (Type(1.0) + exp(s1 - s2 * l_q))) *
                          (Type(1.0) / (Type(1.0) + exp(s1_right - s2_right * l_q)));
                } else if (sel_func == 2) {
                    sel = (Type(1.0) / (Type(1.0) + exp(s1 - s2 * l_q)));
                }
                
                sum_sel_w += sel * w_q;
                sum_w += w_q;
            }
            
            F_mort(i) = catchability * (sum_sel_w / sum_w);
        }
    } else {
        // Point-sampled selectivity
        if(sel_func == 1) {
            for (int i = 0; i < L.size(); i++) {
                Type ilength = L(i);
                Type sel = (Type(1.0) / (Type(1.0) + exp(s1 - s2 * ilength))) *
                    (Type(1.0) / (Type(1.0) + exp(s1_right - s2_right * ilength)));
                F_mort(i) = catchability * sel;
            }

        } else if(sel_func == 2){
            for (int i = 0; i < L.size(); i++) {
                Type ilength = L(i);
                Type sel = (Type(1.0) / (Type(1.0) + exp(s1 - s2 * ilength)));
                F_mort(i) = catchability * sel;
            }
        }
    }

    return F_mort;
}

template<class Type>
vector<Type> calculate_growth(vector<Type> w, vector<Type> ergr,
                              vector<Type> matur, Type m, Type n, Type w_repro_max,
                              int defaults_edition)
{
    // Matches mizer's psi formula: repro_prop = pmin(1, (w/w_repro_max)^(m-n))
    // then psi = maturity * repro_prop, also capped at 1
    vector<Type> repro_prop = pow(w / w_repro_max, m - n);
    for (int i = 0; i < repro_prop.size(); i++) {
        if (repro_prop(i) > Type(1.0)) repro_prop(i) = Type(1.0);
    }
    vector<Type> psi = matur * repro_prop;
    for (int i = 0; i < psi.size(); i++) {
        if (psi(i) > Type(1.0)) psi(i) = Type(1.0);
    }
    
    // Set psi for all w > w_repro_max to 1 if defaults_edition < 2
    if (defaults_edition < 2) {
        for (int i = 0; i < w.size(); i++) {
            if (w(i) > w_repro_max) {
                psi(i) = Type(1.0);
            }
        }
    }
    
    vector<Type> growth = (Type(1.0) - psi) * ergr;
    for (int i = 0; i < growth.size(); i++) {
        if (growth(i) < Type(0.0)) growth(i) = Type(0.0);
    }

    return growth;
}

template<class Type>
vector<Type> calculate_N(vector<Type> mort, vector<Type> growth,
                         vector<Type> dw, vector<Type> d)
{
    // Solve the diffusion-aware single-species steady-state spectrum, matching
    // mizer's default upwind scheme (`get_transport_coefs_upwind()` evaluated at
    // dt = 1 with the fixed egg-density boundary of `get_steady_state_n()`, see
    // mizer's numerical_details vignette, "Steady-State Solution"). `growth` is
    // the growth velocity at the bin boundaries g_j, `mort` the bin-averaged
    // mortality mu_j, `dw` the bin widths and `d` the diffusion rate d_j.
    // With d == 0 the upper diagonal vanishes and this reduces exactly to the
    // pure-advection recursion N(j) = N(j-1) g(j-1) / (g(j) + mu(j) dw(j)).
    int size = dw.size();
    vector<Type> a(size), b(size), c(size), S(size);

    // Egg boundary: N(0) is fixed to 1; the overall scale is normalised later.
    a(0) = Type(0.0);
    b(0) = Type(1.0);
    c(0) = Type(0.0);
    S(0) = Type(1.0);

    for (int j = 1; j < size; ++j) {
        Type inv_dw = Type(1.0) / dw(j);
        a(j) = -inv_dw * (growth(j - 1) + Type(0.5) * d(j - 1) / dw(j - 1));
        b(j) = mort(j) + inv_dw * (growth(j)
                                   + Type(0.5) * d(j) / dw(j)
                                   + Type(0.5) * d(j) / dw(j - 1));
        c(j) = (j < size - 1) ?
            (-inv_dw * Type(0.5) * d(j + 1) / dw(j)) : Type(0.0);
        S(j) = Type(0.0);
    }

    // Thomas algorithm (double sweep) for the tridiagonal system.
    vector<Type> cp(size), sp(size), N(size);
    cp(0) = c(0) / b(0);
    sp(0) = S(0) / b(0);
    for (int j = 1; j < size; ++j) {
        Type denom = b(j) - a(j) * cp(j - 1);
        cp(j) = c(j) / denom;
        sp(j) = (S(j) - a(j) * sp(j - 1)) / denom;
    }
    N(size - 1) = sp(size - 1);
    for (int j = size - 2; j >= 0; --j) {
        N(j) = sp(j) - cp(j) * N(j + 1);
    }

    return N;
}

template<class Type>
vector<Type> calculate_N_second_order(vector<Type> mort, vector<Type> growth,
                                      vector<Type> dw, vector<Type> w, vector<Type> d, vector<Type> chi)
{
    int size = dw.size();
    vector<Type> a(size), b(size), c(size), S(size);

    // Spacing on the log-size grid
    Type h = log(w(1) / w(0));
    Type beta = w(1) / w(0);

    // Egg boundary: N(0) is fixed to 1; the overall scale is normalised later.
    a(0) = Type(0.0);
    b(0) = Type(1.0);
    c(0) = Type(0.0);
    S(0) = Type(1.0);

    for (int j = 1; j < size; ++j) {
        Type inv_dw = Type(1.0) / dw(j);
        
        Type gj = growth(j);
        Type chij = chi(j);
        
        Type gj_plus_1 = (j < size - 1) ? growth(j + 1) : growth(j);
        Type chij_plus_1 = (j < size - 1) ? chi(j + 1) : Type(0.0);
        
        Type wj = w(j);
        Type wj_plus_1 = (j < size - 1) ? w(j + 1) : w(j) * beta;
        
        Type dj = d(j);
        Type dj_minus_1 = d(j - 1);
        Type dj_plus_1 = (j < size - 1) ? d(j + 1) : d(j);

        a(j) = -inv_dw * (gj * (Type(1.0) - Type(0.5) * chij) + dj_minus_1 / (Type(2.0) * h * wj));
        
        b(j) = mort(j) + inv_dw * (gj_plus_1 * (Type(1.0) - Type(0.5) * chij_plus_1) 
                                   - Type(0.5) * chij * gj 
                                   + dj * (Type(1.0) / (Type(2.0) * h * wj) + Type(1.0) / (Type(2.0) * h * wj_plus_1)));
        
        c(j) = (j < size - 1) ?
            (inv_dw * (Type(0.5) * chij_plus_1 * gj_plus_1 - dj_plus_1 / (Type(2.0) * h * wj_plus_1))) : Type(0.0);
            
        S(j) = Type(0.0);
    }

    // Thomas algorithm (double sweep) for the tridiagonal system.
    vector<Type> cp(size), sp(size), N(size);
    cp(0) = c(0) / b(0);
    sp(0) = S(0) / b(0);
    for (int j = 1; j < size; ++j) {
        Type denom = b(j) - a(j) * cp(j - 1);
        cp(j) = c(j) / denom;
        sp(j) = (S(j) - a(j) * sp(j - 1)) / denom;
    }
    N(size - 1) = sp(size - 1);
    for (int j = size - 2; j >= 0; --j) {
        N(j) = sp(j) - cp(j) * N(j + 1);
    }

    return N;
}

template<class Type>
Type objective_function<Type>::operator() ()
{

    DATA_INTEGER(use_counts);
    DATA_MATRIX(counts);
    DATA_IVECTOR(bin_index);
    DATA_IVECTOR(f_index);
    DATA_VECTOR(sel_func);
    DATA_VECTOR(coeff_fj);
    DATA_VECTOR(coeff_fj1);

    DATA_VECTOR(dw);
    DATA_VECTOR(w);
    DATA_VECTOR(l);
    DATA_VECTOR(yield);
    DATA_SCALAR(production);
    DATA_SCALAR(biomass);
    DATA_INTEGER(biomass_cutoff_idx);
    // DATA_VECTOR(growth);
    DATA_SCALAR(w_mat);
    DATA_SCALAR(d);
    DATA_VECTOR(yield_weight);
    DATA_VECTOR(catch_dist_weight);
    DATA_SCALAR(production_weight);
    DATA_VECTOR(matur);
    DATA_VECTOR(ergr);
    DATA_SCALAR(n);
    DATA_SCALAR(w_repro_max);
    DATA_INTEGER(second_order);
    // Whether the summary integrals (biomass, yield, production) trapezoidally
    // bin-average their weight, as mizer's own summary functions and the
    // package's `get...()` functions do when
    // `second_order_w(params)[["bin_average"]]` is TRUE. This is a separate
    // switch from `second_order`, which selects the flux scheme of the steady
    // state solver.
    DATA_INTEGER(bin_average);
    DATA_VECTOR(chi);
    DATA_INTEGER(defaults_edition);
    DATA_SCALAR(b_lw);

    PARAMETER_VECTOR(logit_l50);
    PARAMETER_VECTOR(log_ratio_left);
    PARAMETER_VECTOR(log_l50_right_offset);
    PARAMETER_VECTOR(log_ratio_right);
    PARAMETER_VECTOR(log_catchability);
    PARAMETER(z_ext);
    PARAMETER(m);
    PARAMETER(log_D_ext);

    int n_bins = w.size();
    int n_g = yield.size();

    Type min_len = l(0);
    Type max_len = l(n_bins - 1);

    vector<Type> total_F_mort(n_bins);
    total_F_mort.fill(Type(0.0));

    matrix<Type> F_mort_mat(n_bins, n_g);
    F_mort_mat.fill(Type(0.0));

    for (int g = 0; g < n_g; ++g) {
        vector<Type> F_mort_g = calculate_F_mort(sel_func[g], logit_l50[g],
                                                 log_ratio_left[g], log_l50_right_offset[g],
                                                 log_ratio_right[g], log_catchability[g], l, min_len, max_len,
                                                 w, dw, b_lw, second_order);
        for (int i = 0; i < n_bins; ++i) {
            F_mort_mat(i, g) = F_mort_g(i);
            total_F_mort(i) += F_mort_g(i);
        }
    }

    vector<Type> mort = z_ext * pow(w, d) + total_F_mort;

    vector<Type> growth = calculate_growth(w, ergr, matur, m, n, w_repro_max, defaults_edition);

    // External diffusion rate as a power law d(w) = D_ext * w^(n+1), matching
    // the model's `ext_diffusion` slot (getDiffusion()) when use_predation_diffusion
    // is FALSE.
    Type D_ext = exp(log_D_ext);
    vector<Type> d_diff(n_bins);
    if (second_order == 1) {
        Type d_exp = n + Type(1.0);
        for (int i = 0; i < n_bins; ++i) {
            Type w_next = w(i) + dw(i);
            d_diff(i) = D_ext * (pow(w_next, d_exp + Type(1.0)) - pow(w(i), d_exp + Type(1.0))) / ((d_exp + Type(1.0)) * dw(i));
        }
    } else {
        for (int i = 0; i < n_bins; ++i) {
            d_diff(i) = D_ext * pow(w(i), n + Type(1.0));
        }
    }

    vector<Type> N;
    if (second_order == 1) {
        N = calculate_N_second_order(mort, growth, dw, w, d_diff, chi);
    } else {
        N = calculate_N(mort, growth, dw, d_diff);
    }

    // Composite weight of the biomass integral: the size window given by
    // `biomass_cutoff_idx` times w. The window is part of the weight so that,
    // when bin-averaging, the bin straddling the cutoff contributes only
    // partially, which is what mizer::getBiomass() does.
    vector<Type> w_biomass(n_bins);
    for (int i = 0; i < n_bins; ++i) {
        w_biomass(i) = (i >= biomass_cutoff_idx) ? w(i) : Type(0.0);
    }
    if (bin_average == 1) {
        w_biomass = bin_average_weight(w_biomass);
    }

    Type unscaled_biomass = 0;
    for (int i = 0; i < n_bins; ++i) {
        unscaled_biomass += N(i) * w_biomass(i) * dw(i);
    }

    N *= (biomass / unscaled_biomass);

    vector<Type> catch_dens = N * total_F_mort;

    matrix<Type> catch_per_bin_mat(n_bins,n_g);
    catch_per_bin_mat.fill(Type(0.0));

    matrix<Type> dens_per_bin_mat(n_bins,n_g);
    dens_per_bin_mat.fill(Type(0.0));

    vector<Type> model_yield(n_g);
    model_yield.fill(Type(0.0));

    Type nll = Type(0.0);
    Type multinomial_nll = Type(0.0);
    Type yield_nll = Type(0.0);

    int c_bins = counts.rows();
    int num_segs = bin_index.size();

    if (use_counts == 1) {

        for (int g = 0; g < n_g; ++g) {

            vector<Type> F_mort_g(n_bins);
            for (int i = 0; i < n_bins; i++) {
                F_mort_g(i) = F_mort_mat(i, g);
            }

            vector<Type> catch_per_bin_g = N * F_mort_g * dw;
            vector<Type> dens_per_bin_g = N * F_mort_g;
            catch_per_bin_mat.col(g) = catch_per_bin_g;
            dens_per_bin_mat.col(g) = dens_per_bin_g;
            
            // Composite weight of the yield integral: F_mort_g * w, as in
            // mizer::getYield().
            vector<Type> w_yield = F_mort_g * w;
            if (bin_average == 1) {
                w_yield = bin_average_weight(w_yield);
            }
            Type model_yield_val = Type(0.0);
            for (int i = 0; i < n_bins; ++i) {
                model_yield_val += N(i) * w_yield(i) * dw(i);
            }
            model_yield(g) = model_yield_val;

            vector<Type> counts_g(c_bins);
            counts_g.fill(Type(0.0));
            for (int i = 0; i < c_bins; i++) {
                counts_g(i) = counts(i, g);
            }

            vector<Type> probs(c_bins);
            probs.fill(Type(0.0));

            for (int k = 0; k < num_segs; k++) {
                int i = bin_index(k);
                int j = f_index(k);

                probs(i) += coeff_fj(k) * dens_per_bin_g(j) + coeff_fj1(k) * dens_per_bin_g(j+1);
            }

            // Normalise before applying the floor that guards against log(0), so
            // that this term is invariant under a rescaling of the gear's
            // catchability. An absolute floor swamps the signal for a gear with
            // very small catchability (e.g. a survey gear at 1e-12), leaving the
            // likelihood flat in l50/l25 so that the optimiser never moves away
            // from its starting values.
            probs /= probs.sum();
            probs += Type(1e-10);
            probs /= probs.sum();

            // The species-gear-specific `catch_dist_weight` weights this gear's
            // penalty for the mismatch between the observed and the modelled
            // catch size distribution.
            Type mult_val = catch_dist_weight(g) *
                (-dmultinom(counts_g, probs, true) / counts_g.sum());
            multinomial_nll += mult_val;
            nll += mult_val;

            // The species-gear-specific `yield_weight` weights this gear's
            // yield penalty.
            if (yield_weight(g) > 0) {
                Type yield_val = yield_weight(g) * pow(log(model_yield(g) / yield(g)), Type(2));
                yield_nll += yield_val;
                nll += yield_val;
            }

        }

    }

    // else {
    //
    //   if (yield_weight(0) > 0) {
    //     // **Calculate model yield**
    //     vector<Type> yield_per_bin = catch_dens * w * dw;
    //     Type model_yield = yield_per_bin.sum();
    //     REPORT(model_yield);
    //     // Penalise deviation from observed yield
    //     nll += yield_weight(0) * pow(log(model_yield / yield(0)), Type(2));
    //   }
    //
    // }

    Type production_nll = Type(0.0);
    if (production_weight > 0) {
        // **Calculate production**
        // The weight of the integral is the growth rate, as in
        // getSomaticProduction().
        vector<Type> w_production = growth;
        if (bin_average == 1) {
            w_production = bin_average_weight(w_production);
        }
        vector<Type> production_per_bin = N * w_production * dw;
        Type model_production = production_per_bin.sum();
        REPORT(model_production);
        // The species-specific `production_weight` weights this penalty for
        // deviation from the observed production.
        production_nll = production_weight * pow(log(model_production / production), Type(2));
        nll += production_nll;
    }

    // Final sanity checks
    TMBAD_ASSERT(nll >= 0);
    TMBAD_ASSERT(CppAD::isfinite(nll));
    if (!CppAD::isfinite(nll)) {
        error("nll is not finite");
    }

    // **Reporting**

    REPORT(N);
    REPORT(total_F_mort);
    REPORT(mort);
    REPORT(growth);
    REPORT(d_diff);
    REPORT(model_yield);
    REPORT(multinomial_nll);
    REPORT(yield_nll);
    REPORT(production_nll);

    // Check final biomass again
    Type total_biomass = 0;
    for (int i = 0; i < n_bins; ++i) {
        total_biomass += N(i) * w_biomass(i) * dw(i);
    }
    REPORT(total_biomass);

    return nll;
}
