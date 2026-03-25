#include <TMB.hpp>

template<class Type>
vector<Type> calculate_F_mort(Type sel_func, Type logit_l50, Type log_ratio_left, Type log_l50_right_offset, Type log_ratio_right,
                              Type log_catchability, vector<Type> L, Type min_len, Type max_len)
{
    Type l50 = min_len + (max_len - min_len) * invlogit(logit_l50);
    Type l25 = l50 * (1 - invlogit(log_ratio_left));
    Type sr = l50 - l25;
    Type s1 = l50 * log(Type(3.0)) / sr;
    Type s2 = s1 / l50;

    Type catchability = exp(log_catchability);

    vector<Type> F_mort(L.size());

    if(sel_func == 1) {
        Type l50_right = l50 + exp(log_l50_right_offset);
        Type l25_right = l50_right * (1 + exp(log_ratio_right));
        Type sr_right = l50_right - l25_right;
        Type s1_right = l50_right * log(Type(3.0)) / sr_right;
        Type s2_right = s1_right / l50_right;

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

    return F_mort;
}

template<class Type>
vector<Type> calculate_growth(vector<Type> w, vector<Type> ergr,
                              vector<Type> matur, Type m, Type n, Type w_max)
{
    vector<Type> phi = matur * pow( w/w_max, m-n);
    vector<Type> growth = (1- phi) * ergr;

    return growth;
}

template<class Type>
vector<Type> calculate_N(vector<Type> mort, vector<Type> growth,
                         vector<Type> dw)
{
    int size = dw.size();
    vector<Type> N(size);
    N(0) = Type(1.0);
    for (int i = 1; i < size; ++i) {
        Type denominator = growth(i) + mort(i) * dw(i);
        N(i) = N(i - 1) * growth(i - 1) / denominator;
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
    DATA_SCALAR(yield_lambda);
    DATA_SCALAR(production_lambda);
    DATA_VECTOR(matur);
    DATA_VECTOR(ergr);
    DATA_SCALAR(n);
    DATA_SCALAR(w_max);

    PARAMETER_VECTOR(logit_l50);
    PARAMETER_VECTOR(log_ratio_left);
    PARAMETER_VECTOR(log_l50_right_offset);
    PARAMETER_VECTOR(log_ratio_right);
    PARAMETER_VECTOR(log_catchability);
    PARAMETER(mu_mat);
    PARAMETER(m);

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
                                                                                        log_ratio_right[g], log_catchability[g], l, min_len, max_len);
        for (int i = 0; i < n_bins; ++i) {
            F_mort_mat(i, g) = F_mort_g(i);
            total_F_mort(i) += F_mort_g(i);
        }
    }

    vector<Type> mort = mu_mat * pow(w / w_mat, d) + total_F_mort;

    vector<Type> growth = calculate_growth(w, ergr, matur, m, n, w_max);

    vector<Type> N = calculate_N(mort, growth, dw);

    Type unscaled_biomass = 0;
    for (int i = biomass_cutoff_idx; i < n_bins; ++i) {
        unscaled_biomass += N(i) * w(i) * dw(i);
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
            vector<Type> model_yield_g = catch_per_bin_g * w;
            model_yield(g) = model_yield_g.sum();

            vector<Type> counts_g(c_bins);
            counts_g.fill(Type(0.0));
            for (int i = 0; i < c_bins; i++) {
                counts_g(i) = counts(i, g);
            }

            vector<Type> probs(c_bins);
            probs.fill(Type(1e-10));

            for (int k = 0; k < num_segs; k++) {
                int i = bin_index(k);
                int j = f_index(k);

                probs(i) += coeff_fj(k) * dens_per_bin_g(j) + coeff_fj1(k) * dens_per_bin_g(j+1);
            }

            probs /= probs.sum();

            nll -= dmultinom(counts_g, probs, true) / counts_g.sum();

            if (yield_lambda > 0) {
                nll += yield_lambda * pow(log(model_yield(g) / yield(g)), Type(2));
            }

        }

    }

    // else {
    //
    //   if (yield_lambda > 0) {
    //     // **Calculate model yield**
    //     vector<Type> yield_per_bin = catch_dens * w * dw;
    //     Type model_yield = yield_per_bin.sum();
    //     REPORT(model_yield);
    //     // Penalise deviation from observed yield
    //     nll += yield_lambda * pow(log(model_yield / yield), Type(2));
    //   }
    //
    // }

    if (production_lambda > 0) {
        // **Calculate production**
        vector<Type> production_per_bin = N * growth * dw;
        Type model_production = production_per_bin.sum();
        REPORT(model_production);
        // **Add penalty for deviation from observed production**
        nll += production_lambda * pow(log(model_production / production), Type(2));
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

    // Check final biomass again
    Type total_biomass = 0;
    for (int i = biomass_cutoff_idx; i < N.size(); ++i) {
        total_biomass += N(i) * w(i) * dw(i);
    }
    REPORT(total_biomass);

    return nll;
}
