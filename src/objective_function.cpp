#include <TMB.hpp>

// Helper function: Fishing mortality
template<class Type>
vector<Type> calculate_F_mort(Type l50, Type ratio, Type catchability,
                              vector<Type> l)
{
    Type c1 = Type(1.0);
    Type sr = l50 * (c1 - ratio);
    Type s1 = l50 * log(Type(3.0)) / sr;
    Type s2 = s1 / l50;

    vector<Type> F_mort = catchability / (c1 + exp(s1 - s2 * l));

    // Ensure all elements are finite and >= 0
    TMBAD_ASSERT((F_mort.array().isFinite() && (F_mort.array() >= 0)).all());

    return F_mort;
}

// Helper function: Steady-state numbers-at-size
template<class Type>
vector<Type> calculate_N(vector<Type> mort, vector<Type> growth,
                         Type biomass, vector<Type> w, vector<Type> dw)
{
    int size = dw.size();
    vector<Type> N(size);
    N(0) = Type(1.0);
    for (int i = 1; i < size; ++i) {
        Type denominator = growth(i) + mort(i) * dw(i);
        N(i) = N(i - 1) * growth(i - 1) / denominator;
    }

    // Ensure all elements are finite and >= 0
    TMBAD_ASSERT((N.array().isFinite() && (N.array() >= 0)).all());

    return N;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
    // **Data Section**
    // Introduce an integer flag to indicate whether 'counts' data are present.
    // Calling code can set 'use_counts = 1' if counts exist, or '0' if they do not.
    DATA_INTEGER(use_counts);
    // The next lines read counts data, but if use_counts = 0, we won't use them.
    DATA_VECTOR(counts);     // Observed count data for each bin
    DATA_IVECTOR(bin_index); // Bin indices (for overlapping segments)
    DATA_IVECTOR(f_index);   // Function indices (for overlapping segments)
    DATA_VECTOR(coeff_fj);
    DATA_VECTOR(coeff_fj1);

    DATA_VECTOR(dw);
    DATA_VECTOR(w);
    DATA_VECTOR(l);                 // lengths corresponding to w
    DATA_SCALAR(yield);             // Observed yield
    DATA_SCALAR(production);        // Observed production
    DATA_SCALAR(biomass);           // Observed biomass
    DATA_VECTOR(growth);            // Growth rate
    DATA_SCALAR(w_mat);             // Maturity size (weight)
    DATA_SCALAR(d);                 // Exponent for mortality power-law
    DATA_SCALAR(yield_lambda);      // Penalty strength for yield deviation
    DATA_SCALAR(production_lambda); // Penalty strength for production deviation

    // **Parameter Section**
    PARAMETER(l50);          // Length at 50% gear selectivity
    PARAMETER(ratio);        // Ratio between l25 and l50
    PARAMETER(mu_mat);       // Mortality at maturity size
    PARAMETER(catchability); // Catchability

    // **Calculate fishing mortality rate**
    vector<Type> F_mort = calculate_F_mort(l50, ratio, catchability, l);

    // **Calculate total mortality rate**
    vector<Type> mort = mu_mat * pow(w / w_mat, d) + F_mort;

    // **Calculate steady-state number density**
    //   This is unscaled (N is not matched to 'biomass' yet).
    vector<Type> N = calculate_N(mort, growth, biomass, w, dw);

    // Rescale to match observed biomass
    vector<Type> biomass_in_bins = N * w * dw;
    Type unscaled_biomass = biomass_in_bins.sum();
    N *= (biomass / unscaled_biomass);  // Now total biomass matches 'biomass'.

    // **Calculate catch density**
    vector<Type> catch_dens = N * F_mort;

    // **Negative Log-Likelihood (NLL)**
    // We initialize nll to zero and only add contributions that apply.
    Type nll = Type(0.0);

    // If 'use_counts == 1', we have counts data and proceed with the multinomial likelihood.
    if (use_counts == 1) {
        int num_bins = counts.size();
        int num_segs = bin_index.size();
        vector<Type> probs(num_bins);
        probs.fill(Type(1e-10)); // Avoid exactly zero probabilities

        // Accumulate the bin probabilities from overlapping segments
        for (int k = 0; k < num_segs; k++) {
            int i = bin_index(k);
            int j = f_index(k);
            probs(i) += coeff_fj(k) * catch_dens(j) + coeff_fj1(k) * catch_dens(j+1);
        }
        // Normalize probabilities
        probs /= probs.sum();

        // Multinomial negative log-likelihood
        // You might or might not want to divide by counts.sum()—that’s up to you.
        nll -= dmultinom(counts, probs, true) / counts.sum();
    }
    // else: skip all calculations involving 'counts'

    // **Add penalty for deviation from observed yield**
    if (yield_lambda > 0) {
        // **Calculate model yield**
        vector<Type> yield_per_bin = catch_dens * w * dw;
        Type model_yield = yield_per_bin.sum();
        REPORT(model_yield);
        // Penalise deviation from observed yield
        nll += yield_lambda * pow(log(model_yield / yield), Type(2));
    }

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
    REPORT(F_mort);
    REPORT(mort);

    // Check final biomass again
    biomass_in_bins = N * w * dw;
    Type total_biomass = biomass_in_bins.sum();
    REPORT(total_biomass);

    return nll;
}
