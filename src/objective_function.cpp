#include <TMB.hpp>

template<class Type>
vector<Type> calculate_F_mort(Type l50, Type ratio, Type catchability,
                              vector<Type> l)
{
    Type c1 = Type(1.0);
    Type sr = l50 * (c1 - ratio);
    Type s1 = l50 * log(Type(3.0)) / sr;
    Type s2 = s1 / l50;
    vector<Type> F_mort = catchability /
        (c1 + exp(s1 - s2 * l));

    // Check that all elements are finite and non-negative
    TMBAD_ASSERT((F_mort.array().isFinite() && (F_mort.array() >= 0)).all());

    return F_mort;
}

template<class Type>
vector<Type> calculate_growth(vector<Type>EReproAndGrowth,
                              vector<Type>repro_prop,
                              Type w_mat, Type U,
                              vector<Type> w)
{
    Type c1 = Type(1.0);
    vector<Type> psi = repro_prop / (c1 + pow(w / w_mat, -U));
    vector<Type> growth = EReproAndGrowth * (c1 - psi);

    // Check that all elements are finite and non-negative
    TMBAD_ASSERT((growth.array().isFinite() && (growth.array() >= 0)).all());
    return growth;
}

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

    // Check that all elements are finite and non-negative
    TMBAD_ASSERT((N.array().isFinite() && (N.array() >= 0)).all());

    return N;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
    // **Data Section**
    DATA_VECTOR(counts);           // Count observations in bins
    DATA_IVECTOR(bin_index);     // Bin indices for overlapping segments
    DATA_IVECTOR(f_index);       // Function indices (j) for overlapping segments
    DATA_VECTOR(coeff_fj);       // Coefficients for f(j)
    DATA_VECTOR(coeff_fj1);      // Coefficients for f(j+1)
    DATA_VECTOR(dw);
    DATA_VECTOR(w);
    DATA_VECTOR(l);                // lengths corresponding to w
    DATA_SCALAR(yield);            // Observed yield
    DATA_SCALAR(biomass);          // Observed biomass
    DATA_VECTOR(EReproAndGrowth);  // The rate at which energy is available for growth
                                   // and reproduction
    DATA_VECTOR(repro_prop);       // Proportion of energy allocated to reproduction
    DATA_SCALAR(w_mat);            // Maturity size is currently not optimised
    DATA_SCALAR(d);                // Exponent of mortality power-law
    DATA_SCALAR(yield_lambda);     // controls the strength of the penalty for
                                   // deviation from the observed yield.

    // **Parameter Section**
    PARAMETER(l50);          // Length at 50% gear selectivity
    PARAMETER(ratio);        // Ratio between l25 and l50
    PARAMETER(mu_mat);       // Mortality at maturity size
    PARAMETER(U);            // Steepness parameter of maturity ogive
    PARAMETER(catchability); // Catchability

    // **Calculate fishing mortality rate**
    vector<Type> F_mort = calculate_F_mort(l50, ratio, catchability, l);

    // **Calculate total mortality rate**
    vector<Type> mort = mu_mat * pow(w / w_mat, d) + F_mort;

    // **Calculate growth rate**
    vector<Type> growth = calculate_growth(EReproAndGrowth, repro_prop,
                                           w_mat, U, w);

    // **Calculate steady-state number density**
    vector<Type> N = calculate_N(mort, growth, biomass, w, dw);

    // Rescale to get observed biomass
    vector<Type> biomass_in_bins = N * w * dw;
    N = N * biomass / biomass_in_bins.sum();

    // **Calculate catch density**
    vector<Type> catch_dens = N * F_mort;

    // **Calculate model yield**
    vector<Type> yield_per_bin = catch_dens * w * dw;
    Type model_yield = yield_per_bin.sum();

    // **Calculate catch probabilities**
    int num_bins = counts.size();    // Number of bins
    int num_segs = bin_index.size(); // Number of overlapping segments

    vector<Type> probs(num_bins);   // Vector to store bin probabilities
    probs.setZero();             // Initialize to zero

    // Compute bin probabilities using precomputed weights
    for (int k = 0; k < num_segs; k++) { // Loop over segments
        int i = bin_index(k);      // Bin index
        int j = f_index(k);        // Function index

        // Accumulate the contributions to the bin probability
        probs(i) += coeff_fj(k) * catch_dens(j) + coeff_fj1(k) * catch_dens(j+1);
    }

    // Avoid 0 probability for non-zero counts
    for (int i = 0; i < probs.size(); ++i) {
        if (probs[i] == 0 && counts[i] == 0) probs[i] = 1e-10; // Small epsilon value
    }

    // Normalize the bin probabilities so they sum to 1
    probs = probs / probs.sum();

    // **Negative Log-Likelihood Calculation**
    // Compute the negative log-likelihood using the multinomial distribution
    Type nll = -dmultinom(counts, probs, true) / counts.sum();

    // **Add penalty for deviation from observed yield**
    nll += yield_lambda * pow(log(model_yield / yield), Type(2));

    TMBAD_ASSERT(nll >= 0);
    TMBAD_ASSERT(CppAD::isfinite(nll));
    if (!CppAD::isfinite(nll)) error("nll is not finite");

    REPORT(probs);
    REPORT(model_yield);
    REPORT(N);
    REPORT(F_mort);
    REPORT(mort);
    REPORT(growth);

    // Check that rescaling worked
    biomass_in_bins = N * w * dw;
    Type total_biomass = biomass_in_bins.sum();
    REPORT(total_biomass);

    return nll;
}
