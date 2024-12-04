#include <TMB.hpp>

template<class Type>
vector<Type> calculate_F_mort(Type l50, Type ratio, Type catchability,
                              vector<Type> bin_boundary_lengths)
{
    Type c1 = Type(1.0);
    Type sr = l50 * (c1 - ratio);
    Type s1 = l50 * log(Type(3.0)) / sr;
    Type s2 = s1 / l50;
    vector<Type> F_mort = catchability /
        (c1 + exp(s1 - s2 * bin_boundary_lengths));

    // Check that all elements are finite and non-negative
    TMBAD_ASSERT((F_mort.array().isFinite() && (F_mort.array() >= 0)).all());

    return F_mort;
}

template<class Type>
vector<Type> calculate_mort(vector<Type> F_mort, Type M, Type d,
                            vector<Type> bin_boundaries)
{
    vector<Type> mort = M * pow(bin_boundaries, d) + F_mort;

    // Check that all elements are finite and positive
    TMBAD_ASSERT((mort.array().isFinite() && (mort.array() > 0)).all());
    return mort;
}

template<class Type>
vector<Type> calculate_growth(vector<Type>EReproAndGrowth,
                              vector<Type>repro_prop,
                              Type w_mat, Type U,
                              vector<Type> bin_boundaries)
{
    Type c1 = Type(1.0);
    vector<Type> psi = repro_prop / (c1 + pow(bin_boundaries / w_mat, -U));
    vector<Type> growth = EReproAndGrowth * (c1 - psi);

    // Check that all elements are finite and non-negative
    TMBAD_ASSERT((growth.array().isFinite() && (growth.array() >= 0)).all());
    return growth;
}

template<class Type>
vector<Type> calculate_N(vector<Type> mort, vector<Type> growth,
                         Type biomass,
                         vector<Type> bin_widths,
                         vector<Type> bin_boundaries)
{
    int size = bin_widths.size();
    vector<Type> N(size + 1);
    N(0) = Type(1.0);
    for (int i = 1; i < size; ++i) {
        Type denominator = growth(i) + mort(i) * bin_widths(i);
        N(i) = N(i - 1) * growth(i - 1) / denominator;
    }
    N(size) = N(size - 1) * growth(size - 1) /
        (mort(size - 1) * bin_widths(size - 1) + growth(size - 1));

    // Rescale to get observed biomass
    Type total_biomass = Type(0.0);
    for (int i = 0; i < size; ++i) {
        total_biomass += N(i) * bin_widths(i) *
            (bin_boundaries(i) + bin_boundaries(i + 1)) / Type(2.0);
    }
    N = N * biomass / total_biomass;

    // Check that all elements are finite and non-negative
    TMBAD_ASSERT((N.array().isFinite() && (N.array() >= 0)).all());

    return N;
}

template<class Type>
vector<Type> calculate_catch_per_bin(vector<Type> N, vector<Type> F_mort,
                                     vector<Type> bin_widths)
{
    // **Calculate catch Densities at Bin Boundaries**
    vector<Type> densities = N * F_mort;
    TMBAD_ASSERT((densities.array() >= 0).all());


    // **Integrate density over bin using Trapezoidal Rule**
    int num_bins = bin_widths.size(); // Number of bins
    TMBAD_ASSERT(densities.size() == num_bins + 1);
    vector<Type> catch_per_bin(num_bins);
    for (int i = 0; i < num_bins; ++i) {
        // Trapezoidal rule
        catch_per_bin[i] = bin_widths[i] *
            (densities[i] + densities[i + 1]) / Type(2.0);
    }
    return catch_per_bin;
}

template<class Type>
Type calculate_yield(vector<Type> catch_per_bin,
                            vector<Type> bin_boundaries)
{
    // **Calculate model yield**
    Type model_yield = Type(0.0);
    for (int i = 0; i < catch_per_bin.size(); ++i) {
        model_yield += catch_per_bin[i] *
            (bin_boundaries[i] + bin_boundaries[i + 1]) / Type(2.0);
    }
    return model_yield;
}



template<class Type>
Type objective_function<Type>::operator() ()
{
    // **Data Section**
    DATA_VECTOR(counts);           // Count observations in bins
    DATA_VECTOR(bin_widths);       // Width of each bin in grams
    DATA_VECTOR(bin_boundaries);   // Boundaries of each bin in grams
    DATA_VECTOR(bin_boundary_lengths);  // Boundaries of each bin in cm
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
    PARAMETER(M);            // Coefficient of natural mortality rate power law
    PARAMETER(U);            // Steepness parameter of maturity ogive
    PARAMETER(catchability); // Catchability

    // Check lengths of data vectors
    TMBAD_ASSERT(bin_widths.size() == bin_boundaries.size() - 1);
    TMBAD_ASSERT(bin_boundaries.size() == bin_boundary_lengths.size());
    TMBAD_ASSERT(counts.size() == bin_widths.size());
    TMBAD_ASSERT(EReproAndGrowth.size() == bin_boundaries.size());
    TMBAD_ASSERT(repro_prop.size() == bin_boundaries.size());

    // **Calculate fishing mortality rate**
    vector<Type> F_mort = calculate_F_mort(l50, ratio, catchability,
                                           bin_boundary_lengths);

    // **Calculate total mortality rate**
    vector<Type> mort = calculate_mort(F_mort, M, d, bin_boundaries);

    // **Calculate growth rate**
    vector<Type> growth = calculate_growth(EReproAndGrowth, repro_prop,
                                           w_mat, U, bin_boundaries);

    // **Calculate steady-state number density**
    vector<Type> N = calculate_N(mort, growth, biomass,
                                 bin_widths, bin_boundaries);

    // **Calculate catch per bin**
    vector<Type> catch_per_bin = calculate_catch_per_bin(N, F_mort, bin_widths);

    // **Calculate model yield**
    Type model_yield = calculate_yield(catch_per_bin, bin_boundaries);

    // **Calculate catch probabilities**
    // Ensure Probabilities Are Positive and Sum to 1
    // Add a small epsilon to avoid log(0) and normalize
    vector<Type> probs = catch_per_bin + Type(1e-10);
    probs = probs / probs.sum();
    REPORT(probs);

    // **Negative Log-Likelihood Calculation**
    // Compute the negative log-likelihood using the multinomial distribution
    Type nll = -dmultinom(counts, probs, true) / counts.sum();

    // **Add penalty for deviation from observed yield**
    nll += yield_lambda * pow(log(model_yield / yield), Type(2));

    TMBAD_ASSERT(nll >= 0);
    TMBAD_ASSERT(CppAD::isfinite(nll));
    if (!CppAD::isfinite(nll)) error("nll is not finite");

    REPORT(model_yield);

    return nll;
}
