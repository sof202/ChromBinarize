#!/usr/bin/awk -f
function min(x,y) {
    return x>y ? y : x
}
function max(x,y) {
    return x>y ? x : y
}
function calculate_agresti_coull(x, n, z_a, is_wgbs) {
    # x -> successes
    # n -> trials
    # za -> quantile of a standard normal distribution
    n_hat = n + (z_a * z_a);
    p_hat = (1 / n_hat) * (x + (z_a * z_a) / 2);
    confidence_interval_radius = z_a * sqrt(p_hat * (1 - p_hat) / n_hat);

    # These coonfidence intervals are calculated to check overlaps. As such
    # we only need to calculate the lower bound for WGBS and upper bound for
    # oxBS to see if WGBS signal is significantly higher than oxBS signal.
    if (is_wgbs) {
        return max(p_hat - confidence_interval_radius, 0);
    }
    return min(p_hat + confidence_interval_radius, 1);
}
BEGIN{OFS="\t"}
{print $1,$2,$3,$4,$5,calculate_agresti_coull($4,$5,za,is_wgbs)}
