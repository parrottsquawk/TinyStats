package common.algorithm;

import org.apache.commons.math3.stat.descriptive.StatisticalSummary;

public class TinyStats implements StatisticalSummary {
    private double n, nᝍ1, nᝍ2, nᝍ3; // number of items, decremented versions to save processing time.
    private double Σ; // Sum
    private double μ_1, μ_2, μ_3, μ_4; //first, second, third, and fourth moments.
    private double Σxlogx; // sum of x* log_2(x)
    private double max = Double.MIN_VALUE, min = Double.MAX_VALUE;

    public TinyStats(double[] data) {
        for (double d : data) {
            this.put(d);
        }
    }

    public TinyStats() {
    }

    public void put(double x) {
        nᝍ3 = nᝍ2;
        nᝍ2 = nᝍ1;
        nᝍ1 = n;
        n++;

        max = x > max ? x : max;
        min = x < min ? x : min;

        Σ += x;
        Σxlogx += x * Math.log(x) / Math.log(2);

        double δ = x - μ_1;
        double δ〳n = δ / n;
        double δᶺ2〳n = δ * δ〳n;
        double δᶺ3〳nᶺ2 = δ〳n * δᶺ2〳n;
        double t1 = δ〳n * 3 * μ_2;

        μ_4 += (δ〳n * δᶺ3〳nᶺ2) * nᝍ1 * (n * nᝍ3 + 3) + δ〳n * t1 * 2 - δ〳n * 4 * μ_3;
        μ_3 += δᶺ3〳nᶺ2 * nᝍ2 * nᝍ1 - t1;
        μ_2 += δᶺ2〳n * nᝍ1;
        μ_1 += δ〳n;

    }

    public double getMean() {
        return μ_1;
    }

    public double getSum() {
        return Σ;
    }

    public double getVariance() {
        return n < 2 ? Double.NaN : μ_2 / nᝍ1;
    }

    public double getStandardDeviation() {
        return Math.sqrt(getVariance());
    }

    @Override
    public double getMax() {
        return max;
    }

    @Override
    public double getMin() {
        return min;
    }

    @Override
    public long getN() {
        return (long) n;
    }

    public double getKurtosis() {
        return ((n * μ_4) / (μ_2 * μ_2));
    }

    public double getSkewness() {
        return (Math.sqrt(n) * μ_3) / Math.sqrt(μ_2 * μ_2 * μ_2);
    }

    public double getEntropy() {
        return (-Σxlogx * Math.log(n) / (n * Math.log(2)));
    }

    public double getExactHistogramBinSize() {
        return 3.49 * getStandardDeviation() * Math.pow(getN(), -(1.0 / 3.0));
    }
