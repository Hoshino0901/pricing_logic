using MathNet.Numerics.Providers.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;

namespace PricingLogic
{
    public class AnalyticBS
    {
        public double S0 { get; private set; }
        public double InteRate { get; private set; }
        public double Vola { get; private set; }
        public double Maturity { get; private set; }
        public double Strike { get; private set; }

        public AnalyticBS(double s0, double r, double sigma, double T, double K)
        {
            S0 = s0;
            InteRate = r;
            Vola = sigma;
            Maturity = T;
            Strike = K;
        }

        public double AnalyticBSFormula()
        {
            double d_plus = 1 / (Vola * Math.Sqrt(Maturity)) * (Math.Log(S0 / Strike) + (InteRate + 0.5 * Vola * Vola) * Maturity);
            double d_minus = d_plus - Vola * Math.Sqrt(Maturity);
            var norm = new Normal(0, 1);
            return S0 * norm.CumulativeDistribution(d_plus) - Strike * Math.Exp(-InteRate * Maturity) * norm.CumulativeDistribution(d_minus);
        } 
    }
}
