using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Optimization;
using MathNet.Numerics.Random;
using MathNet.Numerics.LinearRegression;

namespace PricingLogic
{
    public class LsmEuropean
    {
        public double S0 { get; private set; }
        public double InteRate { get; private set; }
        public double Vola { get; private set; }
        public double Maturity { get; private set; }
        public double Strike { get; private set; }

        public LsmEuropean(double s0, double r, double sigma, double T, double K)
        {
            S0 = s0;
            InteRate = r;
            Vola = sigma;
            Maturity = T;
            Strike = K;
        }
        private double EuropeanPayoff(double x, double t)
        {
            return Math.Max(x - Strike, 0) * Math.Exp(-InteRate * (Maturity - t));
        }
        private double[] GeneratePathByEM(double[] bm)
        {
            double dt = Maturity / bm.Length;
            double[] path = new double[bm.Length + 1];
            path[0] = S0;
            for (int i = 0; i < bm.Length; i++)
            {
                path[i + 1] = path[i] * (1 + InteRate * dt + Vola * Math.Sqrt(dt) * bm[i]);
            }
            return path;
        }

        public double LeastSquareMonte(int N, int M, int seed = 0)
        {
            var x = new double[M, N + 1];
            var objectiveVariable = new double[M];
            var mt = new MersenneTwister(seed);
            var normalDist = new Normal(0, 1, mt);
            for (int m = 0; m < M; m++)
            {
                var bm = new double[N];
                normalDist.Samples(bm);
                var path = GeneratePathByEM(bm);
                for (int n = 0; n < N + 1; n++)
                {
                    x[m, n] = path[n];
                }
            }
            for (int m = 0; m < M; m++)
            {
                objectiveVariable[m] = EuropeanPayoff(x[m, N], 0);
            }
            var explanatoryVariable = new double[M];
            for (int m = 0; m < N; m++)
            {
                explanatoryVariable[m] = S0;
            }
            int degree = 3;
            double[] coeff = Fit.Polynomial(explanatoryVariable, objectiveVariable, degree);
            double expectation = 0;
            for (int m = 0; m < M; m++)
            {
                expectation += coeff[0] + coeff[1] * explanatoryVariable[m] + coeff[2] * explanatoryVariable[m] * explanatoryVariable[m] + coeff[3] * explanatoryVariable[m] * explanatoryVariable[m] * explanatoryVariable[m];
            }
            return expectation / M;
        }
    }
}