using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PricingLogic
{
    internal class Program
    {
        static void Main(string[] args)
        {
            double s0 = 125;
            double r = 0.02;
            double sigma = 0.3;
            double T = 2;
            double K = 100;

            var analytic = new AnalyticBS(s0, r, sigma, T, K);
            double analyticPrice = analytic.AnalyticBSFormula();

            int N = 100;

            int M = (int)1e6;
            var lsm = new Lsm(s0, r, sigma, T, K);
            double lsmPrice = lsm.LeastSquareMonte(N, M);

            double lsmError = Math.Abs(lsmPrice - analyticPrice);
            Console.WriteLine($"analytic price : {analyticPrice}, lsm price : {lsmPrice}");
            Console.WriteLine($"error : {lsmError}");
        }
    }
}
