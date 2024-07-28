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
            double s0 = 75;
            double r = 0.02;
            double sigma = 0.1;
            double T = 2;
            double K = 100;

            var analytic = new AnalyticBS(s0, r, sigma, T, K);
            double analyticPrice = analytic.AnalyticBSFormula();

            int N = 100;

            //int M = (int)1e6;
            //var lsm = new Lsm(s0, r, sigma, T, K);
            //double lsmPrice = lsm.LeastSquareMonte(N, M);

            //double lsmError = Math.Abs(lsmPrice - analyticPrice);
            //Console.WriteLine($"analytic price : {analyticPrice}, lsm price : {lsmPrice}");
            //Console.WriteLine($"error : {lsmError}");

            N = 5000;
            var fdm = new Fdm(s0, r, sigma, T, K);
            double fdmPrice = fdm.ImplicitFdm(N);
            double fdmError = Math.Abs(fdmPrice - analyticPrice);
            Console.WriteLine($"analytic price : {analyticPrice}, fdm price : {fdmPrice}");
            Console.WriteLine($"error : {fdmError}");
        }
    }
}
