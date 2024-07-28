using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PricingLogic
{
    public class Fdm
    {
        public double S0 { get; private set; }
        public double InteRate { get; private set; }
        public double Vola { get; private set; }
        public double Maturity { get; private set; }
        public double Strike { get; private set; }

        public Fdm(double s0, double r, double sigma, double T, double K)
        {
            S0 = s0;
            InteRate = r;
            Vola = sigma;
            Maturity = T;
            Strike = K;
        }

        private double ReturnPayoffPV(double x, double t)
        {
            return Math.Max(x - Strike, 0) * Math.Exp(-InteRate * (Maturity - t));
        }
        private double TerminalCondition(double x)
        {
            return Math.Max(x - Strike, 0);
        }

        private double BoundaryConditionAtZero(double t)
        {
            return 0;
        }

        private double BoundaryConditionAtInfty(double X, double t)
        {
            return X - Math.Exp(-InteRate * (Maturity - t)) * Strike;
        }

        public double ImplicitFdm(int N, int x = 8)
        {
            double X = S0 * x;
            while (true)
            {
                if (X > Strike)
                {
                    break;
                }
                else
                {
                    x += 1;
                    X = S0 * x;
                }
            }

            int M = x * 100;
            double dx = X / M;
            double dt = Maturity / N;

            double[,] fSim = new double[N + 1, M + 1];

            for (int j = 0; j <= M; j++)
            {
                fSim[N, j] = TerminalCondition(j * dx);
            }

            for (int i = 0; i <= N; i++)
            {
                fSim[i, 0] = BoundaryConditionAtZero(i * dt);
                fSim[i, M] = BoundaryConditionAtInfty(X, i * dt);
            }

            double[] a = new double[M]; // a[0]は使わない
            double[] b = new double[M]; // b[0]は使わない
            double[] c = new double[M]; // c[0]は使わない
            double[] d = new double[M]; // d[0]は使わない

            for (int i = N - 1; i >= 0; i--)
            {
                for (int j = 1; j < M; j++)
                {
                    a[j] = -0.5 * dt * (Vola * Vola * j * j + InteRate * j);
                    b[j] = 1 + dt * (Vola * Vola * j * j + InteRate);
                    c[j] = -0.5 * dt * (Vola * Vola * j * j - InteRate * j);
                    d[j] = fSim[i + 1, j];
                }

                d[1] -= c[1] * fSim[i, 0];
                d[M - 1] -= a[M - 1] * fSim[i, M];

                double[] P = new double[M]; // P[0]は使わない
                double[] Q = new double[M]; // Q[0]は使わない

                P[1] = a[1] / b[1];
                Q[1] = d[1] / b[1];

                for (int j = 2; j < M; j++)
                {
                    P[j] = a[j] / (b[j] - c[j] * P[j - 1]);
                    Q[j] = (d[j] - c[j] * Q[j - 1]) / (b[j] - c[j] * P[j - 1]);
                }

                fSim[i, M - 1] = Q[M - 1];

                for (int j = M - 2; j > 0; j--)
                {
                    fSim[i, j] = -P[j] * fSim[i, j + 1] + Q[j];
                }
                for (int j = 0; j < M; j++)
                {
                    fSim[i, j] = Math.Max(fSim[i, j], ReturnPayoffPV(dx * j, dt * i));
                }
            }
            int index = (int)(S0 / dx);
            return fSim[0, index];
        }
    }
}
