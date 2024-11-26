using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;

namespace AOA
{
    class AquilaOptimizer: IOptimizationAlgorithm
    {
        public double[] XBest { get; set; }
        public string Name { get; set; }
        public double FBest{ get; set; }
        public int NumberOfEvaluationFitnessFunction { get; set; }
        private Func<double[], double> function;
        private int dim;
        private double[][] limits;
        private int population_size;
        private int iterations;
        private double alpha;
        private double delta;
        private double beta;
        private double[][] population;
        private Random rand;
        public AquilaOptimizer(Func<double[], double> function, int dim, double[][] limits, int population_size, int iterations, double alpha, double delta, double beta)
        {
            this.function = function;
            this.dim = dim;
            this.limits = limits;
            this.population_size = population_size;
            this.iterations = iterations;
            this.alpha = alpha;
            this.delta = delta;
            this.beta = beta;
            this.rand = new Random();
        }

        private void CreatePopulation()
        {
            population = new double[population_size][];
            for (int i = 0; i < population_size; i++)
            {
                population[i] = new double[dim];
                population[i] = Enumerable.Range(0, dim).Select(j => rand.NextDouble() * (limits[j][1] - limits[j][0]) + limits[j][0]).ToArray();
            }
        }

        private void FindBest()
        {
            for (int i = 0; i < population_size; i++)
            {
                if (function(population[i]) < FBest)
                {
                    FBest = function(population[i]);
                    XBest = (double[])population[i].Clone();
                }
            }
        }

        private double[] Levy(int dim)
        {
            double sigma = Math.Pow(
            (SpecialFunctions.Gamma(1 + beta) * Math.Sin(Math.PI * beta / 2))
            / (SpecialFunctions.Gamma((1 + beta) / 2) * beta * Math.Pow(2, (beta - 1) / 2)),
            1 / beta);
            double[] u = new double[dim];
            double[] v = new double[dim];

            for (int i = 0; i < dim; i++)
            {
                u[i] = Math.Sqrt(-2.0 * Math.Log(rand.NextDouble())) * Math.Cos(2.0 * Math.PI * rand.NextDouble()) * sigma;
                v[i] = Math.Sqrt(-2.0 * Math.Log(rand.NextDouble())) * Math.Cos(2.0 * Math.PI * rand.NextDouble());
            }

            double[] step = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                step[i] = u[i] / Math.Pow(Math.Abs(v[i]), 1 / beta);
            }

            return step;
        }

        private void updatePositions(int iteration)
        {
            double g2 = 2 * rand.NextDouble() - 1;
            double g1 = 2 * (1 - iteration / iterations);
            int[] d1 = Enumerable.Range(1, dim).ToArray();
            double u = 0.0265;
            int r0 = 10;
            double[] r = d1.Select(d => r0 + u * d).ToArray();
            double omega = 0.005;
            double phi0 = 3 * Math.PI / 2;
            double[] phi = d1.Select(d => -omega * d + phi0).ToArray();
            double[] x = phi.Zip(r, (pi, ri) => ri * Math.Sin(pi)).ToArray();
            double[] y = phi.Zip(r, (pi, ri) => ri * Math.Cos(pi)).ToArray();
            double qf = Math.Pow(iteration, ((2 * rand.NextDouble() - 1) / (Math.Pow(1 - iterations, 2))));
            double[][] new_population = new double[population_size][];
            for (int i = 0; i < population_size; i++)
            {
                new_population[i] = new double[dim];
                if (iteration < iterations * 2 / 3)
                {
                    if (rand.NextDouble() < 0.5)
                    {
                        new_population[i] = XBest.Select(xi => xi * (1 - iteration / iterations) + (population[i].Average() - xi) * rand.NextDouble()).ToArray();
                    }
                    else
                    {
                        new_population[i] = Enumerable.Range(0, XBest.Length).Select(j => XBest[j] * Levy(dim)[j] + population[rand.Next(population_size)][j] +
                        (x[j] - y[j]) * rand.NextDouble()).ToArray();

                    }
                }
                else
                {
                    if (rand.NextDouble() < 0.5)
                    {
                        new_population[i] = limits.Zip(XBest, (li, xi) => (xi - population[i].Average()) * alpha - rand.NextDouble() +
                        ((li[1] - li[0]) * rand.NextDouble() - li[0]) * delta).ToArray();
                    }
                    else
                    {
                        new_population[i] = XBest.Zip(population[i].Zip(Levy(dim), (pi, li) => g2 * pi * rand.NextDouble() + g1 * li).ToArray()
                            , (xi, lp) => qf * xi - lp + rand.NextDouble() * g2).ToArray();
                    }
                }

            }
            for (int k = 0; k < population_size; k++)
            {
                if (function(population[k]) > function(new_population[k]))
                {
                    population[k] = (double[])new_population[k].Clone();
                }
            }

        }
        public double Solve()
        {
            CreatePopulation();
            XBest = (double[])population[0].Clone();
            FBest = function(XBest);
            for (int iter = 0; iter < iterations; iter++)
            {
                updatePositions(iter);
                FindBest();
            }
            return FBest;
        }
    }
}
