using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;

namespace AOA
{
    class Aquila : IOptimizationAlgorithm
    {
        public string Name { get; set; }

        public double[] XBest { get; set; }

        public double FBest { get; set; }

        public int NumberOfEvaluationFitnessFunction { get; set; }
        private double[][] range;
        private double[][] population;
        private double[] functionValue;
        private int populationSize;
        private int dimentions;
        private int iterations;
        private double alpha;
        private double beta;
        private double delta;
        private double omega;
        private int D;
        private double s;
        private double U;
        private int r1;
        public int OOO = 0;
        private double x;
        private double y;
        private double sigma;
        private Random rand;
        Func<double[], double> testFunction;

        private int rcount;
        private double r;

        public Aquila(Func<double[], double> myFuncion, int dim_in, int iterations_in, double[][] range_in, int population_size_in, double alpha_in, double beta_in,
            double delta_in, double omega_in, int D_in, double U_in, double s_in, int r1_in)
        {
            Name = "AOA";
            populationSize = population_size_in;
            dimentions = dim_in;
            iterations = iterations_in;
            testFunction = myFuncion;
            range = range_in;
            alpha = alpha_in;
            beta = beta_in;
            delta = delta_in;
            omega = omega_in;
            rcount = 0;
            U = U_in;
            s = s_in;
            r1 = r1_in;
            OOO = 0;
            sigma = ((SpecialFunctions.Gamma(1 + beta)) * Math.Sin(Math.PI * beta / 2)) / ((SpecialFunctions.Gamma((1 + beta) / 2)) * beta * Math.Pow(2, (beta - 1) / 2));
            rand = new Random();

        }

        private void CreatePopulation()
        {
            population = new double[populationSize][];
            for (int i = 0; i < populationSize; i++)
            {
                population[i] = new double[dimentions];
            }

            for (int i = 0; i < populationSize; i++)
            {
                for (int j = 0; j < dimentions; j++)
                {
                    population[i][j] = rand.NextDouble() * (range[j][1] - range[j][0]) + range[j][0];
                }
            }
        }

        private void UpdatePositions(int iteration)
        {
            double[] newPosition = new double[dimentions];

            double[] meanDimentions = new double[dimentions];
            for (int i = 0; i < dimentions; i++)
            {
                meanDimentions[i] = 0.0;
                for (int j = 0; j < populationSize; j++)
                {
                    meanDimentions[i] += population[j][i];
                }
                meanDimentions[i] /= populationSize;
            }

            if (iteration < iterations * 2 / 3)
            {
                if (rand.NextDouble() <= 0.5)
                {
                   
                    for (int i = 0; i < populationSize; i++)
                    {
                        for (int j = 0; j < dimentions; j++)
                        {
                            int do_control = 0;
                            do
                            {
                                newPosition[j] = XBest[j] * (1 - (iteration / iterations)) + (meanDimentions[j] - XBest[j] * rand.NextDouble());
                                do_control++;
                            } while (!(range[j][0] <= newPosition[j] && newPosition[j] <= range[j][1]) && do_control < 1000);
                            if (do_control == 1000)
                            {
                                newPosition[j] = meanDimentions[j];
                                //OOO += 1;
                            }



                        }

                        if (testFunction(newPosition) < functionValue[i])
                        {
                            population[i] = newPosition;
                            functionValue[i] = testFunction(newPosition);


                            if (testFunction(newPosition) < FBest)
                            {
                                XBest = (double[])newPosition.Clone();
                                FBest = testFunction(XBest);
                            }
                        }
                    }


                }
                else
                {
                    
                    int r1 = rand.Next(20);
                    double[] levyD = new double[dimentions];
                    for (int j = 0; j < dimentions; j++)
                    {
                         levyD[j] = s * rand.NextDouble() * sigma / Math.Pow(rand.NextDouble(), 1 / beta);
                    }
                    if (rcount == 0)
                    {
                        r = r1 + U * rand.Next(dimentions);
                        rcount = rand.Next(10);
                    }
                    rcount--;
                    for (int i = 0; i < populationSize; i++)
                    {
                        double theta = -omega * rand.Next(dimentions) + ((3 * Math.PI) / 2);
                        
                        x = r * Math.Sin(theta);
                        y = r * Math.Cos(theta);
                       // difrent for every dim ???? 
                        int randomPopulation = rand.Next(populationSize);
                        for (int j = 0; j < dimentions; j++)
                        {
                            int do_control = 0;
                            do
                            {
                                newPosition[j] = XBest[j] * levyD[j] + population[randomPopulation][j] * ((rand.Next(20) / 10) - 1) +  (y - x) * ((rand.Next(20)/10)-1);
                                do_control++;
                            } while (!(range[j][0] <= newPosition[j] && newPosition[j] <= range[j][1]) && do_control < 1000);
                            if (do_control == 1000)
                            {
                                newPosition[j] = meanDimentions[j];
                                OOO += 1;
                            }
                        }

                        if (testFunction(newPosition) < functionValue[i])
                        {
                            population[i] = newPosition;
                            functionValue[i] = testFunction(newPosition);


                            if (testFunction(newPosition) < FBest)
                            {
                                XBest = (double[])newPosition.Clone();
                                FBest = testFunction(XBest);
                            }
                        }
                    }
                }
            }
            else
            {
                if (rand.NextDouble() <= 0.5)
                {
                    

                    for (int i = 0; i < populationSize; i++)
                    {
                        for (int j = 0; j < dimentions; j++)
                        {
                            int do_control = 0;
                            do
                            {
                                do_control++;
                                newPosition[j] = (XBest[j] - meanDimentions[j]) * alpha - rand.NextDouble() + ((range[j][1] - range[j][0]) * rand.NextDouble() + range[j][0]) * delta;
                            } while (!(range[j][0] <= newPosition[j] && newPosition[j] <= range[j][1]) && do_control < 1000);
                            if (do_control == 1000)
                            {
                                newPosition[j] = meanDimentions[j];
                                //OOO += 1;
                            }
                        }

                        if (testFunction(newPosition) < functionValue[i])
                        {
                            population[i] = newPosition;
                            functionValue[i] = testFunction(newPosition);


                            if (testFunction(newPosition) < FBest)
                            {
                                XBest = (double[])newPosition.Clone();
                                FBest = testFunction(XBest);
                            }
                        }
                    }
                }
                else
                {
                    double qf = Math.Pow(iteration, (2 * rand.NextDouble() - 1) / Math.Pow(1 - iterations, 2));
                    int r1 = rand.Next(20);

                    double g2 = 2 * (1 - iteration / iterations);

                    for (int i = 0; i < populationSize; i++)
                    {
                        for (int j = 0; j < dimentions; j++)
                        {
                            int do_control = 0;
                            do
                            {
                                double g1 = 2 * rand.NextDouble() - 1;

                                double theta = -omega * rand.Next(dimentions) + ((3 * Math.PI) / 2);
                                double r = r1 + U * rand.Next(dimentions);
                                x = r * Math.Sin(theta);
                                y = r * Math.Cos(theta);
                                double levy = s * rand.NextDouble() * sigma / Math.Pow(rand.NextDouble(), 1 / beta); // difrent for every dim ???? 
                                newPosition[j] = qf * XBest[j] - g1 * population[i][j] * rand.NextDouble() - g2 * levy + rand.NextDouble() * g1;
                                do_control++;
                            } while (!(range[j][0] <= newPosition[j] && newPosition[j] <= range[j][1]) && do_control < 1000);
                            if (do_control == 1000)
                            {
                                newPosition[j] = meanDimentions[j];
                                //OOO += 1;
                            }
                        }

                        if (testFunction(newPosition) < functionValue[i])
                        {
                            population[i] = newPosition;
                            functionValue[i] = testFunction(newPosition);


                            if (testFunction(newPosition) < FBest)
                            {
                                XBest = (double[])newPosition.Clone();
                                FBest = testFunction(XBest);
                            }
                        }
                    }
                }
            }
        }

        public double Solve()
        {
            CreatePopulation();
            XBest = population[0];
            FBest = testFunction(population[0]);
            functionValue = new double[populationSize];
            for (int i = 0; i < populationSize; i++)
            {
                functionValue[i] = testFunction(population[i]);
            }
            for (int i = 0; i < iterations; i++)
            {
                //Console.WriteLine(i < iterations * 2 / 3);
                //Console.WriteLine(FBest);
                //foreach (double d in XBest)
                {
                  //  Console.Write(" " + d.ToString());
                }
                //Console.WriteLine();
                UpdatePositions(i);


            }
            return FBest;
        }
    }
}