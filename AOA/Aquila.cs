using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


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
        private double x;
        private double y;
        private double sigma;
        private Random rand;
        Func<double[], double> testFunction;

        static double Gamma(double z)
        {
            // Lanczos' approximation coefficients
            double[] p = {
            0.99999999999980993,  676.5203681218851, -1259.1392167224028,
            771.32342877765313, -176.61502916214059,  12.507343278686905,
            -0.1385710331296526,  9.9843695780195716e-6,  1.5056327351493116e-7
        };

            if (z < 0.5)
                return Math.PI / (Math.Sin(Math.PI * z) * Gamma(1 - z));

            z -= 1;
            double x = p[0];
            for (int i = 1; i < p.Length; i++)
            {
                x += p[i] / (z + i);
            }

            double t = z + p.Length - 0.5;
            return Math.Sqrt(2 * Math.PI) * Math.Pow(t, z + 0.5) * Math.Exp(-t) * x;
        }

        public Aquila(Func<double[], double> myFuncion, int dim_in, int iterations_in, double[][] range_in, int population_size_in, double alpha_in, double beta_in,
            double delta_in, double omega_in, int D_in, double U_in, double s_in,int r1_in)
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
            U = U_in;
            D = D_in;
            s = s_in;
            r1 = r1_in;
            double theta = -omega_in * D + (3 * Math.PI / 2);
            double r = r1 + U * D;
            x = r * Math.Sin(theta);
            y = r * Math.Cos(theta);
            sigma = ((Gamma(1+ beta)) * Math.Sin(Math.PI * beta / 2)) / ((Gamma((1+beta/2))) * beta * Math.Pow(2,(beta - 1)/2));
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
          
                if (iteration < iterations*2/3)
            {
                if (rand.NextDouble() <= 0.5)
                {
                    double[] meanDimentions = new double[dimentions];
                    for (int i = 0; i < dimentions; i++)
                    {
                        meanDimentions[i] = 0;
                        for (int j = 0; j < populationSize; j++)
                        {
                            meanDimentions[i] += population[j][i];
                        }
                        meanDimentions[i] /= populationSize;
                    }
                    for (int i = 0; i < populationSize; i++)
                    {
                        for (int j = 0; j < dimentions; j++)
                        {

                            newPosition[j] = XBest[j] * (1 - iteration / iterations) + (meanDimentions[j] - XBest[j] * rand.NextDouble());
                        }

                        if (testFunction(newPosition) < functionValue[i])
                        {
                            population[i] = newPosition;
                            functionValue[i] = testFunction(newPosition);


                            if (testFunction(newPosition) < FBest)
                            {
                                XBest = newPosition;
                                FBest = testFunction(XBest);
                            }
                        }    
                    }
                    

                }
                else
                {
                    int randomPopulation = rand.Next(populationSize);
                   
                    for (int i = 0; i < populationSize; i++)
                    {
                        for (int j = 0; j < dimentions; j++)
                        {
                            
                            double levy = s * rand.NextDouble() * sigma / Math.Pow(rand.NextDouble(),1/beta); // difrent for every dim ???? 
                            newPosition[j] = XBest[j] * levy + population[randomPopulation][j] + (y-x) * rand.NextDouble();
                        }

                        if (testFunction(newPosition) < functionValue[i])
                        {
                            population[i] = newPosition;
                            functionValue[i] = testFunction(newPosition);


                            if (testFunction(newPosition) < FBest)
                            {
                                XBest = newPosition;
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
                    double[] meanDimentions = new double[dimentions];
                    for (int i = 0; i < dimentions; i++)
                    {
                        meanDimentions[i] = 0;
                        for (int j = 0; j < populationSize; j++)
                        {
                            meanDimentions[i] += population[j][i];
                        }
                        meanDimentions[i] /= populationSize;
                    }

                    for (int i = 0; i < populationSize; i++)
                    {
                        for (int j = 0; j < dimentions; j++)
                        {

                            newPosition[j] = (XBest[j] - meanDimentions[j]) * alpha - rand.NextDouble() + ((range[j][1] - range[j][0]) * rand.NextDouble() + range[j][0]) * delta;
                        }

                        if (testFunction(newPosition) < functionValue[i])
                        {
                            population[i] = newPosition;
                            functionValue[i] = testFunction(newPosition);


                            if (testFunction(newPosition) < FBest)
                            {
                                XBest = newPosition;
                                FBest = testFunction(XBest);
                            }
                        }
                    }
                }
                else
                {
                    double g1 = 2 * rand.NextDouble() - 1;
                    double g2 = 2 * (1 - iteration/iterations);
                    double qf = Math.Pow(iteration, (2 * rand.NextDouble() - 1) / Math.Pow(1 - iterations, 2));
                    for (int i = 0; i < populationSize; i++)
                    {
                        for (int j = 0; j < dimentions; j++)
                        {
                            double levy = s * rand.NextDouble() * sigma / Math.Pow(rand.NextDouble(), 1 / beta); // difrent for every dim ???? 
                            newPosition[j] = qf * XBest[j] - g1 * XBest[j] * rand.NextDouble() - g2 * levy + rand.NextDouble() * g1;
                        }

                        if (testFunction(newPosition) < functionValue[i])
                        {
                            population[i] = newPosition;
                            functionValue[i] = testFunction(newPosition);


                            if (testFunction(newPosition) < FBest)
                            {
                                XBest = newPosition;
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
            for(int i =0; i < populationSize;i++)
            {
                functionValue[i] = testFunction(population[i]);
            }
            for (int i = 0; i < iterations; i++)
            {
                Console.WriteLine(i < iterations * 2 / 3);
                Console.WriteLine(FBest);
                foreach(double d in XBest)
                {
                    Console.Write(" " + d.ToString());
                }
                Console.WriteLine();
                UpdatePositions(i);
                
               
            }
            return FBest;
        }
    }
}
