using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AOA
{
    public class OptimizationResult
    {
        public string name { get; set; }
        public int testNumber { get; set; }
        public double a { get; set; }
        public double c { get; set; }
        public double p { get; set; }
        public int dimentions { get; set; }
        public int populationSize { get; set; }
        public int iterationNumber { get; set; }
        public double MeanResult { get; set; }
        public double StdResult { get; set; }
        public double[] MeanPos { get; set; }
        public double[] StdPos { get; set; }
        public double[] BestPositions { get; set; }

        public OptimizationResult()
        {
            BestPositions = new double[2];
            MeanPos = new double[2];
            StdPos = new double[2];
        }
    }


    class TestsFunctions
    { 

        public static double Sphere(double[] xi)
        {
            double ret = 0;
            foreach (double x in xi)
            {
                ret += Math.Pow(x, 2);
            }
            return ret;
        }

        public static double Rastrigin(double[] xi)
        {
            double ret = 0;
            foreach (double x in xi)
            {
                ret += (Math.Pow(x, 2) - 10 * Math.Cos(2 * Math.PI * x) + 10);
            }
            return ret;
        }

        public static double Rosenbrock(double[] xi)
        {
            return Math.Pow((1.5 - xi[0] + xi[0] * xi[1]), 2) + Math.Pow((2.25 - xi[0] + xi[0] * Math.Pow(xi[1], 2)), 2) +
                Math.Pow((2.625 - xi[0] + xi[0] * Math.Pow(xi[1], 3)), 2);
        }

        public static double Beale(double[] xi)
        {
            double ret = 0;
            for (int i = 0; i < xi.Length - 1; i++)
            {
                ret += (100 * Math.Pow(xi[i + 1] - Math.Pow(xi[i], 2), 2) + Math.Pow(xi[i] - 1, 2));
            }
            return ret;
        }

        public static double BukinN6(double[] xi)
        {
            return 100 * Math.Sqrt(Math.Abs(xi[1] - 0.01 * xi[0] * xi[0])) + 0.01 * Math.Abs(xi[0] + 10);
        }

        public static double HimmelblauN6(double[] xi)
        {
            return Math.Pow(xi[0] * xi[0] + xi[1] - 11, 2) + Math.Pow(xi[0] + xi[1] * xi[1] - 7, 2);
        }

    //    public static OptimizationResult FindBest(Func<double[], double> myFuncion, int dim_in, int iterations_in, double[][] range_in, int population_size_in,
    //        double[] a_list, double[] c_list, double[] p_list, string name,int test_nb)
    //    {
    //        OptimizationResult opt = new OptimizationResult();
    //        opt.MeanResult = 1000000;
    //        opt.name = name;
    //        opt.populationSize = population_size_in;
    //        opt.iterationNumber = iterations_in;
    //        opt.dimentions = dim_in;
    //        opt.testNumber = test_nb;
    //        foreach (double a in a_list)
    //        {
    //            foreach (double c in c_list)
    //            {
    //                foreach (double p in p_list)
    //                {
    //                    double[] results = new double[test_nb];
    //                    double[][] results_pos = new double[test_nb][];
    //                    for (int i =0; i < test_nb;i++)
    //                    {
    //                        results_pos[i] = new double[dim_in];
    //                    }

    //                    double mean = 0;
    //                    double[] mean_pos = new double[dim_in];
    //                    for (int i = 0; i < test_nb; i++)
    //                    {
    //                        Aquila but = new Aquila(myFuncion,dim_in,iterations_in,range_in,population_size_in,a,c,p);
    //                        but.Solve();
    //                        results[i] = but.FBest;
    //                        mean += but.FBest;
    //                        for (int j = 0; j < dim_in; j++)
    //                        {
    //                            results_pos[i][j] = but.XBest[j];
    //                            mean_pos[j] += but.XBest[j];
    //                        }
    //                    }
    //                    mean /= test_nb;
    //                    for (int j = 0; j < dim_in; j++)
    //                    {
    //                        mean_pos[j] /= test_nb;
    //                    }
    //                    double std = 0;
    //                    double[] std_pos = new double[dim_in];
    //                    for (int i = 0; i < test_nb; i++)
    //                    {
    //                        std += Math.Pow(results[i] - mean, 2);
    //                        for (int j = 0; j < dim_in; j++)
    //                        {
    //                            std_pos[j] = Math.Pow(results_pos[i][j] - mean_pos[j], 2);
    //                        }
    //                    }
    //                    std = Math.Sqrt(std / test_nb);
    //                    for (int j = 0; j < dim_in; j++)
    //                    {
    //                        std_pos[j] = Math.Sqrt(std_pos[j] / test_nb);
    //                    }
    //                    if (mean < opt.MeanResult)
    //                    {
    //                        opt.a = a;
    //                        opt.c = c;
    //                        opt.p = p;
    //                        opt.MeanResult = mean;
    //                        opt.StdResult = std;
    //                        opt.MeanPos = mean_pos;
    //                        opt.StdPos = std_pos;
    //                    }

    //                }
    //            }
    //        }
    //        Console.WriteLine("Name: " + opt.name.ToString());
    //        Console.WriteLine("Test number: " + opt.testNumber.ToString());
    //        Console.WriteLine("Population size:" + opt.populationSize.ToString());
    //        Console.WriteLine("Dimentions:" + opt.dimentions.ToString());
    //        Console.WriteLine("Iteration number:" + opt.iterationNumber.ToString());
    //        Console.WriteLine("A:" + opt.a.ToString());
    //        Console.WriteLine("C:" + opt.c.ToString());
    //        Console.WriteLine("P:" + opt.p.ToString());
    //        Console.WriteLine("Mean: " + opt.MeanResult.ToString());
    //        Console.WriteLine("Std: " + opt.StdResult.ToString());
    //        for (int j = 0; j < dim_in; j++)
    //        {
    //            Console.Write(" Mean position: " + opt.MeanPos[j].ToString());
    //        }
    //        Console.WriteLine();
    //        for (int j = 0; j < dim_in; j++)
    //        {
    //            Console.Write(" Std positon: " + opt.StdPos[j].ToString());
    //        }
    //        Console.WriteLine();
    //        Console.WriteLine();
    //        Console.WriteLine();
    //        Console.WriteLine();
    //        return opt;
    //    }

    }
}
