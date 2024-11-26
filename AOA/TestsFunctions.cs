using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AOA
{
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

        public static double[][] SphereLimits(int dim)
        {
            double[][] limits = new double[dim][];
            for (int i = 0; i < dim; i++)
            {
                limits[i] = new double[] { -10, 10 };
            }
            return limits;
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

        public static double[][] RastriginLimits(int dim)
        {
            double[][] limits = new double[dim][];
            for (int i = 0; i < dim; i++)
            {
                limits[i] = new double[] { -4.5, 4.5 };
            }
            return limits;
        }

        public static double Beale(double[] xi)
        {
            return Math.Pow((1.5 - xi[0] + xi[0] * xi[1]), 2) + Math.Pow((2.25 - xi[0] + xi[0] * Math.Pow(xi[1], 2)), 2) +
                Math.Pow((2.625 - xi[0] + xi[0] * Math.Pow(xi[1], 3)), 2);
        }
        public static double[][] BealeLimits()
        {
            double[][] limits = new double[2][];
            for (int i = 0; i < 2; i++)
            {
                limits[i] = new double[] { -5.12, 5.12 };
            }
            return limits;
        }


        public static double Rosenbrock(double[] xi)
        {
            double ret = 0;
            for (int i = 0; i < xi.Length - 1; i++)
            {
                ret += (100 * Math.Pow(xi[i + 1] - Math.Pow(xi[i], 2), 2) + Math.Pow(xi[i] - 1, 2));
            }
            return ret;
        }

        public static double[][] RosenbrockLimits(int dim)
        {
            double[][] limits = new double[dim][];
            for (int i = 0; i < dim; i++)
            {
                limits[i] = new double[] { -10, 10 };
            }
            return limits;
        }

        public static double BukinN6(double[] xi)
        {
            return 100 * Math.Sqrt(Math.Abs(xi[1] - 0.01 * xi[0] * xi[0])) + 0.01 * Math.Abs(xi[0] + 10);
        }

        public static double[][] BukinN6Limits()
        {
            double[][] limits = new double[2][];
            limits[0] = new double[] { -15, -5 };
            limits[1] = new double[] { -3, 3 };
            return limits;
        }

        public static double HimmelblauN6(double[] xi)
        {
            return Math.Pow(xi[0] * xi[0] + xi[1] - 11, 2) + Math.Pow(xi[0] + xi[1] * xi[1] - 7, 2);
        }

        public static double[][] HimmelblauN6Limits()
        {
            double[][] limits = new double[2][];
            for (int i = 0; i < 2; i++)
            {
                limits[i] = new double[] { -5, 5 };
            }
            return limits;
        }
    }
}
