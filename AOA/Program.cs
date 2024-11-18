using System;
namespace AOA
{
    class Program
    {
        static void Main(string[] args)
        {

            double[][] range = new double[5][];
            for (int t = 0; t < 5; t++)
                range[t] = new double[] { -5.12, 5.12 };

            Aquila a = new Aquila(myFuncion: TestsFunctions.Sphere, dim_in: 5, iterations_in: 100, range_in: range, population_size_in: 100,
                alpha_in: 0.1, beta_in: 1.0, delta_in: 0.1, omega_in: 0.001, D_in: 3, U_in: 0.00565, s_in: 0.01, r1_in: 1);
            a.Solve();
            Console.WriteLine(a.FBest);
            foreach (double e in a.XBest)
            {
                Console.WriteLine(e);
            }
        }
    }
}
