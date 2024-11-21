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
            double[] xbest = new double[5];
            double fbest = 0.0;
            double ooo = 0;
            for (int j = 0; j < 5; j++)
            {
                xbest[j] = 0.0;
            }
            int tests = 1000;
            int ttt = tests;
            for (int iter = 0; iter < tests; iter++)
            {
                Aquila a = new Aquila(myFuncion: TestsFunctions.Beale, dim_in: 2, iterations_in: 100, range_in: range, population_size_in: 100,
                    alpha_in: 0.1, beta_in: 1.5, delta_in: 0.1, omega_in: 0.005, D_in: 3, U_in: 0.00565, s_in: 0.01, r1_in: 1);
                a.Solve();
                fbest += a.FBest;
                ooo += a.OOO;
                if (a.OOO <= 0) ttt--;
                for (int j = 0; j < 2; j++)
                {
                    xbest[j] += a.XBest[j];
                }

            }
            Console.WriteLine("Arg Fbest from " + tests.ToString()+" tests :" + (fbest / tests).ToString());
            //Console.WriteLine("OOO: " + (ooo / tests));
            for (int j = 0; j < 5; j++)
            {
                Console.WriteLine("Arg Xbest["+j.ToString()+"] from " + tests.ToString() + " tests :" + (xbest[j] / tests).ToString());
            }

            //Console.WriteLine(a.FBest);
            //foreach(double e in a.XBest)
            {
                //  Console.WriteLine(e);
            }
            Console.Read();
        }
    }
}
