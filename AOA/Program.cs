using System;

namespace AOA
{
    class Program
    {
        static void Main(string[] args)
        {
            AquilaOptimizer ao = new AquilaOptimizer(TestsFunctions.Sphere, 20, TestsFunctions.RastriginLimits(20), 50, 50, 0.05, 0.05, 1.5);
            ao.Solve();
            Console.WriteLine(ao.FBest);
            foreach (double x in ao.XBest)
            {
                Console.WriteLine(x);
            }
            Console.ReadLine();
        }
    }
}
