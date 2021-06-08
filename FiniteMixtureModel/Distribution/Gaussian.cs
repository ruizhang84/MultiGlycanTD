using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteMixtureModel.Distribution
{
    public class Gaussian
    {
        public double Mean { get; set; }
        public double Std { get; set; }
        private Random rand = new Random();

        public Gaussian(double mean, double std)
        {
            Mean = mean;
            Std = std;
        }

        private double Square(double x)
        { return x * x; }

        public double ProbabilityDensity(double x)
        {
            return 1 / (Std * Math.Sqrt(2 * Math.PI))
                * Math.Exp(-0.5 * Square((x - Mean) / Std));
        }

        // Box-Muller transform 
        // cite: https://stackoverflow.com/questions/218060/random-gaussian-variables
        public double Sample()
        {
            double u1 = 1.0 - rand.NextDouble(); 
            double u2 = 1.0 - rand.NextDouble();
            double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                         Math.Sin(2.0 * Math.PI * u2);
            double randNormal =
                         Mean + Std * randStdNormal; //random normal(mean,stdDev^2)
            return randNormal;
        }

    }
}
