using FiniteMixtureModel.LanczosApproximation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteMixtureModel.Distribution
{
    public class Gamma
    {
        private double alpha;
        private double gamma;
        public double Alpha
        {
            get
            {
                return alpha;
            }
            set
            {
                gamma = GammaCalculator.Compute(alpha);
                alpha = value;
            }
        }
        public double Beta { get; set; }

        public Gamma(double alpha, double beta)
        {
            Alpha = alpha;
            Beta = beta;
            gamma = GammaCalculator.Compute(alpha);
        }

        public double ProbabilityDensity(double x)
        {
            double prob = Math.Pow(x, Alpha - 1)
                * Math.Exp(-x / Beta) / (Math.Pow(Beta, Alpha) * gamma);
            if (double.IsNaN(prob))
                return 0;
            return prob;
        }

    }
}
