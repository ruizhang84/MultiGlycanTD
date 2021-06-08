using FiniteMixtureModel.Distribution;
using FiniteMixtureModel.EM;
using FiniteMixtureModel.RootFinding;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteMixtureModel.ECM
{
    public class ECMAlgorithm
    {
        double[,] Z;
        double[] lambda;
        List<Gamma> gammas;
        List<double> data;
        EMAlgorithm em;

        public ECMAlgorithm(List<double> data, int component)
        {
            this.data = data;
            lambda = new double[component];
            Z = new double[data.Count, component];
            gammas = new List<Gamma>();
            
            Init(data, component);
        }

        private void Init(List<double> data, int component)
        {
            // init lambda
            em = new EMAlgorithm(
                data.
                Select(x => Math.Pow(x, 1.0 / 3))
                .ToList(), component);
            em.Run();
            Z = em.posterior.Clone() as double[,];
            lambda = em.weights.Clone() as double[];
            // weighted method-of-moments(MOM) estimators
            for (int j = 0; j < component; j++)
            {
                double mean = 0;
                double sum = 0;
                double variance = 0;
                for (int i = 0; i < data.Count; i++)
                {
                    mean += data[i] * Z[i, j];
                    sum += Z[i, j];
                }
                mean /= sum;
                for (int i = 0; i < data.Count; i++)
                {
                    variance += (data[i] - mean) * (data[i] - mean) * Z[i, j];
                }
                variance = data.Count / ((data.Count - 1) * sum);
                gammas.Add(new Gamma(
                    mean * mean / variance,
                    variance / mean
                    ));
            }
        }


        public List<Gamma> Models()
        { return gammas; }

        public bool Run(int maxIter = 100, double tol = 0.01)
        {
            double ll = 0;
            for (int i = 0; i < maxIter; i++)
            {
                EStep();
                CMStep();
                double ll_new = LogLikelihood();
                if (Math.Abs(ll - ll_new) < tol)
                    break;
                ll = ll_new;
                for (int j = 0; j < gammas.Count; j++)
                {
                    // the number of component is setting too much.
                    if (lambda[j] < 0.001)
                        return false;
                }
            }
            return true;
        }

        public void EStep()
        {
            for(int i = 0; i < data.Count; i++)
            {
                double sum = 0;
                for(int j = 0; j < gammas.Count; j++)
                {
                    Z[i, j] = lambda[j] * gammas[j].ProbabilityDensity(data[i]);
                    sum += Z[i, j];
                }
                for(int j = 0; j < gammas.Count; j++)
                {
                    Z[i, j] /= sum;
                }
            }
        }

        public void CMStep()
        {
            // update lambda_j
            for(int j = 0; j < gammas.Count; j++)
            {
                double sum = 0;
                for(int i = 0; i < data.Count; i++)
                {
                    sum += Z[i, j];
                }
                lambda[j] = sum / data.Count;
            }
            // solve for alpha_j
            for(int j = 0; j < gammas.Count; j++)
            {
                double alpha = NewtonRaphson.Find(
                    x =>
                    {
                        double sum = 0;
                        double sumx = 0;
                        double sumlogx = 0;
                        for (int i = 0; i < data.Count; i++)
                        {
                            sum += Z[i, j];
                            sumx += Z[i, j] * data[i];
                            sumlogx += Z[i, j] * Math.Log(data[i]);
                        }
                        return Math.Log(x) - Digamma.Value(x)
                            - Math.Log(sumx / sum) + sumlogx / sum;
                    },
                    x =>
                    {
                        return 1.0 / x - Trigamma.Value(x);
                    },
                    gammas[j].Alpha, 0.0001);
                gammas[j].Alpha = alpha;
            }
            // update beta_j
            for(int j = 0; j < gammas.Count; j++)
            {
                double mean = 0;
                double sum = 0;
                for(int i = 0; i < data.Count; i++)
                {
                    mean += Z[i, j] * data[i];
                    sum += Z[i, j];
                }
                mean /= sum;
                gammas[j].Beta = mean / gammas[j].Alpha;
            }
        }

        public double LogLikelihood()
        {
            double ll = 0;
            for (int i = 0; i < data.Count; i++)
            {
                double sum = 0;
                for (int j = 0; j < gammas.Count; j++)
                {
                    sum += lambda[j] * gammas[j].ProbabilityDensity(data[i]);
                }
                ll += Math.Log(sum);
            }
            return ll;
        }

        public double BIC()
        {
            return gammas.Count * Math.Log(data.Count) - 2 * LogLikelihood();
        }

        public double ProbabilityDensity(double x)
        {
            double prob = 0;
            for (int j = 0; j < gammas.Count; j++)
            {
                prob += lambda[j] * gammas[j].ProbabilityDensity(x);
            }
            return prob;
        }

    }
}
