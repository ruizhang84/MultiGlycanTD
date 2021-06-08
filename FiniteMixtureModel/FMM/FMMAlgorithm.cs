using FiniteMixtureModel.Distribution;
using FiniteMixtureModel.ECM;
using FiniteMixtureModel.EM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteMixtureModel.FMM
{
    public class FMMAlgorithm
    {
        List<Gaussian> gaussians;
        ECMAlgorithm gamma;
        double[] weights;
        double weight0;
        double[,] posterior;
        double[] posterior0;
        List<double> data;

        public FMMAlgorithm(List<double> data, int component, ECMAlgorithm gamma, double weight0)
        {
            this.data = data;
            this.gamma = gamma;
            this.weight0 = weight0;

            weights = new double[component];
            gaussians = new List<Gaussian>();
            posterior = new double[data.Count, component];
            posterior0 = new double[data.Count];

            Init(data, component);
        }

        private void Init(List<double> data, int component)
        {
            // k-means on data
            KMeans means = new KMeans(component);
            means.Run(data);

            // update weights;
            double total = 0;
            for (int j = 0; j < component; j++)
            {
                if (means.Clusters.ContainsKey(j))
                {
                    weights[j] = means.Clusters[j].Count + 1;
                }
                else
                {
                    weights[j] = 1;
                }
                total += weights[j];
            }
            for (int j = 0; j < component; j++)
            {
                weights[j] /= total * (1 - weight0);
            }
            for (int j = 0; j < component; j++)
            {
                for (int i = 0; i < data.Count; i++)
                {
                    if (means.Map[i] == j)
                    {
                        posterior[i, j] = 1 - weight0;
                    }
                    else
                    {
                        posterior[i, j] = 0;
                    }
                }
            }
            for (int i = 0; i < data.Count; i++)
            {
                posterior0[i] = weight0;
            }

            // initial models
            for (int j = 0; j < component; j++)
            {
                gaussians.Add(new Gaussian(0, 1));
            }
            MStep();
        }

        public List<Gaussian> Models()
        { return gaussians; }

        public void Run(int maxIter = 100, double tol = 0.01)
        {
            double ll = 0;
            for (int i = 0; i < maxIter; i++)
            {
                EStep();
                MStep();
                double ll_new = LogLikelihood();
                if (Math.Abs(ll - ll_new) < tol)
                    break;
                ll = ll_new;
            }
        }

        public void EStep()
        {
            for (int i = 0; i < data.Count; i++)
            {
                double total = 0;
                for (int j = 0; j < gaussians.Count; j++)
                {
                    posterior[i, j] = weights[j] * gaussians[j].ProbabilityDensity(data[i]);
                    total += posterior[i, j];
                }
                posterior0[i] = weight0 * gamma.ProbabilityDensity(data[i]);
                total += posterior0[i];
                for (int j = 0; j < gaussians.Count; j++)
                {
                    posterior[i, j] /= total;
                }
                posterior0[i] /= total;
            }
        }

        public void MStep()
        {
            for (int j = 0; j < gaussians.Count; j++)
            {
                double mean = 0;
                double z = 0;
                for (int i = 0; i < data.Count; i++)
                {
                    mean += posterior[i, j] * data[i];
                    z += posterior[i, j];
                }
                gaussians[j].Mean = mean / z;
                weights[j] = z / data.Count;

                double std = 0;
                for (int i = 0; i < data.Count; i++)
                {
                    std += posterior[i, j] * (data[i] - gaussians[j].Mean) * (data[i] - gaussians[j].Mean);
                }
                gaussians[j].Std = Math.Sqrt(std / z);
            }
            weight0 = posterior0.Average();
        }

        public double LogLikelihood()
        {
            double ll = 0;
            for (int i = 0; i < data.Count; i++)
            {
                double sum = weight0 * gamma.ProbabilityDensity(data[i]);
                for (int j = 0; j < gaussians.Count; j++)
                {
                    sum += weights[j] * gaussians[j].ProbabilityDensity(data[i]);
                }
                ll += Math.Log(sum);
            }
            return ll;
        }

        public double BIC()
        {
            return gaussians.Count * Math.Log(data.Count) - 2 * LogLikelihood();
        }

        public double Probability(double x)
        {
            double prob = weight0 * gamma.ProbabilityDensity(x);
            for (int j = 0; j < gaussians.Count; j++)
            {
                prob += weights[j] * gaussians[j].ProbabilityDensity(x);
            }
            return weight0 * gamma.ProbabilityDensity(x) / prob;
        }
    }
}
