using FiniteMixtureModel.Distribution;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteMixtureModel.EM
{
    public class EMAlgorithm
    {
        List<Gaussian> gaussians;
        public double[] weights { get; set; }
        public double[,] posterior { get; set; }
        List<double> data;

        public EMAlgorithm(List<double> data, int component)
        {
            this.data = data;
            weights = new double[component];
            gaussians = new List<Gaussian>();
            posterior = new double[data.Count, component];

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
                weights[j] /= total;
            }
            for (int j = 0; j < component; j++)
            {
                for (int i = 0; i < data.Count; i++)
                {
                    if (means.Map[i] == j)
                    {
                        posterior[i, j] = 1;
                    }
                    else
                    {
                        posterior[i, j] = 0;
                    }
                }
            }

            // initial models
            for (int j = 0; j < component; j++)
            {
                gaussians.Add(new Gaussian(0, 1));
            }
            MStep();

        }

        public List<Gaussian> Models()
        { return gaussians;  }

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
                for (int j = 0; j < gaussians.Count; j++)
                {
                    posterior[i, j] /= total;
                }
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
        }

        public double LogLikelihood()
        {
            double ll = 0;
            for (int i = 0; i < data.Count; i++)
            {
                double sum = 0;
                for(int j = 0; j < gaussians.Count; j++)
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
    }
}
