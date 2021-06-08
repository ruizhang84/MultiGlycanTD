using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteMixtureModel.EM
{
    // Kmeans for 1-d Data
    public class KMeans
    {
        public int K { get; set; }
        public Dictionary<int, int> Map { get; set; } // x ==> cluster
        public Dictionary<int, List<double>> Clusters {get; set;}  // cluster => [x1, x2,...]

        double[] cluster;
        int maxIter = 1000;
        double tol = 0.01;

        public KMeans(int k)
        {
            K = k;
            cluster = new double[K];
            Map = new Dictionary<int, int>();
            Clusters = new Dictionary<int, List<double>>();
        }

        public void Run(List<double> data)
        {
            // init center
            data.Sort();
            int gap = data.Count / K;
            for(int i = 0; i < K; i++)
            {
                cluster[i] = data[gap * i];
            }

            // iteration
            for (int i = 0; i < maxIter; i++)
            {
                double diff = Iteration(data);
                if (diff < tol)
                    break;
            }
        }

        private double Iteration(List<double> data)
        {
            Clusters.Clear();
            for (int i = 0; i < data.Count; i++)
            {
                double minDistance = int.MaxValue;
                int index = 0;
                for(int j = 0; j < cluster.Count(); j++)
                {
                    if (Math.Abs(data[i] - cluster[j]) < minDistance)
                    {
                        minDistance = Math.Abs(data[i] - cluster[j]);
                        index = j;
                    }
                }
                Map[i] = index;
                if (!Clusters.ContainsKey(index))
                {
                    Clusters[index] = new List<double>();
                }
                Clusters[index].Add(data[i]);
            }

            // update cluster
            double diff = 0;
            for (int i = 0; i < cluster.Count(); i++)
            {
                if (Clusters.ContainsKey(i))
                {
                    double avg = Clusters[i].Average();
                    diff += Math.Abs(cluster[i] - avg);
                    cluster[i] = avg;
                }    
            }
            return diff;
        }
    }
}
