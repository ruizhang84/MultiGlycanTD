using SpectrumProcess.algorithm;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.score
{
    public class ClusterKMeans<T>
    {
        public int K { get; set; }
        double[] Center { get; }

        // x_i => cluster_index, x_i is original index
        public Dictionary<int, int> Index { get; }
        // y_i => x_i, find index before sorting
        Dictionary<int, int> map_;

        // cluster_index => [y1, y2,...]
        public Dictionary<int, List<Point<T>>> Clusters { get; set; }

        public int MaxIter { get; set; }
        public double Tol { get; set; }

        public ClusterKMeans(int k = 3, int maxIter = 1000, double tol = 0.01)
        {
            K = k;
            MaxIter = maxIter;
            Tol = tol;
            Center = new double[K];
            Index = new Dictionary<int, int>();
            map_ = new Dictionary<int, int>();
            Clusters = new Dictionary<int, List<Point<T>>>();
        }

        public void Run(List<Point<T>> data)
        {
            // init center
            Index.Clear();
            map_.Clear();
            List<KeyValuePair<Point<T>, int>> sorted = data
                .Select((x, i) => new KeyValuePair<Point<T>, int>(x, i))
                .OrderBy(x => x.Key)
                .ToList();

            data = sorted.Select(x => x.Key).ToList();
            List<int> idx = sorted.Select(x => x.Value).ToList();
            for (int i = 0; i < sorted.Count; i++)
            {
                map_[i] = idx[i];
            }

            int gap = data.Count / K;
            for (int i = 0; i < K; i++)
            {
                Center[i] = data[gap * i].Value();
            }

            // iteration
            for (int i = 0; i < MaxIter; i++)
            {
                double diff = Iteration(data);
                if (diff < Tol)
                    break;
            }
        }

        private double Iteration(List<Point<T>> data)
        {
            Clusters.Clear();
            for (int i = 0; i < data.Count; i++)
            {
                double minDistance = int.MaxValue;
                int index = 0;
                for (int j = 0; j < Center.Count(); j++)
                {
                    double distance = Math.Abs(data[i].Value() - Center[j]);
                    if (distance < minDistance)
                    {
                        minDistance = distance;
                        index = j;
                    }
                }
                Index[map_[i]] = index;
                if (!Clusters.ContainsKey(index))
                {
                    Clusters[index] = new List<Point<T>>();
                }
                Clusters[index].Add(data[i]);
            }

            // update cluster
            double diff = 0;
            for (int i = 0; i < Center.Count(); i++)
            {
                if (Clusters.ContainsKey(i))
                {
                    double avg = Clusters[i].Average(p => p.Value());
                    diff += Math.Abs(Center[i] - avg);
                    Center[i] = avg;
                }
            }
            return diff;
        }
    }
}
