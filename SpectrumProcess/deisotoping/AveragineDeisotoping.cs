using SpectrumData;
using SpectrumProcess.algorithm;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpectrumProcess.deisotoping
{
    

    public class AveragineDeisotoping 
    {
        ISearch<int> searcher;
        Averagine averagine;
        double ion;
        int maxCharge;
        double cutoff = 0.8;

        public AveragineDeisotoping(double ion = 1.0078, 
            int maxCharge = 4, 
            ToleranceBy by = ToleranceBy.Dalton, 
            double tol = 0.2)
        {
            searcher = new BucketSearch<int>(by, tol);
            this.maxCharge = maxCharge;
            this.ion = ion;
            averagine = new Averagine();
        }

        protected List<List<int>> Cluster(int current, 
            List<IPeak> peaks, int charge)
        {
            List<List<int>> cluster = new List<List<int>>();

            double steps =  1.0 / charge;
            IPeak peak = peaks[current];
            int index = 0;
            while (current + index < peaks.Count)
            {
                double target = peak.GetMZ() + steps * index;
                List<int> isotopics = searcher.SearchContent(target);
                if (isotopics.Count == 0)
                    break;
                cluster.Add(isotopics);
                index++;
            }
           
            return cluster;
        }

        public List<DeisotopingPeak> Process(List<IPeak> peaks)
        {
            // init search
            List<Point<int>> points = new List<Point<int>>();
            for(int i = 0; i < peaks.Count; i++)
            {
                Point<int> point = new Point<int>(peaks[i].GetMZ(), i);
                points.Add(point);
            }
            searcher.Init(points);

            // process peaks
            List<DeisotopingPeak> deisotopingPeaks =
                new List<DeisotopingPeak>();
            HashSet<int> processed = new HashSet<int>();
            for (int i = 0; i < peaks.Count; i++)
            {
                if (processed.Contains(i))
                    continue;

                // find cluster
                double bestScore = 0;
                List<int> bestFitted = new List<int>();
                for (int charge = 1; charge <= maxCharge; charge++)
                {
                    List<List<int>> clusters = Cluster(i, peaks, charge);
                    if (clusters.Count < 2)
                        continue;

                    // fit the score
                    Tuple<double, List<int>> fitted = AveragineDeisotopingHelper.Search(
                        clusters, peaks, averagine, charge, ion);
                    if (fitted.Item1 > bestScore)
                    {
                        bestScore = fitted.Item1;
                        bestFitted = fitted.Item2;
                    }
                }

                // set up new peaks
                DeisotopingPeak deisotoping = new DeisotopingPeak(peaks[i]);
                deisotopingPeaks.Add(deisotoping);

                // replace the peaks
                if (bestFitted.Count == 0 || bestScore < cutoff)
                    continue;

                processed.UnionWith(bestFitted);
                foreach (int matched in bestFitted)
                {
                    peaks[matched].SetIntensity(0);
                }
                deisotoping
                    .SetIntensity(bestFitted.Select(index => peaks[index].GetIntensity()).Sum());
            }


            return deisotopingPeaks;
        }

        
    }
}
