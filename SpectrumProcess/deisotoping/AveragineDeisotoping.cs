using SpectrumData;
using SpectrumProcess.algorithm;
using System;
using System.Collections.Generic;
using System.Linq;

namespace SpectrumProcess.deisotoping
{
    public class AveragineDeisotoping
    {
        ISearch<int> searcher;
        Averagine averagine;
        int maxExtend = 7;
        double cutoff = 0.8;

        public AveragineDeisotoping(Averagine averagine, 
            ToleranceBy by = ToleranceBy.Dalton,
            double tol = 0.1)
        {
            searcher = new BucketSearch<int>(by, tol);
            this.averagine = averagine;
        }

        protected List<List<int>> Cluster(int current,
            List<IPeak> peaks, int charge)
        {
            List<List<int>> cluster = new List<List<int>>();

            double steps = 1.0 / charge;
            IPeak peak = peaks[current];
            int index = 0;
            while (current + index < peaks.Count && index < maxExtend)
            {
                double target = peak.GetMZ() + steps * index;
                List<int> isotopics =
                    searcher.SearchContent(target)
                    .Where(index => peaks[index].GetIntensity() > 0).ToList();
                if (isotopics.Count == 0)
                    break;
                cluster.Add(isotopics);
                index++;
            }

            return cluster;
        }

        public List<IPeak> Process(List<IPeak> peaks, int maxCharge, double ion)
        {
            // init search
            List<Point<int>> points = new List<Point<int>>();
            for (int i = 0; i < peaks.Count; i++)
            {
                Point<int> point = new Point<int>(peaks[i].GetMZ(), i);
                points.Add(point);
            }
            searcher.Init(points);

            // process peaks
            List<IPeak> deisotopingPeaks = new List<IPeak>();
            HashSet<int> processed = new HashSet<int>();
            Dictionary<int, double> Restored
                = new Dictionary<int, double>();
            for (int i = 0; i < peaks.Count; i++)
            {
                if (processed.Contains(i))
                    continue;

                // find cluster
                double bestScore = 0;
                List<int> bestFitted = new List<int>();
                int bestShift = 0;
                int bestCharge = 0;
                for (int charge = 1; charge <= maxCharge; charge++)
                {
                    List<List<int>> clusters = Cluster(i, peaks, charge);
                    if (clusters.Count < 2)
                        continue;

                    // fit the score
                    Tuple<double, List<int>, int> fitted = AveragineDeisotopingHelper.Search(
                        clusters, peaks, averagine, charge, ion);
                    if (fitted.Item1 > bestScore)
                    {
                        bestScore = fitted.Item1;
                        bestFitted = fitted.Item2;
                        bestShift = fitted.Item3;
                        bestCharge = charge;
                    }
                }

                // set up new peaks
                DeisotopingPeak deisotoping = new DeisotopingPeak(peaks[i]);
                deisotopingPeaks.Add(deisotoping);

                // replace the peaks
                if (bestFitted.Count == 0 || bestScore < cutoff)
                    continue;

                deisotoping.SetMZ(peaks[i].GetMZ() + 1.0 / bestCharge * bestShift);
                deisotoping
                    .SetIntensity(bestFitted.Select(index => peaks[index].GetIntensity()).Sum());
                deisotoping.Charge = bestCharge;

                // remove the isotopic peaks
                processed.UnionWith(bestFitted);
                foreach (int matched in bestFitted)
                {
                    Restored[matched] = peaks[matched].GetIntensity();
                    peaks[matched].SetIntensity(0);
                }
            }

            // restore original peaks
            foreach (int matched in Restored.Keys)
            {
                peaks[matched].SetIntensity(Restored[matched]);
            }

            return deisotopingPeaks.OrderBy(p => p.GetMZ()).ToList();
        }


    }
}
