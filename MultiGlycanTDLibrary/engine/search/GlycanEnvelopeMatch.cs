using MultiGlycanTDLibrary.model;
using SpectrumData;
using SpectrumProcess.algorithm;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.search
{
    public class GlycanEnvelopeMatch
    {
        EnvelopeProcessor processor;
        Dictionary<string, List<double>> distr_map;
        Dictionary<string, List<double>> mass_map;
        double cutoff_ = 0.8;
        protected double tolerance_;
        protected ToleranceBy type_;

        public GlycanEnvelopeMatch(EnvelopeProcessor processor,
            CompdJson compdJson, ToleranceBy by, double tol, double cutoff = 0.8)
        {
            this.processor = processor;
            distr_map = compdJson.DistrMap;
            mass_map = compdJson.MassMap;
            cutoff_ = cutoff;
            type_ = by;
            tolerance_ = tol;
        }

        public List<SortedDictionary<int, IPeak>> Combinator(
            SortedDictionary<int, List<IPeak>> cluster)
        {
            List<SortedDictionary<int, IPeak>> results = 
                new List<SortedDictionary<int, IPeak>>();
            Queue<Tuple<int, SortedDictionary<int, IPeak>>> queue =
                new Queue<Tuple<int, SortedDictionary<int, IPeak>>>();
            int start = cluster.Keys.Min();
            int end = cluster.Keys.Max();
            foreach (IPeak peak in cluster[start])
            {
                SortedDictionary<int, IPeak> dict = new SortedDictionary<int, IPeak>();
                dict[start] = peak;
                int next = start + 1;
                while (!cluster.ContainsKey(next) && next < end)
                {
                    next++;
                }
                queue.Enqueue(Tuple.Create(next, dict));
            }

            while (queue.Count > 0)
            {
                Tuple<int, SortedDictionary<int, IPeak>> tuple = queue.Dequeue();
                int curr = tuple.Item1;
                SortedDictionary<int, IPeak> dict = tuple.Item2;
                if (curr > end)
                {
                    results.Add(dict);
                    continue;
                }

                foreach (IPeak peak in cluster[curr])
                {
                    SortedDictionary<int, IPeak> newDict =
                        new SortedDictionary<int, IPeak>(dict);

                    newDict[curr] = peak;
                    int next = curr + 1;
                    while (!cluster.ContainsKey(next) && next < end)
                    {
                        next++;
                    }
                    queue.Enqueue(Tuple.Create(next, newDict));
                }


            }
            return results;
        }

        private double Score(List<double> alignedDistr, List<double> alignedIntensity)
        {
            // compute correlation
            double distrMean = alignedDistr.Average();
            double intensityMean = alignedIntensity.Average();

            double norminator = 0.0;
            for (int i = 0; i < alignedDistr.Count; i++)
            {
                norminator += (alignedDistr[i] - distrMean) * (alignedIntensity[i] - intensityMean);
            }
            double denominator1 = 0.0;
            for (int i = 0; i < alignedDistr.Count; i++)
            {
                denominator1 += (alignedDistr[i] - distrMean) * (alignedDistr[i] - distrMean);
            }
            double denominator2 = 0.0;
            for (int i = 0; i < alignedDistr.Count; i++)
            {
                denominator2 += (alignedIntensity[i] - intensityMean)
                    * (alignedIntensity[i] - intensityMean);
            }
            return norminator / Math.Sqrt(denominator1 * denominator2);
        }

        private void AlignData(List<double> distr, SortedDictionary<int, IPeak> isotopics,
            int distrIndex, int isotopicIndex,
            out List<double> alignedDistr, out List<IPeak> alignedPeaks)
        {
            alignedDistr = new List<double>();
            alignedPeaks = new List<IPeak>();

            int alignedDistrIndex = distrIndex;
            int alignedIsotopicIndex = isotopicIndex;
            int start = isotopics.Keys.Min();
            int end = isotopics.Keys.Max();

            // align the data
            while (alignedDistrIndex >= 0 && alignedIsotopicIndex >= start)
            {
                if (isotopics.ContainsKey(alignedIsotopicIndex))
                {
                    alignedDistr.Add(distr[alignedDistrIndex]);
                    alignedPeaks.Add(isotopics[alignedIsotopicIndex]);
                }
                alignedDistrIndex--;
                alignedIsotopicIndex--;
            }

            alignedDistr.Reverse();
            alignedPeaks.Reverse();

            alignedDistrIndex = distrIndex + 1;
            alignedIsotopicIndex = isotopicIndex + 1;

            while (alignedDistrIndex < distr.Count && alignedIsotopicIndex <= end)
            {
                if (isotopics.ContainsKey(alignedIsotopicIndex))
                {
                    alignedDistr.Add(distr[alignedDistrIndex]);
                    alignedPeaks.Add(isotopics[alignedIsotopicIndex]);
                }
                alignedDistrIndex++;
                alignedIsotopicIndex++;
            }
        }

        public double Fit(List<double> distr, int index, SortedDictionary<int, IPeak> isotopics)
        {
            // find the peaks intensity and align data
            List<double> alignedDistr;
            List<IPeak> alignedPeaks;
            AlignData(distr, isotopics, index, 0,
                out alignedDistr, out alignedPeaks);

            // compute correlation
            if (alignedPeaks.Count == 0)
                return 0;
            return Score(alignedDistr, alignedPeaks.Select(p => p.GetIntensity()).ToList());
        }

        private bool IsMatch(double expect, double observe)
        {
            if (type_ == ToleranceBy.PPM)
            {
                return Math.Abs(expect - observe) / expect * 1000000.0 < tolerance_;
            }
            return Math.Abs(expect - observe) < tolerance_;
        }

        private int SearchIndex(List<double> data, double expect)
        {
            int start = 0;
            int end = data.Count - 1;

            while (start <= end)
            {
                int mid = end + (start - end) / 2;
                if (IsMatch(expect, data[mid]))
                {
                    return mid;
                }
                else if (data[mid] > expect)
                {
                    end = mid - 1;
                }
                else
                {
                    start = mid + 1;
                }
            }

            return -1;
        }

        public List<string> Match(List<string> precursorMatched,
            List<IPeak> peaks, double mz, int charge, double ion = 1.0078)
        {
            processor.Init(peaks);
            SortedDictionary<int, List<IPeak>> cluster = processor.Cluster(mz, charge);
            if (cluster.Count == 0)
                return new List<string>();

            double mass = util.mass.Spectrum.To.Compute(mz, ion, charge);

            List<string> results = new List<string>();
            List<SortedDictionary<int, IPeak>> clustered = Combinator(cluster);
            foreach (string compose in precursorMatched)
            {
                List<double> distr = distr_map[compose];
                List<double> massList = mass_map[compose];
                int index = SearchIndex(massList, mass);

                // score by fitting the distr
                double bestScore = -1;
                foreach (SortedDictionary<int, IPeak> sequence in clustered)
                {
                    double score = Fit(distr, index, sequence);
                    if (score > bestScore)
                    {
                        bestScore = score;
                    }
                }

                // score
                if (bestScore < cutoff_)
                    continue;

                results.Add(compose);
            }

            return results;
        }

    }
}
