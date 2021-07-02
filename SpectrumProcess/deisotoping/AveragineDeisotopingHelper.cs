using SpectrumData;
using SpectrumProcess.brain;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpectrumProcess.deisotoping
{
    public class AveragineDeisotopingHelper
    {
        public static Tuple<double, List<int>, int> Search(List<List<int>> cluster, List<IPeak> peaks,
            Averagine averagine, int charge, double ion = 1.0078)
        {
            List<List<int>> clustered = Combinator(cluster);

            double maxScore = 0;
            int bestShift = 0;
            List<int> best = new List<int>();
            foreach (List<int> sequence in clustered)
            {
                Tuple<int, double> shiftScore = Fit(peaks, sequence, averagine, charge, ion);
                if (shiftScore.Item2 >= maxScore)
                {
                    maxScore = shiftScore.Item2;
                    bestShift = shiftScore.Item1;
                    best = sequence;
                }
            }
            return Tuple.Create(maxScore, best, bestShift);
        }

        // Create a isotopic cluster: list of peak lists -> list of peak sequences
        public static List<List<int>> Combinator(List<List<int>> cluster)
        {
            
            List<List<int>> results = new List<List<int>>();
            if (cluster.Count == 0)
                return results;

            Queue<List<int>> queue = new Queue<List<int>>();
            foreach (int peak in cluster[0])
            {
                List<int> seq = new List<int>();
                seq.Add(peak);
                queue.Enqueue(seq);
            }

            while (queue.Count > 0)
            {
                List<int> seq = queue.Dequeue();
                int curr = seq.Count;

                if (curr == cluster.Count)
                {
                    results.Add(seq);
                    continue;
                }

                foreach (int peak in cluster[curr])
                {
                    List<int> newSeq = new List<int>(seq);
                    newSeq.Add(peak);
                    queue.Enqueue(newSeq);
                }

            }
            return results;
        }

        // Compute a fit score
        public static double Score(List<double> alignedDistr, List<double> alignedIntensity)
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

        public static void Align(List<double> distr, List<double> isotopicPeaks,
            int distrIndex, int isotopicIndex,
            List<double> alignedDistr, List<double> alignedIntensity)
        {
            int alignedDistrIndex = distrIndex;
            int alignedIsotopicIndex = isotopicIndex;

            // align the data
            while (alignedDistrIndex >= 0 && alignedIsotopicIndex >= 0)
            {
                alignedDistr.Add(distr[alignedDistrIndex]);
                alignedIntensity.Add(isotopicPeaks[alignedIsotopicIndex]);
                alignedDistrIndex--;
                alignedIsotopicIndex--;
            }

            alignedDistr.Reverse();
            alignedIntensity.Reverse();

            alignedDistrIndex = distrIndex + 1;
            alignedIsotopicIndex = isotopicIndex + 1;

            while (alignedDistrIndex < distr.Count && alignedIsotopicIndex <= isotopicPeaks.Count)
            {

                alignedDistr.Add(distr[alignedDistrIndex]);
                alignedIntensity.Add(isotopicPeaks[alignedIsotopicIndex]);
                
                alignedDistrIndex++;
                alignedIsotopicIndex++;
            }
        }

        // Fit the isotopic cluster
        public static Tuple<int, double> Fit(List<IPeak> peaks, List<int> isotopics,
            Averagine averagine, int charge, double ion = 1.0078)
        {
            // find the average mass by most abundant peak
            double maxIntensity = isotopics.Max(index => peaks[index].GetIntensity());
            double mz = isotopics.Where(index => peaks[index].GetIntensity() == maxIntensity)
                .Average(index => peaks[index].GetMZ());

            // find distribution
            double mass = (mz - ion) * charge;
            Compound compound = averagine.Fit(mass);
            List<double> distr = Brain.Run.Distribute(compound, isotopics.Count);

            // find the max distribution
            int maxDistrIndex = 0;
            double maxProb = 0;
            for (int i = 0; i < distr.Count; i++)
            {
                if (distr[i] > maxProb)
                {
                    maxProb = distr[i];
                    maxDistrIndex = i;
                }
            }

            // find the max intensity peak
            int maxValue = isotopics.OrderByDescending(index => peaks[index].GetIntensity()).First();
            int maxIndex = isotopics.IndexOf(maxValue);


            // corner case
            if (maxIndex > maxDistrIndex)
                return Tuple.Create(0, 0.0);

            //align data
            List<double> alignedDistr = new List<double>();
            List<double> alignedIntensity = new List<double>();
            Align(distr, isotopics.Select(index => peaks[index].GetIntensity()).ToList(),
                maxDistrIndex, maxIndex,
                alignedDistr, alignedIntensity);

            // compute correlation
            int shift = maxIndex - maxDistrIndex;
            return Tuple.Create(shift, Score(alignedDistr, alignedIntensity)); 
        }



    }
}
