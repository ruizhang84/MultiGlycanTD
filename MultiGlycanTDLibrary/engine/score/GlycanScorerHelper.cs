using MultiGlycanTDLibrary.engine.search;
using SpectrumData;
using SpectrumProcess.algorithm;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MultiGlycanTDLibrary.engine.score
{
    public class GlycanScorerHelper
    {
        public static double CosineSim(List<IPeak> p1, List<IPeak> p2, double binWidth = 0.1)
        {
            if (p1.Count == 0 || p2.Count == 0)
                return 0;

            double lowerBound = Math.Min(p1.Min(x => x.GetMZ()), p2.Min(x => x.GetMZ()));
            double upperBound = Math.Max(p1.Max(x => x.GetMZ()), p2.Max(x => x.GetMZ()));
            int bucketNums = (int)Math.Ceiling((upperBound - lowerBound + 1) / binWidth);

            double[] q1 = new double[bucketNums];
            double[] q2 = new double[bucketNums];

            foreach (IPeak pk in p1)
            {
                int index = (int)Math.Ceiling((pk.GetMZ() - lowerBound) / binWidth);
                q1[index] = Math.Max(q1[index], pk.GetIntensity());
            }
            foreach (IPeak pk in p2)
            {
                int index = (int)Math.Ceiling((pk.GetMZ() - lowerBound) / binWidth);
                q2[index] = Math.Max(q2[index], pk.GetIntensity());
            }

            double numerator = 0;
            double denominator1 = 0;
            double denominator2 = 0;
            for (int i = 0; i < bucketNums; i++)
            {
                numerator += q1[i] * q2[i];
                denominator1 += q1[i] * q1[i];
                denominator2 += q2[i] * q2[i];
            }

            double denominator = Math.Sqrt(denominator1) * Math.Sqrt(denominator2);
            return numerator / denominator;
        }

        public static double ComputeScore(
            SearchResult result, double sum)
        {
            double score = 0;

            foreach (int index in result.Matches.Keys)
            {
                PeakMatch match = result.Matches[index];
                score += Math.Sqrt(match.Peak.GetIntensity());
            }

            return score / sum;
        }

        public static double ComputeScore(
           SearchResult result, List<IPeak> peaks)
        {
            double score = 0;
            double sum = peaks.Select(p => Math.Sqrt(p.GetIntensity())).Sum();

            foreach (int index in result.Matches.Keys)
            {
                PeakMatch match = result.Matches[index];
                score += Math.Sqrt(match.Peak.GetIntensity());
            }
            return score / sum;
        }

        public static double ComputeFit(SearchResult result, List<IPeak> peaks)
        {
            double sum = peaks.Select(p => p.GetIntensity() * p.GetIntensity()).Sum();
            return result.Matches
                .Select(r => r.Value)
                .Select(m => m.Peak.GetIntensity() * m.Peak.GetIntensity()).Sum() / sum;
        }

        public static double ComputeFit(SearchResult result, double sum)
        {
            return result.Matches
                .Select(r => r.Value)
                .Select(m => m.Peak.GetIntensity() * m.Peak.GetIntensity()).Sum() / sum;
        }     

    }
}
