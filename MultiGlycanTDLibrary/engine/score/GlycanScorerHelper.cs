using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.engine.search;
using SpectrumData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.score
{
    public class GlycanScorerHelper
    {
        public static double CosineSim(List<IPeak> p1, List<IPeak> p2, double tol = 0.1)
        {
            if (p1.Count == 0 || p2.Count == 0)
                return 0;

            double lowerBound = Math.Min(p1.Min(x => x.GetMZ()), p2.Min(x => x.GetMZ()));
            double upperBound = Math.Max(p1.Max(x => x.GetMZ()), p2.Max(x => x.GetMZ()));
            int bucketNums = (int) Math.Ceiling((upperBound - lowerBound + 1) / tol);

            double[] q1 = new double[bucketNums];
            double[] q2 = new double[bucketNums];

            foreach (IPeak pk in p1)
            {
                int index = (int)Math.Ceiling((pk.GetMZ() - lowerBound) / tol);
                q1[index] = Math.Max(q1[index], pk.GetIntensity());
            }
            foreach (IPeak pk in p2)
            {
                int index = (int)Math.Ceiling((pk.GetMZ() - lowerBound) / tol);
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


        public static double Difference(double expect, double obs, ToleranceBy by)
        {
            if (by == ToleranceBy.PPM)
            {
                return util.mass.Spectrum.To.ComputePPM(expect, obs);
            }
            return Math.Abs(expect - obs);
        }

        public static double ComputeScore(SearchResult result, double sum, ToleranceBy by, double tol)
        {

            double score = 0;
            foreach (int index in result.Matches.Keys)
            {
                PeakMatch match = result.Matches[index];
                double weight = 1 - Math.Pow(Difference(match.TheoreticMZ, match.Peak.GetMZ(), by), 4);
                //weight *= 1 - Math.Pow(match.Potentials / 1000, 4);
                score += Math.Sqrt(match.Peak.GetIntensity()) * weight;

            }
            return score / sum;
        }


    }
}
