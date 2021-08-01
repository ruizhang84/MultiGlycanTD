using MultiGlycanTDLibrary.engine.search;
using SpectrumData;
using SpectrumProcess.algorithm;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.score
{
    public class GlycanScorerClusterShare : GlycanScorerCluster, IGlycanScorer
    {
        public GlycanScorerClusterShare(int thread = 4, double similar = 0.9,
            double binWidth = 1.0, int k = 4,
            int maxIter = 1000, double tol = 0.01) : 
            base(thread, similar, binWidth, k, maxIter, tol)
        {
        }
        public GlycanScorerClusterShare(Dictionary<string, List<double>> diagnosticPeaks,
           ToleranceBy by, double tolerance,
           int thread = 4, double similar = 0.9,
           double binWidth = 1.0, int k = 4,
           int maxIter = 1000, double tol = 0.01) : 
            base(diagnosticPeaks, by, tolerance, thread, similar, binWidth,
                k, maxIter, tol)
        {
        }

        public override void AssignScore()
        {

            Dictionary<double, List<double>> peakRatioList = 
                new Dictionary<double, List<double>>();
            foreach (int scan in SpectrumResults.Keys)
            {
                List<IPeak> peaks = Spectra[scan].GetPeaks();
                ComputeCoverageScore(scan);

                double sum = peaks.Select(p => Math.Sqrt(p.GetIntensity())).Sum();
                Dictionary<double, double> ratios = new Dictionary<double, double>();
                foreach (SearchResult result in SpectrumResults[scan])
                {
                    // finding common peaks and intensity ratio;
                    foreach(PeakMatch match in result.Matches.Values)
                    {
                        double mz = Math.Round(match.TheoreticMZ, 2);
                        double ratio = Math.Sqrt(match.Peak.GetIntensity()) / sum;
                        if (!ratios.ContainsKey(mz))
                        {
                            ratios[mz] = ratio;
                        }
                        else
                        {
                            ratios[mz] = Math.Max(ratio, ratios[mz]);
                        }
                    }

                    // compute fit scores
                    result.Fit = GlycanScorerHelper.ComputeFit(result, peaks);
                }

                foreach(double mz in ratios.Keys)
                {
                    if (!peakRatioList.ContainsKey(mz))
                    {
                        peakRatioList[mz] = new List<double>();
                    }
                    peakRatioList[mz].Add(ratios[mz]);
                }

            };

            // Averaging the peak intensity ratios
            Dictionary<double, double> peakRatios =
                new Dictionary<double, double>();

            foreach(double mz in peakRatioList.Keys)
            {
                // at least find two spectrum
                if (peakRatioList[mz].Count > 1)
                    peakRatios[mz] = peakRatioList[mz].Average();
            }

            // compute score
            foreach (int scan in SpectrumResults.Keys)
            {
                List<IPeak> peaks = Spectra[scan].GetPeaks();
                double sum = peaks.Select(p => Math.Sqrt(p.GetIntensity())).Sum();
                foreach (SearchResult result in SpectrumResults[scan])
                {
                    result.Score = 0;

                    List<double> scores = new List<double>();
                    double score = GlycanScorerHelper.ComputeScore(result, peaks);
                    foreach(PeakMatch match in result.Matches.Values)
                    {
                        double mz = Math.Round(match.TheoreticMZ, 2);
                        double ratio = Math.Sqrt(match.Peak.GetIntensity()) / sum;
                        if (peakRatios.ContainsKey(mz))
                        {
                            scores.Add(peakRatios[mz] / ratio * score);
                        }
                    }
                    if (scores.Count > 0)
                        result.Score = scores.Average();
                }
            }

        }

    }
}
