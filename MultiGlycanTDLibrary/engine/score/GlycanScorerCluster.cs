using MultiGlycanTDLibrary.engine.search;
using SpectrumData;
using SpectrumProcess.algorithm;
using System.Collections.Generic;
using System.Linq;

namespace MultiGlycanTDLibrary.engine.score
{
    public class GlycanScorerCluster : GlycanScorer, IGlycanScorer
    {
        ClusterKMeans<IPeak> cluster;
        public GlycanScorerCluster(int thread = 4, double similar = 0.9,
            double binWidth = 1.0, int k = 4,
            int maxIter = 1000, double tol = 0.01) : base(thread, similar, binWidth)
        {
            cluster = new ClusterKMeans<IPeak>(k, maxIter, tol);
        }

        protected override void AssignScore()
        {
            foreach (int scan in SpectrumResults.Keys)
            {
                List<IPeak> peaks = Spectra[scan].GetPeaks();
                List<Point<IPeak>> points =
                    peaks.Select(p => new Point<IPeak>(p.GetIntensity(), p)).ToList();
                cluster.Run(points);
                double minClusterIntensity = int.MaxValue;
                int minClusterIndex = 0;
                foreach (int index in cluster.Clusters.Keys)
                {
                    double average =
                        cluster.Clusters[index].Average(peaks => peaks.Content().GetIntensity());
                    if (average < minClusterIntensity)
                    {
                        minClusterIntensity = average;
                        minClusterIndex = index;
                    }
                }

                foreach (SearchResult result in SpectrumResults[scan])
                {
                    int nTotal = result.Matches.Count;
                    int nMatched = 0;
                    foreach (int index in result.Matches.Keys)
                    {
                        int clusterIndex = cluster.Index[index];
                        if (clusterIndex == minClusterIndex)
                            continue;
                        nMatched++;
                    }
                    result.Coverage = nMatched * 1.0 / nTotal;
                }

                foreach (SearchResult result in SpectrumResults[scan])
                {

                    result.Score = GlycanScorerHelper.ComputeScore(result, peaks);
                    result.Fit = GlycanScorerHelper.ComputeFit(result, peaks);
                }
            }
        }

        protected override List<SearchResult> BestResultsFromSpectrum(int scan)
        {
            double bestScore = 0;
            List<SearchResult> tempBestResults = new List<SearchResult>();
            foreach (SearchResult result in SpectrumResults[scan])
            {
                if (result.Coverage > bestScore)
                {
                    bestScore = result.Coverage;
                    tempBestResults.Clear();
                }
                if (result.Coverage == bestScore)
                {
                    tempBestResults.Add(result);
                }

            }

            bestScore = 0;
            List<SearchResult> bestResults = new List<SearchResult>();
            foreach (SearchResult result in tempBestResults)
            {
                if (result.Score > bestScore)
                {
                    bestScore = result.Score;
                    bestResults.Clear();
                }
                if (result.Score == bestScore)
                {
                    bestResults.Add(result);
                }

            }

            return bestResults;
        }

    }
}
