using MultiGlycanTDLibrary.engine.search;
using SpectrumData;
using SpectrumProcess.algorithm;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.score
{
    public class GlycanScorerCluster : GlycanScorer, IGlycanScorer
    {
        protected ClusterKMeans<IPeak> cluster;
        protected Dictionary<string, List<double>> glycanDiagnosticPeak;
        protected ISearch<IPeak> searcher_;

        public GlycanScorerCluster(int thread = 4, double similar = 0.9,
            double binWidth = 1.0, int k = 4,
            int maxIter = 1000, double tol = 0.01 ) : base(thread, similar, binWidth)
        {
            cluster = new ClusterKMeans<IPeak>(k, maxIter, tol);
            glycanDiagnosticPeak = new Dictionary<string, List<double>>();
        }

        public GlycanScorerCluster(Dictionary<string, List<double>> diagnosticPeaks,
            ToleranceBy by, double tolerance,
            int thread = 4, double similar = 0.9,
            double binWidth = 1.0, int k = 4,
            int maxIter = 1000, double tol = 0.01) : base(thread, similar, binWidth)
        {
            cluster = new ClusterKMeans<IPeak>(k, maxIter, tol);
            glycanDiagnosticPeak = diagnosticPeaks;
            searcher_ = new BucketSearch<IPeak>(by, tolerance);
        }

        protected void ComputeCoverageScore(int scan)
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

            if (searcher_ is not null)
                searcher_.Init(points);
            foreach (SearchResult result in SpectrumResults[scan])
            {
                int nTotal = 0;
                for (int index = 0; index < peaks.Count; index++)
                {
                    int clusterIndex = cluster.Index[index];
                    if (clusterIndex == minClusterIndex)
                        continue;
                    nTotal++;
                }

                int nMatched = 0;
                foreach (int index in result.Matches.Keys)
                {
                    int clusterIndex = cluster.Index[index];
                    if (clusterIndex == minClusterIndex)
                        continue;
                    nMatched++;
                }

                if (glycanDiagnosticPeak.ContainsKey(result.Glycan))
                {
                    foreach (double mz in glycanDiagnosticPeak[result.Glycan])
                    {
                        if (searcher_.Match(mz))
                        {
                            nMatched++;
                        }
                        nTotal++;
                    }
                }

                result.Coverage = nMatched * 1.0 / nTotal;

            }
        }

        public override void AssignScore()
        {
            foreach(int scan in SpectrumResults.Keys)
            {
                if (scan == 15900)
                    System.Console.WriteLine("here");
                List<IPeak> peaks = Spectra[scan].GetPeaks();
                ComputeCoverageScore(scan);

                foreach (SearchResult result in SpectrumResults[scan])
                {
                    result.Score = GlycanScorerHelper.ComputeScore(result, peaks);
                    result.Fit = GlycanScorerHelper.ComputeFit(result, peaks);
                }
            };
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

        protected override List<SearchResult> BestResultsFromGlycan(string glycan)
        {
            List<SearchResult> BestResultsCandid = GlycanResults[glycan];

            if (glycanDiagnosticPeak.ContainsKey(glycan))
            {
                BestResultsCandid = new List<SearchResult>();

                foreach (SearchResult result in GlycanResults[glycan])
                {
                    List<Point<IPeak>> points = result.Matches
                        .Values.Select(p => new Point<IPeak>(p.Peak.GetMZ(), p.Peak)).ToList();
                    searcher_.Init(points);
                    foreach (double mz in glycanDiagnosticPeak[glycan])
                    {
                        if (searcher_.Match(mz))
                        {
                            BestResultsCandid.Add(result);
                            break;
                        }
                    }
                }

                if (BestResultsCandid.Count == 0)
                    BestResultsCandid = GlycanResults[glycan];

            }

            double bestScore = 0;
            List<SearchResult> bestResults = new List<SearchResult>();
            foreach (SearchResult result in BestResultsCandid)
            {
                if (result.Fit > bestScore)
                {
                    bestScore = result.Fit;
                    bestResults.Clear();
                }
                if (result.Fit == bestScore)
                {
                    bestResults.Add(result);
                }
            }
            return bestResults;
        }



    }
}
