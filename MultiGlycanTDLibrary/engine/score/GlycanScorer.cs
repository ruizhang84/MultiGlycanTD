using MultiGlycanTDLibrary.engine.search;
using SpectrumData;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.score
{
    public class GlycanScorer : IGlycanScorer
    {
        protected ConcurrentDictionary<int, ISpectrum> Spectra;
        protected Dictionary<int, List<SearchResult>> SpectrumResults;
        protected Dictionary<string, List<SearchResult>> GlycanResults;
        protected Dictionary<int, List<SearchResult>> ScoreResults;

        protected double Similar = 0.9;
        protected int Thread = 4;
        protected double BinWidth = 1.0;

        public GlycanScorer(int thread = 4, double similar = 0.9, double binWidth = 1.0)
        {
            Thread = thread;
            Similar = similar;
            BinWidth = binWidth;
        }

        public void Init(ConcurrentDictionary<int, ISpectrum> spectra,
            List<SearchResult> results)
        {
            Spectra = spectra;
            SpectrumResults = new Dictionary<int, List<SearchResult>>();
            GlycanResults = new Dictionary<string, List<SearchResult>>();
            ScoreResults = new Dictionary<int, List<SearchResult>>();

            foreach (SearchResult result in results)
            {
                int scan = result.Scan;
                string glycan = result.Glycan;
                if (!SpectrumResults.ContainsKey(scan))
                {
                    SpectrumResults[scan] = new List<SearchResult>();
                }
                if (!GlycanResults.ContainsKey(glycan))
                {
                    GlycanResults[glycan] = new List<SearchResult>();
                }
                SpectrumResults[scan].Add(result);
                GlycanResults[glycan].Add(result);
            }
        }

        public void Run()
        {
            AssignScore();

            AssignSpectrumResults();
            while (SpectrumResults.Count > 0)
            {
                AssignGlycanResults();
                AssignSpectrumResults();
            }
        }

        public List<SearchResult> Result()
        {
            return ScoreResults.SelectMany(p => p.Value).OrderBy(r => r.Scan).ToList();
        }

        public virtual void AssignScore()
        {
            Parallel.ForEach(SpectrumResults.Keys,
                new ParallelOptions { MaxDegreeOfParallelism = Thread },
                scan =>
                {
                    ISpectrum spectrum = Spectra[scan];

                    double sum = spectrum.GetPeaks().Select(p => Math.Sqrt(p.GetIntensity())).Sum();
                    double sum2 = spectrum.GetPeaks()
                        .Select(p => p.GetIntensity() * p.GetIntensity()).Sum();
                    foreach (SearchResult result in SpectrumResults[scan])
                    {
                        result.Score = GlycanScorerHelper.ComputeScore(result, sum);
                        result.Fit = GlycanScorerHelper.ComputeFit(result, sum2);
                    }
                });
        }

        protected virtual List<SearchResult> BestResultsFromGlycan(string glycan)
        {
            double bestScore = 0;
            List<SearchResult> bestResults = new List<SearchResult>();
            foreach (SearchResult result in GlycanResults[glycan])
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

        protected virtual void AssignSpectrumResults()
        {
            // Assign glycan to spectrum
            SpectrumResults.Clear();

            Dictionary<string, HashSet<int>> AssignedGlycanSpectrurmResults
                = new Dictionary<string, HashSet<int>>();
            foreach (string glycan in GlycanResults.Keys)
            {
                // find best score
                List<SearchResult> bestResults =
                    BestResultsFromGlycan(glycan);
                if (bestResults.Count == 0)
                    continue;

                // record assigned results
                AssignedGlycanSpectrurmResults[glycan] = new HashSet<int>();
                foreach (SearchResult bestResult in bestResults)
                {
                    AssignedGlycanSpectrurmResults[glycan].Add(bestResult.Scan);
                }

                // assign spectrum
                foreach (SearchResult bestResult in bestResults)
                {
                    if (!SpectrumResults.ContainsKey(bestResult.Scan))
                    {
                        SpectrumResults[bestResult.Scan] = new List<SearchResult>();
                    }
                    SpectrumResults[bestResult.Scan].Add(bestResult);

                    // compute similarity and assign more
                    foreach (SearchResult result in GlycanResults[glycan])
                    {
                        if (result.Score == 0)
                            continue;
                        if (result.Scan == bestResult.Scan)
                            continue;
                        if (ScoreResults.ContainsKey(result.Scan))
                            continue;

                        double cosine = GlycanScorerHelper.CosineSim(
                            Spectra[bestResult.Scan].GetPeaks(),
                            Spectra[result.Scan].GetPeaks(), BinWidth);
                        if (cosine >= Similar)
                        {
                            if (!SpectrumResults.ContainsKey(result.Scan))
                            {
                                SpectrumResults[result.Scan] = new List<SearchResult>();
                            }
                            SpectrumResults[result.Scan].Add(result);
                            AssignedGlycanSpectrurmResults[glycan].Add(result.Scan);
                        }
                    }
                }
            }

            // Remove assigned glycan spectrum result
            Dictionary<string, List<SearchResult>> newGlycanResults
                = new Dictionary<string, List<SearchResult>>();
            foreach (string glycan in GlycanResults.Keys)
            {
                if (!AssignedGlycanSpectrurmResults.ContainsKey(glycan))
                    continue;
                HashSet<int> assigned = AssignedGlycanSpectrurmResults[glycan];
                newGlycanResults[glycan] = new List<SearchResult>();
                foreach (SearchResult result in GlycanResults[glycan])
                {
                    if (assigned.Contains(result.Scan))
                        continue;
                    newGlycanResults[glycan].Add(result);
                }
            }
            GlycanResults = newGlycanResults;
        }

        protected virtual List<SearchResult> BestResultsFromSpectrum(int scan)
        {
            double bestScore = 0;
            List<SearchResult> bestResults = new List<SearchResult>();
            foreach (SearchResult result in SpectrumResults[scan])
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

        protected virtual void AssignGlycanResults()
        {
            // assign the result with highest score in the same spectrum
            foreach (int scan in SpectrumResults.Keys)
            {

                List<SearchResult> bestResults = BestResultsFromSpectrum(scan);
                if (bestResults.Count == 0)
                    continue;

                // Remove the assigned glycans
                foreach (SearchResult bestResult in bestResults)
                {
                    if (!ScoreResults.ContainsKey(scan))
                        ScoreResults[scan] = new List<SearchResult>();
                    else
                    {
                        // make score always highest
                        double score = ScoreResults[scan].First().Score;
                        if (score < bestResult.Score)
                            ScoreResults[scan].Clear();
                    }
                    ScoreResults[scan].Add(bestResult);

                    if (GlycanResults.ContainsKey(bestResult.Glycan))
                    {
                        GlycanResults[bestResult.Glycan]
                            = GlycanResults[bestResult.Glycan].Where(
                                r => r.Charge != bestResult.Charge).ToList();
                        if (GlycanResults[bestResult.Glycan].Count == 0)
                            GlycanResults.Remove(bestResult.Glycan);
                    }
                }
            }
        }





    }
}
