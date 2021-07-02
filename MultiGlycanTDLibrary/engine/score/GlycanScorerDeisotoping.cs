using MultiGlycanTDLibrary.engine.search;
using SpectrumData;
using SpectrumProcess.algorithm;
using SpectrumProcess.deisotoping;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.score
{
    public class GlycanScorerDeisotoping : IGlycanScorer
    {
        ConcurrentDictionary<int, ISpectrum> Spectra;
        Dictionary<int, List<SearchResult>> SpectrumResults;
        Dictionary<string, List<SearchResult>> GlycanResults;
        Dictionary<int, List<SearchResult>> ScoreResults;
        Averagine Averagine;
        int MaxCharge;
        ToleranceBy By;
        double Tol;

        double Similar = 0.9;
        int Thread = 4;

        public GlycanScorerDeisotoping(ConcurrentDictionary<int, ISpectrum> spectra,
             Averagine averagine, int maxCharge, ToleranceBy by, double tol, 
             int thread = 4, double similar = 0.9)
        {
            Spectra = spectra;
            Averagine = averagine;
            MaxCharge = maxCharge;
            By = by;
            Tol = tol;
            Thread = thread;
            Similar = similar;
        }

        public void Init(List<SearchResult> results)
        {
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

            AssignGlycanResults();
            while (GlycanResults.Count > 0)
            {
                AssignSpectrumResults();
                AssignGlycanResults();
            }
        }

        public List<SearchResult> Result()
        {
            return ScoreResults.SelectMany(p => p.Value).OrderBy(r => r.Scan).ToList();
        }

        protected void LocalAssignScore(int scan, 
            AveragineDeisotoping deisotoping)
        {
            foreach (SearchResult result in SpectrumResults[scan])
            {
                List<IPeak> peaks = deisotoping
                    .Process(Spectra[scan].GetPeaks(), result.Ion)
                    .Select(p => p as IPeak).ToList();
                result.Score = GlycanScorerHelper.ComputeScore(result, peaks);
                result.Fit = GlycanScorerHelper.ComputeFit(result, peaks);
            }
        }

        protected void AssignScoreTask(ConcurrentQueue<int> ScanQueue)
        {
            AveragineDeisotoping deisotoping =
                new AveragineDeisotoping(Averagine, MaxCharge, By, Tol);
            while (ScanQueue.TryDequeue(out int scan))
            {
                LocalAssignScore(scan, deisotoping);
            }
        }

        public void AssignScore()
        {
            ConcurrentQueue<int> ScanQueue =
                new ConcurrentQueue<int>(SpectrumResults.Keys);

            List<Task> scoer = new List<Task>();
            for (int i = 0; i < Thread; i++)
            {
                Task LastTask = new Task(() => AssignScoreTask(ScanQueue));
                LastTask.Start();
                scoer.Add(LastTask);
            }
            Task.WaitAll(scoer.ToArray());
        }

        public void AssignSpectrumResults()
        {
            // Assign glycan to spectrum
            SpectrumResults.Clear();

            Dictionary<string, HashSet<int>> AssignedGlycanSpectrurmResults
                = new Dictionary<string, HashSet<int>>();
            foreach (string glycan in GlycanResults.Keys)
            {
                // find best score
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
                            Spectra[result.Scan].GetPeaks());
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

        public void AssignGlycanResults()
        {
            // assign the result with highest score in the same spectrum
            foreach (int scan in SpectrumResults.Keys)
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
