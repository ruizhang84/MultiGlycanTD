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
    public class GlycanScorer
    {
        Dictionary<int, ISpectrum> Spectra;
        Dictionary<int, List<SearchResult>> SpectrumResults;
        Dictionary<string, List<SearchResult>> GlycanResults;
        Dictionary<int, SearchResult> ScoreResults;
        ToleranceBy by;
        double tol;
        double similar = 0.9;

        public GlycanScorer(Dictionary<int, ISpectrum> spectra, 
            List<SearchResult> results, ToleranceBy by, double tol)
        {
            Spectra = spectra;
            SpectrumResults = new Dictionary<int, List<SearchResult>>();
            GlycanResults = new Dictionary<string, List<SearchResult>>();
            ScoreResults = new Dictionary<int, SearchResult>();
            this.by = by;
            this.tol = tol;

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
            return ScoreResults.Select(p => p.Value).OrderBy(r => r.Scan).ToList();
        }

        public void AssignScore()
        {
            foreach (int scan in SpectrumResults.Keys)
            {
                ISpectrum spectrum = Spectra[scan];
                double sum = spectrum.GetPeaks().Select(p => Math.Sqrt(p.GetIntensity())).Sum();
                double sum2 = spectrum.GetPeaks()
                    .Select(p => p.GetIntensity() * p.GetIntensity()).Sum();
                foreach (SearchResult result in SpectrumResults[scan])
                {
                    result.Score = GlycanScorerHelper.ComputeScore(result, sum, by, tol);
                    result.Fit = GlycanScorerHelper.ComputeFit(result, sum2);
                }
            }
        }

        public void AssignSpectrumResults()
        {
            // Assign glycan to spectrum
            SpectrumResults.Clear();

            foreach (string glycan in GlycanResults.Keys)
            {
                // find best score
                double bestScore = 0;
                List<SearchResult> bestResults = new List<SearchResult>();
                foreach (SearchResult result in GlycanResults[glycan])
                {
                    // already assigned
                    if (ScoreResults.ContainsKey(result.Scan))
                        continue;

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
                        if (cosine >= similar)
                        {
                            if (!SpectrumResults.ContainsKey(result.Scan))
                            {
                                SpectrumResults[result.Scan] = new List<SearchResult>();
                            }
                            SpectrumResults[result.Scan].Add(result);
                        }
                    }
                }

            }
        }

        public void AssignGlycanResults()
        {
            // assign the result with highest score
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

                foreach (SearchResult bestResult in bestResults)
                {
                    ScoreResults[scan] = bestResult;
                    if (GlycanResults.ContainsKey(bestResult.Glycan))
                        GlycanResults.Remove(bestResult.Glycan);
                }
            }
        }

       
        


    }
}
