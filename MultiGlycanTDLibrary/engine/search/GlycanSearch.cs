using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.engine.analysis;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model;
using MultiGlycanTDLibrary.model.glycan;
using SpectrumData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.search
{
    public class GlycanSearch
    {
        ISearch<string> searcher_;
        Dictionary<string, List<string>> id_map_;
        Dictionary<double, List<string>> fragments_map_;
        //private readonly int maxCharge = 2;
        private Random random;
        private readonly int lower = 1;
        private readonly int upper = 30;
        public int Seed { get; set; } = 2;

        public GlycanSearch(ISearch<string> searcher, GlycanJson glycanJson)
        {
            searcher_ = searcher;
            id_map_ = glycanJson.IDMap;
            fragments_map_ = glycanJson.Fragments;
            List<Point<string>> points = new List<Point<string>>();
            foreach (double mass in fragments_map_.Keys)
            {
                foreach (string glycan in fragments_map_[mass])
                {
                    points.Add(new Point<string>(mass, glycan));
                }
            }
            searcher_.Init(points);
            random = new Random(Seed);
        }

        public List<SearchResult> Search(List<IPeak> peaks, int precursorCharge,
            List<string> candidates, bool decoy = false, double ion = 1.0078)
        {
            // process composition
            Dictionary<string, string> glycanCandid = new Dictionary<string, string>();
            foreach (string composition in candidates)
            {
                foreach (string glycan in id_map_[composition])
                {
                    glycanCandid[glycan] = composition;
                }
            }

            // search peaks
            Dictionary<string, HashSet<int>> matched =
                new Dictionary<string, HashSet<int>>();
            for (int i = 0; i < peaks.Count; i++)
            {
                IPeak peak = peaks[i];
                double randomMass = 0;
                if (decoy)
                    randomMass = random.NextDouble() * (upper - lower) + lower;
                for (int charge = 1; charge <= precursorCharge; charge++)
                {
                    double mass = util.mass.Spectrum.To.Compute(peak.GetMZ(),
                       ion, charge);
                    if (decoy)
                        mass += randomMass;

                    List<string> glycans = searcher_.Search(mass, mass);
                    foreach (string glycan in glycans)
                    {
                        if (!glycanCandid.ContainsKey(glycan))
                            continue;

                        if (!matched.ContainsKey(glycan))
                            matched[glycan] = new HashSet<int>();
                        matched[glycan].Add((i));
                    }
                }
            }

            // compute score
            double maxScore = 0;
            List<SearchResult> results = new List<SearchResult>();
            double sum = peaks.Select(p => Math.Sqrt(p.GetIntensity())).Sum();
            foreach (string composition in candidates)
            {
                SearchResult result = new SearchResult();
                result.set_glycan(composition);

                // score the result
                double bestScore = 0;
                List<string> isomers = new List<string>();
                foreach (string glycan in matched.Keys)
                {
                    double score = matched[glycan].Select(
                        index => Math.Sqrt(peaks[index].GetIntensity())).Sum();

                    // compare score
                    if (score > bestScore)
                    {
                        bestScore = score;
                        isomers.Clear();
                        isomers.Add(glycan);
                    }
                    else if (score == bestScore)
                    {
                        isomers.Add(glycan);
                    }
                }
                result.set_isomers(isomers);
                result.set_score(bestScore / sum);

                // set up results
                if (bestScore > maxScore)
                {
                    maxScore = bestScore;
                    results.Clear();
                    results.Add(result);
                }
                else if (bestScore == maxScore)
                {
                    results.Add(result);
                }
            }

            return results;
        }

    }
}