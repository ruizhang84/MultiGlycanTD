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
    public class GlycanSearchForTest
    {
        ISearch<string> searcher_;
        ISearch<string> searcher2_;
        Dictionary<string, List<string>> id_map_;
        Dictionary<double, List<string>> fragments_map_;
        private readonly int maxCharge = 2;
        private Random random;
        private readonly int lower = 1;
        private readonly int upper = 30;
        public int Seed { get; set; } = 2;

        public GlycanSearchForTest(ISearch<string> searcher, ISearch<string> searcher2, 
            GlycanJson glycanJson, GlycanJson glycanJson2)
        {
            searcher_ = searcher;
            searcher2_ = searcher2;
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

            fragments_map_ = glycanJson2.Fragments;
            List<Point<string>> points2 = new List<Point<string>>();
            foreach (double mass in fragments_map_.Keys)
            {
                foreach (string glycan in fragments_map_[mass])
                {
                    points2.Add(new Point<string>(mass, glycan));
                }
            }
            searcher2_.Init(points2);
            random = new Random(Seed);
        }

        public List<SearchResult> Search(List<IPeak> peaks, int precursorCharge,
            List<string> candidates, bool decoy = false, int testCharge=-1)
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
            Dictionary<string, HashSet<int>> matched2 =
                new Dictionary<string, HashSet<int>>();
            for (int i = 0; i < peaks.Count; i++)
            {
                IPeak peak = peaks[i];
                double randomMass = 0;
                if (decoy)
                    randomMass = random.NextDouble() * (upper - lower) + lower;
                for (int charge = 1; charge <= maxCharge; charge++)
                {
                    double mass = util.mass.Spectrum.To.Compute(peak.GetMZ(),
                       util.mass.Spectrum.Proton, charge);
                    if (decoy)
                        mass += randomMass;

                    List<string> glycans = searcher_.Search(mass, mass);
                    List<string> glycans2 = searcher2_.Search(mass, mass);
                    if (testCharge > 0)
                    {
                        if (charge != testCharge)
                        {
                            glycans2.Clear();
                        }
                    }
                    foreach (string glycan in glycans)
                    {
                        if (!glycanCandid.ContainsKey(glycan))
                            continue;

                        if (!matched.ContainsKey(glycan))
                            matched[glycan] = new HashSet<int>();
                        matched[glycan].Add((i));
                    }
                    foreach (string glycan in glycans2)
                    {
                        if (!glycanCandid.ContainsKey(glycan))
                            continue;

                        if (!matched2.ContainsKey(glycan))
                            matched2[glycan] = new HashSet<int>();
                        matched2[glycan].Add((i));
                    }
                }
            }

            // compute score
            double maxScore = 0;
            List<SearchResult> results = new List<SearchResult>();
            double sum = peaks.Select(p => Math.Log(p.GetIntensity())).Sum();
            foreach(string composition in candidates)
            {
                SearchResult result = new SearchResult();
                result.set_glycan(composition);

                // score the result
                double bestScore = 0;
                double count = 0;
                List<string> isomers = new List<string>();
                foreach (string glycan in matched.Keys)
                {
                    double score = matched[glycan].Select(
                        index => Math.Log(peaks[index].GetIntensity())).Sum();
                    Console.WriteLine(glycan + " " + score.ToString() + ": peak, " +
                        string.Join("-", matched[glycan].Select(
                        index => peaks[index].GetMZ().ToString() + " " + peaks[index].GetIntensity().ToString())));

                    // compare score
                    if (score > bestScore)
                    {
                        bestScore = score;
                        if (matched2.ContainsKey(glycan))
                            count = matched2[glycan].Count;
                        isomers.Clear();
                        isomers.Add(glycan);
                    }
                    else if (score == bestScore)
                    {
                        if (matched2.ContainsKey(glycan))
                            count = Math.Max(count, matched2[glycan].Count);
                        isomers.Add(glycan);
                    }
                }
                result.set_isomers(isomers);
                result.set_score(count);

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
