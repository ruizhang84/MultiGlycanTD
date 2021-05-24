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
        ISearch<int> searcher_;
        Dictionary<string, List<string>> id_map_;
        Dictionary<string, List<double>> fragments_map_;

        public GlycanSearch(ISearch<int> searcher, GlycanJson glycanJson)
        {
            searcher_ = searcher;
            id_map_ = glycanJson.IDMap;
            fragments_map_ = glycanJson.Fragments;
        }

        void Init(List<IPeak> peaks, int precursorCharge)
        {
            List<Point<int>> points = new List<Point<int>>();
            for(int i = 0; i < peaks.Count; i++)
            {
                IPeak peak = peaks[i];
                for (int charge = 1; charge <= precursorCharge; charge++)
                {
                    double mass = util.mass.Spectrum.To.Compute(peak.GetMZ(),
                       util.mass.Spectrum.Proton, charge);
                    points.Add(new Point<int>(mass, i));
                }
            }
            searcher_.Init(points);
        }

        public List<SearchResult> Search(List<IPeak> peaks, int precursorCharge,
            List<string> candidates)
        {
            Init(peaks, precursorCharge);

            List<SearchResult> results = new List<SearchResult>();

            double maxScore = 0;
            foreach(string composition in candidates)
            {
                SearchResult result = new SearchResult();
                result.set_glycan(composition);

                // score the result
                double bestScore = 0;
                List<string> isomers = new List<string>();
                foreach (string glycan in id_map_[composition])
                {
                    HashSet<int> matchedIndex = new HashSet<int>();
                    List<double> fragments = fragments_map_[glycan];
                    foreach (double mass in fragments)
                    {
                        List<int> matched = searcher_.Search(mass, mass);
                        matchedIndex.UnionWith(matched);
                    }
                    double score = matchedIndex.Select(
                        index => Math.Log(peaks[index].GetIntensity())).Sum();

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
                result.set_score(bestScore);

                // set up results
                if (bestScore > maxScore)
                {
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
