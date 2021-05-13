using MultiGlycanTDLibrary.model.glycan;
using SpectrumData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.analysis
{
    public class SearchAnalyzer
    {
        int cutoff = 4;

        public SearchAnalyzer(int cutoff=4)
        {
            this.cutoff = cutoff;
        }

        public List<SearchResult> Analyze(int scan, List<IPeak> peaks,
            Dictionary<string, HashSet<int>> searchResults)
        {
            List<SearchResult> results = new List<SearchResult>();
            double best = 0;
            foreach (string glycan_id in searchResults.Keys)
            {
                double score = ComputePeakScore(peaks, searchResults[glycan_id]);
                if (score > best)
                {
                    best = score;
                    results.Clear();
                    SearchResult r = new SearchResult();
                    r.set_scan(scan);
                    r.set_glycan(glycan_id);
                    r.set_score(score);
                    results.Add(r);
                }
                else if (score == best)
                {
                    SearchResult r = new SearchResult();
                    r.set_scan(scan);
                    r.set_glycan(glycan_id);
                    r.set_score(score);
                    results.Add(r);
                }
            }

            return results;
        }

        public double ComputePeakScore(List<IPeak> peaks, 
            HashSet<int> indexList)
        {
            // if number of matches is too smaller
            if (indexList.Count < cutoff) return 0;

            return indexList.Select(index => Math.Log(peaks[index].GetIntensity())).Sum();
        }
    }
}
