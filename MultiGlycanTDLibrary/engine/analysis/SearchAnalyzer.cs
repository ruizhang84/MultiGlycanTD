using MultiGlycanTDLibrary.model.glycan;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.analysis
{
    public class SearchAnalyzer
    {
        public List<SearchResult> Analyze(Tuple<double, List<IGlycan>> searched,
            double mz, int scan, double retention)
        {
            Dictionary<string, SearchResult> results = new Dictionary<string, SearchResult>();
            double score = searched.Item1;
            foreach(IGlycan glycan in searched.Item2)
            {
                string name = glycan.Name();
                if (!results.ContainsKey(name))
                {
                    results[name] = new SearchResult();
                    results[name].set_mz(mz);
                    results[name].set_scan(scan);
                    results[name].set_glycan(name);
                    results[name].set_score(score);
                    results[name].set_retention(retention);
                }
                results[name].Add(glycan.ID());
            }
            return results.Values.ToList();
        }
    }
}
