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
        public List<SearchResult> Analyze(List<SearchResult> results,
            double mz, int scan, double retention)
        {
            foreach (SearchResult r in results)
            {
                r.set_scan(scan);
                r.set_mz(mz);
                r.set_retention(retention);
            }
            return results;
        }
    }
}
