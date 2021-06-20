using MultiGlycanTDLibrary.model.glycan;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.analysis
{
    public class SearchMetaData
    {
        public List<SearchResult> Commit(List<SearchResult> results,
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
